// Collecting, dealing with, calculating systematic uncertainties:
//   qlcms
//	 rhofitmax
//	 nevt_avg
//	 mT choice
// Then, finally plotting them with result points on m_T vs Lévy parameter and sqrt(sNN) plots
// Do not forget that the point markers here should be almost the same size as the lines connecting them; the uncertainty bars should be filled areas instead of normal error bars

#include "header_for_all_emissionsource.h"
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPad.h>
#include <TColor.h>
#include <TSystem.h>
#include <vector>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <TF1.h>
#include <TLatex.h>

bool debug = true; // set to true to print parameter results for each systematic variation
bool plot_mtavg = false; // set to true to plot m_T-averaged param vs sNN as well
bool mtchoice_syst = false;
bool use_larger_same_side_shift = false; // if true, use max(|diff1|,|diff2|) when both variations are on the same side of default
bool comparewithSTAR = true; // set to true to include analytic alpha vs sNN function fitted to STAR data

const int NIKT_SNN_PLOTS = 2;
int ikt_indices_to_plot[NIKT_SNN_PLOTS] = {0, 4}; // choose which m_T bins are shown on parameter-vs-sNN plots

const char* levy_params[3] = {"alpha","R","N"};
const char* levy_params_3d[3] = {"R_out","R_side","R_long"};

// To be set
//int NEVT_AVG_DEFAULT = 5; // single index

void correct_syserr_direction(double* param_uncert_up, double* param_uncert_dn, 
                              TGraphAsymmErrors* systcheckgraph1, TGraphAsymmErrors* systcheckgraph2,
                              TGraphAsymmErrors* defaultgraph, int ikt)
{
    // This function will check the direction of the deviation of the systcheckgraph from the defaultgraph.
    // If the systcheck point is above the default, it will add the squared difference to the upward uncertainty; if below, to the downward uncertainty.
    // If both systcheck points are on the same side, then add the average of the two differences.
    // Sample:
    // func(*param*_syserr_up[ikt], *param*_syserr_dn[ikt], *param*_*systcheck1[ienergy][1], *param*_*systcheck*[ienergy][2], *param*_default[ienergy])
    double def = defaultgraph->GetY()[ikt];
    double y1 = systcheckgraph1->GetY()[ikt];
    double y2 = systcheckgraph2->GetY()[ikt];
    
    // Taking care of potential zero or near-zero default values to avoid misleading large percent differences; if default is zero, we cannot determine the direction, so we skip this point
    if(def < 1e-12)
    {
        // skip tiny/invalid points
        cerr << "Zero value in default graph! ikt: " << ikt << endl;
        return;
    }

    bool valid1 = std::isfinite(y1) && (y1 >= 1e-12);
    bool valid2 = std::isfinite(y2) && (y2 >= 1e-12);
    if(!valid1) cerr << "Invalid/tiny value in systcheckgraph1! Skipping this variation. ikt: " << ikt << endl;
    if(!valid2) cerr << "Invalid/tiny value in systcheckgraph2! Skipping this variation. ikt: " << ikt << endl;
    if(!valid1 && !valid2)
    {
        return;
    }

    if(valid1 && !valid2)
    {
        double diff1 = y1 - def;
        if(diff1 >= 0.) param_uncert_up[ikt] += pow(diff1, 2);
        else param_uncert_dn[ikt] += pow(diff1, 2);
        return;
    }
    if(!valid1 && valid2)
    {
        double diff2 = y2 - def;
        if(diff2 >= 0.) param_uncert_up[ikt] += pow(diff2, 2);
        else param_uncert_dn[ikt] += pow(diff2, 2);
        return;
    }

    // Both valid from here
    double diff1 = y1 - def;
    double diff2 = y2 - def;
    double same_side_shift = 0.5 * (diff1 + diff2); // default behavior: average
    if(use_larger_same_side_shift)
    {
        same_side_shift = (fabs(diff1) >= fabs(diff2)) ? diff1 : diff2;
    }

    if((diff1 >= 0. && diff2 >= 0.) || (diff1 < 0. && diff2 < 0.))
    {
        // both differences on the same side -> use average or larger-shift mode
        if(same_side_shift >= 0.)
        {
            param_uncert_up[ikt] += pow(same_side_shift, 2);
        }
        else
        {
            param_uncert_dn[ikt] += pow(same_side_shift, 2);
        }
    }
    else
    {
        // differences on opposite sides -> use each difference separately
        if(diff1 >= 0.)
        {
            param_uncert_up[ikt] += pow(diff1, 2);
        }
        else
        {
            param_uncert_dn[ikt] += pow(diff1, 2);
        }
        if(diff2 >= 0.)
        {
            param_uncert_up[ikt] += pow(diff2, 2);
        }
        else
        {
            param_uncert_dn[ikt] += pow(diff2, 2);
        }
    }
    //return diffavg; // change from void to double if you want to return the average difference
}

// I do not know if this makes any sense, this will return the average syst. difference in percent
double calc_syst_contribution(TGraphAsymmErrors* defaultgraph, double* errups, double* errdns, bool up=true)
{
    if(!defaultgraph) return 0.;
    double percent = 0.;
    int points = 0;
    for(int ikt=0; ikt<NKT; ikt++)
    {
        double defy = defaultgraph->GetY()[ikt];
        if(defy < 1e-12) continue; // skip tiny/invalid points
        if(up)
        {
            // errups is an absolute upward uncertainty (magnitude)
            double mag = errups ? errups[ikt] : 0.;
            percent += mag / defy;
            points++;
        }
        else
        {
            // errdns is an absolute downward uncertainty (magnitude)
            double mag = errdns ? errdns[ikt] : 0.;
            percent += mag / defy;
            points++;
        }
    }
    if(points == 0) return 0.;
    percent /= points;
    percent *= 100.; // to percent
    return percent;
}

// Return median percent contribution (same inputs as calc_syst_contribution)
double calc_syst_contrib_median(TGraphAsymmErrors* defaultgraph, double* errups, double* errdns, bool up=true)
{
    if(!defaultgraph) return 0.;
    std::vector<double> vals;
    vals.reserve(NKT);
    for(int ikt=0; ikt<NKT; ikt++){
        double defy = defaultgraph->GetY()[ikt];
        if(defy < 1e-12) continue;
        double mag = 0.;
        if(up){ mag = errups ? errups[ikt] : 0.; }
        else { mag = errdns ? errdns[ikt] : 0.; }
        // store percent for this bin
        vals.push_back( (mag / defy) * 100. );
    }
    if(vals.empty()) return 0.;
    std::sort(vals.begin(), vals.end());
    size_t n = vals.size();
    if(n % 2 == 1) return vals[n/2];
    // even -> average two middle
    return 0.5 * (vals[n/2 - 1] + vals[n/2]);
}

Double_t starfunc(double x)
{
    return 0.85+pow(x,-0.14);
}


// -------- MAIN FUNCTION --------
int calc_and_plot_syserr(int energy_to_plot=-1)
{
    // Get m_T bin centers
    double mtbin_centers[NKT];
    for(int ii=0; ii<NKT; ii++)
    {
        //mtbin_centers[ii] = 0.5 * (ktbins[ii] + ktbins[ii+1]);
        mtbin_centers[ii] = kT_center[ii];
        mtbin_centers[ii] = sqrt(mtbin_centers[ii]*mtbin_centers[ii] + Mass2_pi);
    }
    // x-error arrays (half-bin widths) - set to zero for now
    double xerr_low[NKT] = {0};
    double xerr_high[NKT] = {0};

    auto make_graph_if_populated = [&](double* yvals, double* yerrdn, double* yerrup) -> TGraphAsymmErrors*
    {
        bool has_nonzero = false;
        for(int ikt=0; ikt<NKT; ikt++) {
            if(std::isfinite(yvals[ikt]) && std::fabs(yvals[ikt]) > 1e-12) {
                has_nonzero = true;
                break;
            }
        }
        if(!has_nonzero) return nullptr;
        return new TGraphAsymmErrors(NKT, mtbin_centers, yvals, xerr_low, xerr_high, yerrdn, yerrup);
    };

    // Build m_T-choice uncertainty from neighboring bins around a chosen center bin.
    auto calc_neighbor_uncert = [](TGraphAsymmErrors* g, int center_idx, double& up, double& dn)
    {
        up = 0.;
        dn = 0.;
        if(!g) return;
        if(center_idx < 0 || center_idx >= g->GetN()) return;
        double center = g->GetY()[center_idx];
        if(!std::isfinite(center)) return;
        double up_sq = 0., dn_sq = 0.;
        for(int j : {center_idx-1, center_idx+1})
        {
            if(j < 0 || j >= g->GetN()) continue;
            double nb = g->GetY()[j];
            if(!std::isfinite(nb)) continue;
            double diff = nb - center;
            if(diff >= 0.) up_sq += diff*diff;
            else dn_sq += diff*diff;
        }
        up = sqrt(up_sq);
        dn = sqrt(dn_sq);
    };

    auto fill_3d_radii_from_file = [](TFile* f, int ikt,
                                      double& Rout, double& dRout,
                                      double& Rside, double& dRside,
                                      double& Rlong, double& dRlong) -> bool
    {
        if(!f) return false;
        TH1F* hRout = (TH1F*)f->Get(Form("R_out_hist_ikt%i", ikt));
        TH1F* hRside = (TH1F*)f->Get(Form("R_side_hist_ikt%i", ikt));
        TH1F* hRlong = (TH1F*)f->Get(Form("R_long_hist_ikt%i", ikt));
        if(!hRout || !hRside || !hRlong) return false;
        if(hRout->GetEntries() < 1e-5 || hRside->GetEntries() < 1e-5 || hRlong->GetEntries() < 1e-5) return false;

        Rout = hRout->GetMean();
        dRout = hRout->GetStdDev();
        Rside = hRside->GetMean();
        dRside = hRside->GetStdDev();
        Rlong = hRlong->GetMean();
        dRlong = hRlong->GetStdDev();
        return true;
    };

    // CONTAINERS
    // Default systematics -> Lévy params
    TGraphAsymmErrors* alpha_default[NENERGIES] = {nullptr};
    TGraphAsymmErrors* R_default[NENERGIES] = {nullptr};
    TGraphAsymmErrors* N_default[NENERGIES] = {nullptr};
    TGraphAsymmErrors* R_out_default[NENERGIES] = {nullptr};
    TGraphAsymmErrors* R_side_default[NENERGIES] = {nullptr};
    TGraphAsymmErrors* R_long_default[NENERGIES] = {nullptr};
    // qLCMS systematics
    TGraphAsymmErrors* alpha_qlcms[NENERGIES][3] = {{nullptr}}; // 0: default, 1: strict, 2: loose
    TGraphAsymmErrors* R_qlcms[NENERGIES][3] = {{nullptr}};
    TGraphAsymmErrors* N_qlcms[NENERGIES][3] = {{nullptr}};
    TGraphAsymmErrors* R_out_qlcms[NENERGIES][3] = {{nullptr}};
    TGraphAsymmErrors* R_side_qlcms[NENERGIES][3] = {{nullptr}};
    TGraphAsymmErrors* R_long_qlcms[NENERGIES][3] = {{nullptr}};
    // rhofitmax systematics
    TGraphAsymmErrors* alpha_rhofitmax[NENERGIES][3] = {{nullptr}}; // 0: default, 1: strict, 2: loose
    TGraphAsymmErrors* R_rhofitmax[NENERGIES][3] = {{nullptr}};
    TGraphAsymmErrors* N_rhofitmax[NENERGIES][3] = {{nullptr}};
    TGraphAsymmErrors* R_out_rhofitmax[NENERGIES][3] = {{nullptr}};
    TGraphAsymmErrors* R_side_rhofitmax[NENERGIES][3] = {{nullptr}};
    TGraphAsymmErrors* R_long_rhofitmax[NENERGIES][3] = {{nullptr}};
    // NEVT_AVG systematics
    TGraphAsymmErrors* alpha_nevtavg[NENERGIES][3] = {{nullptr}}; // 0: default, 1: reasonably smaller, 2: reasonably larger
    TGraphAsymmErrors* R_nevtavg[NENERGIES][3] = {{nullptr}};
    TGraphAsymmErrors* N_nevtavg[NENERGIES][3] = {{nullptr}};
    TGraphAsymmErrors* R_out_nevtavg[NENERGIES][3] = {{nullptr}};
    TGraphAsymmErrors* R_side_nevtavg[NENERGIES][3] = {{nullptr}};
    TGraphAsymmErrors* R_long_nevtavg[NENERGIES][3] = {{nullptr}};

    // Main mT vs param graph
    TGraphAsymmErrors* alpha_syserr[NENERGIES] = {nullptr};
    TGraphAsymmErrors* R_syserr[NENERGIES] = {nullptr};
    TGraphAsymmErrors* N_syserr[NENERGIES] = {nullptr};
    TGraphAsymmErrors* R_out_syserr[NENERGIES] = {nullptr};
    TGraphAsymmErrors* R_side_syserr[NENERGIES] = {nullptr};
    TGraphAsymmErrors* R_long_syserr[NENERGIES] = {nullptr};

    // CSV lines to collect percent breakdown per energy
    std::vector<std::string> csv_lines;

    // Arrays to store per-energy percent contributions (average/median) for reuse in sqrt(sNN) plots
    double alpha_qlcms_pct_up_arr[NENERGIES] = {0}, alpha_qlcms_pct_dn_arr[NENERGIES] = {0};
    double alpha_rhofit_pct_up_arr[NENERGIES] = {0}, alpha_rhofit_pct_dn_arr[NENERGIES] = {0};
    double alpha_nevt_pct_up_arr[NENERGIES] = {0}, alpha_nevt_pct_dn_arr[NENERGIES] = {0};
    double R_qlcms_pct_up_arr[NENERGIES] = {0}, R_qlcms_pct_dn_arr[NENERGIES] = {0};
    double R_rhofit_pct_up_arr[NENERGIES] = {0}, R_rhofit_pct_dn_arr[NENERGIES] = {0};
    double R_nevt_pct_up_arr[NENERGIES] = {0}, R_nevt_pct_dn_arr[NENERGIES] = {0};
    double N_qlcms_pct_up_arr[NENERGIES] = {0}, N_qlcms_pct_dn_arr[NENERGIES] = {0};
    double N_rhofit_pct_up_arr[NENERGIES] = {0}, N_rhofit_pct_dn_arr[NENERGIES] = {0};
    double N_nevt_pct_up_arr[NENERGIES] = {0}, N_nevt_pct_dn_arr[NENERGIES] = {0};
    double R_out_qlcms_pct_up_arr[NENERGIES] = {0}, R_out_qlcms_pct_dn_arr[NENERGIES] = {0};
    double R_out_rhofit_pct_up_arr[NENERGIES] = {0}, R_out_rhofit_pct_dn_arr[NENERGIES] = {0};
    double R_out_nevt_pct_up_arr[NENERGIES] = {0}, R_out_nevt_pct_dn_arr[NENERGIES] = {0};
    double R_side_qlcms_pct_up_arr[NENERGIES] = {0}, R_side_qlcms_pct_dn_arr[NENERGIES] = {0};
    double R_side_rhofit_pct_up_arr[NENERGIES] = {0}, R_side_rhofit_pct_dn_arr[NENERGIES] = {0};
    double R_side_nevt_pct_up_arr[NENERGIES] = {0}, R_side_nevt_pct_dn_arr[NENERGIES] = {0};
    double R_long_qlcms_pct_up_arr[NENERGIES] = {0}, R_long_qlcms_pct_dn_arr[NENERGIES] = {0};
    double R_long_rhofit_pct_up_arr[NENERGIES] = {0}, R_long_rhofit_pct_dn_arr[NENERGIES] = {0};
    double R_long_nevt_pct_up_arr[NENERGIES] = {0}, R_long_nevt_pct_dn_arr[NENERGIES] = {0};
    double alpha_total_pct_up_arr[NENERGIES] = {0}, alpha_total_pct_dn_arr[NENERGIES] = {0};
    double R_total_pct_up_arr[NENERGIES] = {0}, R_total_pct_dn_arr[NENERGIES] = {0};
    double N_total_pct_up_arr[NENERGIES] = {0}, N_total_pct_dn_arr[NENERGIES] = {0};
    double R_out_total_pct_up_arr[NENERGIES] = {0}, R_out_total_pct_dn_arr[NENERGIES] = {0};
    double R_side_total_pct_up_arr[NENERGIES] = {0}, R_side_total_pct_dn_arr[NENERGIES] = {0};
    double R_long_total_pct_up_arr[NENERGIES] = {0}, R_long_total_pct_dn_arr[NENERGIES] = {0};

    // Also store mT-averaged default values and the ikt=2 etc. default points per energy
    double avg_alpha[NENERGIES] = {0}, avg_R[NENERGIES] = {0}, avg_N[NENERGIES] = {0};
    double avg_R_out[NENERGIES] = {0}, avg_R_side[NENERGIES] = {0}, avg_R_long[NENERGIES] = {0};
    double val_alpha_ikt2[NENERGIES] = {0}, val_R_ikt2[NENERGIES] = {0}, val_N_ikt2[NENERGIES] = {0};
    double val_R_out_ikt2[NENERGIES] = {0}, val_R_side_ikt2[NENERGIES] = {0}, val_R_long_ikt2[NENERGIES] = {0};
    // And the m_T choice syst uncert
    double alpha_mTchoice_syst_up[NENERGIES] = {0}, alpha_mTchoice_syst_dn[NENERGIES] = {0};
    double R_mTchoice_syst_up[NENERGIES] = {0}, R_mTchoice_syst_dn[NENERGIES] = {0};
    double N_mTchoice_syst_up[NENERGIES] = {0}, N_mTchoice_syst_dn[NENERGIES] = {0};
    double R_out_mTchoice_syst_up[NENERGIES] = {0}, R_out_mTchoice_syst_dn[NENERGIES] = {0};
    double R_side_mTchoice_syst_up[NENERGIES] = {0}, R_side_mTchoice_syst_dn[NENERGIES] = {0};
    double R_long_mTchoice_syst_up[NENERGIES] = {0}, R_long_mTchoice_syst_dn[NENERGIES] = {0};

    for(int ienergy = 0; ienergy < NENERGIES; ienergy++)
    {
        if(debug) cout << "Processing energy: " << energies[ienergy] << endl;
        // Default
        TFile* f_default;
        f_default = new TFile(Form("levyfit/results/UrQMD_onedfitresults_lcms_cent0-10_%s_AVG%d.root", 
                                        energies[ienergy], NEVT_AVGsyst[NEVT_AVG_DEFAULT[ienergy]]), "READ");
        if(!f_default || f_default->IsZombie())
        {
            printf("Could not open default file for energy %s: %s\n", energies[ienergy], f_default ? f_default->GetName() : "NULL");
            continue;
        }
        double alpha_vals[NKT] = {0};
        double alpha_errdn[NKT] = {0};
        double alpha_errup[NKT] = {0};
        double R_vals[NKT] = {0};
        double R_errup[NKT] = {0};
        double R_errdn[NKT] = {0};
        double N_vals[NKT] = {0};
        double N_errup[NKT] = {0};
        double N_errdn[NKT] = {0};
        double R_out_vals[NKT] = {0};
        double R_out_errup[NKT] = {0};
        double R_out_errdn[NKT] = {0};
        double R_side_vals[NKT] = {0};
        double R_side_errup[NKT] = {0};
        double R_side_errdn[NKT] = {0};
        double R_long_vals[NKT] = {0};
        double R_long_errup[NKT] = {0};
        double R_long_errdn[NKT] = {0};
        for(int ikt=0; ikt<NKT; ikt++)
        {
            TH2F* alphaR = (TH2F*)f_default->Get(Form("alpha_vs_R_ikt%i", ikt));
            if(!alphaR)
            {
                printf("Warning: [default] missing alpha_vs_R_ikt%i for energy %s. Skipping this point.\n", ikt, energies[ienergy]);
                continue;
            }
            if(alphaR->GetEntries() < 1e-5)
            {
                printf("Warning: [default] alpha_vs_R_ikt%i histogram for energy %s has no entries. Skipping this point.\n", ikt, energies[ienergy]);
                continue;
            } // Could do the same for N as well, but does not really matter since the more important point is done in correct_syserr_direction() where we compare the default to the syst checks, so if default is zero there, we skip the point and do not include it in the uncertainty calculation at all
            alpha_vals[ikt] = alphaR->GetMean(1);
            if(alpha_vals[ikt] < 1e-12)
            {
                printf("Warning: [default] alpha value for energy %s, ikt %i is very small (%g). Skipping this point.\n", energies[ienergy], ikt, alpha_vals[ikt]);
                continue;
            }
            alpha_errdn[ikt] = alphaR->GetStdDev(1);
            alpha_errup[ikt] = alphaR->GetStdDev(1);
            R_vals[ikt] = alphaR->GetMean(2);
            if(R_vals[ikt] < 1e-12)
            {
                printf("Warning: [default] R value for energy %s, ikt %i is very small (%g). Skipping this point.\n", energies[ienergy], ikt, R_vals[ikt]);
                continue;
            }
            R_errdn[ikt] = alphaR->GetStdDev(2);
            R_errup[ikt] = alphaR->GetStdDev(2);
            TH1F* hN = (TH1F*)f_default->Get(Form("Nhist_ikt%i", ikt));
            if(!hN)
            {
                printf("Warning: [default] missing Nhist_ikt%i for energy %s. Skipping this point.\n", ikt, energies[ienergy]);
                continue;
            }
            N_vals[ikt] = hN->GetMean();
            N_errdn[ikt] = hN->GetStdDev();
            N_errup[ikt] = hN->GetStdDev();
            fill_3d_radii_from_file(f_default, ikt,
                                    R_out_vals[ikt], R_out_errup[ikt],
                                    R_side_vals[ikt], R_side_errup[ikt],
                                    R_long_vals[ikt], R_long_errup[ikt]);
            R_out_errdn[ikt] = R_out_errup[ikt];
            R_side_errdn[ikt] = R_side_errup[ikt];
            R_long_errdn[ikt] = R_long_errup[ikt];
            if(debug)
            {
                cout << "Default." << endl;
                cout << Form("ikt %i: alpha = %g +%g -%g; R = %g +%g -%g; N = %g +%g -%g", ikt, alpha_vals[ikt], alpha_errup[ikt], alpha_errdn[ikt], R_vals[ikt], R_errup[ikt], R_errdn[ikt], N_vals[ikt], N_errup[ikt], N_errdn[ikt]) << endl;
                cout << Form("       R_out = %g +%g -%g; R_side = %g +%g -%g; R_long = %g +%g -%g", 
                                R_out_vals[ikt], R_out_errup[ikt], R_out_errdn[ikt],
                                R_side_vals[ikt], R_side_errup[ikt], R_side_errdn[ikt],
                                R_long_vals[ikt], R_long_errup[ikt], R_long_errdn[ikt]) << endl;
            }
        }
        alpha_default[ienergy] = new TGraphAsymmErrors(NKT, mtbin_centers, alpha_vals, xerr_low, xerr_high, alpha_errdn, alpha_errup);
        R_default[ienergy] = new TGraphAsymmErrors(NKT, mtbin_centers, R_vals, xerr_low, xerr_high, R_errdn, R_errup);
        N_default[ienergy] = new TGraphAsymmErrors(NKT, mtbin_centers, N_vals, xerr_low, xerr_high, N_errdn, N_errup);
        R_out_default[ienergy] = make_graph_if_populated(R_out_vals, R_out_errdn, R_out_errup);
        R_side_default[ienergy] = make_graph_if_populated(R_side_vals, R_side_errdn, R_side_errup);
        R_long_default[ienergy] = make_graph_if_populated(R_long_vals, R_long_errdn, R_long_errup);
        f_default->Close();
        // Reset arrays for next use
        for(int ikt=0; ikt<NKT; ikt++)
        {
            alpha_vals[ikt] = 0;
            alpha_errdn[ikt] = 0;
            alpha_errup[ikt] = 0;
            R_vals[ikt] = 0;
            R_errup[ikt] = 0;
            R_errdn[ikt] = 0;
            N_vals[ikt] = 0;
            N_errup[ikt] = 0;
            N_errdn[ikt] = 0;
            R_out_vals[ikt] = 0;
            R_out_errup[ikt] = 0;
            R_out_errdn[ikt] = 0;
            R_side_vals[ikt] = 0;
            R_side_errup[ikt] = 0;
            R_side_errdn[ikt] = 0;
            R_long_vals[ikt] = 0;
            R_long_errup[ikt] = 0;
            R_long_errdn[ikt] = 0;
        }

        // qLCMS
        const char* qlcms_labels[3] = {"","_strictqLCMS","_looseqLCMS"};
        for(int iq=0; iq<3; iq++)
        {
            TFile* f_qlcms;
            f_qlcms = new TFile(Form("levyfit/results/UrQMD_onedfitresults_lcms_cent0-10_%s_AVG%d%s.root", 
                                        energies[ienergy], NEVT_AVGsyst[NEVT_AVG_DEFAULT[ienergy]], qlcms_labels[iq]), "READ");
            if(!f_qlcms || f_qlcms->IsZombie())
            {
                printf("Could not open qLCMS file for energy %s: %s\n", energies[ienergy], f_qlcms ? f_qlcms->GetName() : "NULL");
                continue;
            }
            for(int ikt=0; ikt<NKT; ikt++)
            {
                TH2F* alphaR = (TH2F*)f_qlcms->Get(Form("alpha_vs_R_ikt%i", ikt));
                if(!alphaR)
                {
                    printf("Warning: [qLCMS %s] missing alpha_vs_R_ikt%i for energy %s. Skipping this point.\n", qlcms_labels[iq], ikt, energies[ienergy]);
                    continue;
                }
                if(alphaR->GetEntries() < 1e-5)
                {
                    printf("Warning: [qLCMS %s] alpha_vs_R_ikt%i histogram for energy %s has no entries. Skipping this point.\n", qlcms_labels[iq], ikt, energies[ienergy]);
                    continue;
                }
                alpha_vals[ikt] = alphaR->GetMean(1);
                alpha_errdn[ikt] = alphaR->GetStdDev(1);
                alpha_errup[ikt] = alphaR->GetStdDev(1);
                R_vals[ikt] = alphaR->GetMean(2);
                R_errdn[ikt] = alphaR->GetStdDev(2);
                R_errup[ikt] = alphaR->GetStdDev(2);
                TH1F* hN = (TH1F*)f_qlcms->Get(Form("Nhist_ikt%i", ikt));
                if(!hN)
                {
                    printf("Warning: [qLCMS %s] missing Nhist_ikt%i for energy %s. Skipping this point.\n", qlcms_labels[iq], ikt, energies[ienergy]);
                    continue;
                }
                N_vals[ikt] = hN->GetMean();
                N_errdn[ikt] = hN->GetStdDev();
                N_errup[ikt] = hN->GetStdDev();
                fill_3d_radii_from_file(f_qlcms, ikt,
                                        R_out_vals[ikt], R_out_errup[ikt],
                                        R_side_vals[ikt], R_side_errup[ikt],
                                        R_long_vals[ikt], R_long_errup[ikt]);
                R_out_errdn[ikt] = R_out_errup[ikt];
                R_side_errdn[ikt] = R_side_errup[ikt];
                R_long_errdn[ikt] = R_long_errup[ikt];

                if(debug)
                {
                    cout << "qLCMS variation: " << qlcms_labels[iq] << endl;
                    cout << Form("ikt %i: alpha = %g +%g -%g; R = %g +%g -%g; N = %g +%g -%g", ikt, alpha_vals[ikt], alpha_errup[ikt], alpha_errdn[ikt], R_vals[ikt], R_errup[ikt], R_errdn[ikt], N_vals[ikt], N_errup[ikt], N_errdn[ikt]) << endl;
                    cout << Form("       R_out = %g +%g -%g; R_side = %g +%g -%g; R_long = %g +%g -%g", 
                                    R_out_vals[ikt], R_out_errup[ikt], R_out_errdn[ikt],
                                    R_side_vals[ikt], R_side_errup[ikt], R_side_errdn[ikt],
                                    R_long_vals[ikt], R_long_errup[ikt], R_long_errdn[ikt]) << endl;
                }
            }
            alpha_qlcms[ienergy][iq] = new TGraphAsymmErrors(NKT, mtbin_centers, alpha_vals, xerr_low, xerr_high, alpha_errdn, alpha_errup);
            R_qlcms[ienergy][iq] = new TGraphAsymmErrors(NKT, mtbin_centers, R_vals, xerr_low, xerr_high, R_errdn, R_errup);
            N_qlcms[ienergy][iq] = new TGraphAsymmErrors(NKT, mtbin_centers, N_vals, xerr_low, xerr_high, N_errdn, N_errup);
            R_out_qlcms[ienergy][iq] = make_graph_if_populated(R_out_vals, R_out_errdn, R_out_errup);
            R_side_qlcms[ienergy][iq] = make_graph_if_populated(R_side_vals, R_side_errdn, R_side_errup);
            R_long_qlcms[ienergy][iq] = make_graph_if_populated(R_long_vals, R_long_errdn, R_long_errup);
            f_qlcms->Close();
            // Reset arrays for next use
            for(int ikt=0; ikt<NKT; ikt++)
            {
                alpha_vals[ikt] = 0;
                alpha_errdn[ikt] = 0;
                alpha_errup[ikt] = 0;
                R_vals[ikt] = 0;
                R_errup[ikt] = 0;
                R_errdn[ikt] = 0;
                N_vals[ikt] = 0;
                N_errup[ikt] = 0;
                N_errdn[ikt] = 0;
                R_out_vals[ikt] = 0;
                R_out_errup[ikt] = 0;
                R_out_errdn[ikt] = 0;
                R_side_vals[ikt] = 0;
                R_side_errup[ikt] = 0;
                R_side_errdn[ikt] = 0;
                R_long_vals[ikt] = 0;
                R_long_errup[ikt] = 0;
                R_long_errdn[ikt] = 0;
            }
        }

        // rhofitmax
        const char* rhofitmax_labels[3] = {"","_strictrhoFitMax","_looserhoFitMax"};
        for(int ir=0; ir<3; ir++)
        {
            TFile* f_rhofitmax;
            f_rhofitmax = new TFile(Form("levyfit/results/UrQMD_onedfitresults_lcms_cent0-10_%s_AVG%d%s.root", 
                                            energies[ienergy], NEVT_AVGsyst[NEVT_AVG_DEFAULT[ienergy]], rhofitmax_labels[ir]), "READ");
            if(!f_rhofitmax || f_rhofitmax->IsZombie())
            {
                printf("Could not open rhofitmax file for energy %s: %s\n", energies[ienergy], f_rhofitmax ? f_rhofitmax->GetName() : "NULL");
                continue;
            }
            for(int ikt=0; ikt<NKT; ikt++)
            {
                TH2F* alphaR = (TH2F*)f_rhofitmax->Get(Form("alpha_vs_R_ikt%i", ikt));
                if(!alphaR)
                {
                    printf("Warning: [rhofitmax %s] missing alpha_vs_R_ikt%i for energy %s. Skipping this point.\n", rhofitmax_labels[ir], ikt, energies[ienergy]);
                    continue;
                }
                if(alphaR->GetEntries() < 1e-5)
                {
                    printf("Warning: [rhofitmax %s] alpha_vs_R_ikt%i histogram for energy %s has no entries. Skipping this point.\n", rhofitmax_labels[ir], ikt, energies[ienergy]);
                    continue;
                }
                alpha_vals[ikt] = alphaR->GetMean(1);
                alpha_errdn[ikt] = alphaR->GetStdDev(1);
                alpha_errup[ikt] = alphaR->GetStdDev(1);
                R_vals[ikt] = alphaR->GetMean(2);
                R_errdn[ikt] = alphaR->GetStdDev(2);
                R_errup[ikt] = alphaR->GetStdDev(2);
                TH1F* hN = (TH1F*)f_rhofitmax->Get(Form("Nhist_ikt%i", ikt));
                if(!hN)
                {
                    printf("Warning: [rhofitmax %s] missing Nhist_ikt%i for energy %s. Skipping this point.\n", rhofitmax_labels[ir], ikt, energies[ienergy]);
                    continue;
                }
                N_vals[ikt] = hN->GetMean();
                N_errdn[ikt] = hN->GetStdDev();
                N_errup[ikt] = hN->GetStdDev();
                fill_3d_radii_from_file(f_rhofitmax, ikt,
                                        R_out_vals[ikt], R_out_errup[ikt],
                                        R_side_vals[ikt], R_side_errup[ikt],
                                        R_long_vals[ikt], R_long_errup[ikt]);
                R_out_errdn[ikt] = R_out_errup[ikt];
                R_side_errdn[ikt] = R_side_errup[ikt];
                R_long_errdn[ikt] = R_long_errup[ikt];
                if(debug)
                {
                    cout << "rhofitmax variation: " << rhofitmax_labels[ir] << endl;
                    cout << Form("ikt %i: alpha = %g +%g -%g; R = %g +%g -%g; N = %g +%g -%g", ikt, alpha_vals[ikt], alpha_errup[ikt], alpha_errdn[ikt], R_vals[ikt], R_errup[ikt], R_errdn[ikt], N_vals[ikt], N_errup[ikt], N_errdn[ikt]) << endl;
                    cout << Form("       R_out = %g +%g -%g; R_side = %g +%g -%g; R_long = %g +%g -%g", 
                                    R_out_vals[ikt], R_out_errup[ikt], R_out_errdn[ikt],
                                    R_side_vals[ikt], R_side_errup[ikt], R_side_errdn[ikt],
                                    R_long_vals[ikt], R_long_errup[ikt], R_long_errdn[ikt]) << endl;
                }
            }
            alpha_rhofitmax[ienergy][ir] = new TGraphAsymmErrors(NKT, mtbin_centers, alpha_vals, xerr_low, xerr_high, alpha_errdn, alpha_errup);
            R_rhofitmax[ienergy][ir] = new TGraphAsymmErrors(NKT, mtbin_centers, R_vals, xerr_low, xerr_high, R_errdn, R_errup);
            N_rhofitmax[ienergy][ir] = new TGraphAsymmErrors(NKT, mtbin_centers, N_vals, xerr_low, xerr_high, N_errdn, N_errup);
            R_out_rhofitmax[ienergy][ir] = make_graph_if_populated(R_out_vals, R_out_errdn, R_out_errup);
            R_side_rhofitmax[ienergy][ir] = make_graph_if_populated(R_side_vals, R_side_errdn, R_side_errup);
            R_long_rhofitmax[ienergy][ir] = make_graph_if_populated(R_long_vals, R_long_errdn, R_long_errup);
            f_rhofitmax->Close();
            // Reset arrays for next use
            for(int ikt=0; ikt<NKT; ikt++)
            {
                alpha_vals[ikt] = 0;
                alpha_errdn[ikt] = 0;
                alpha_errup[ikt] = 0;
                R_vals[ikt] = 0;
                R_errup[ikt] = 0;
                R_errdn[ikt] = 0;
                N_vals[ikt] = 0;
                N_errup[ikt] = 0;
                N_errdn[ikt] = 0;
                R_out_vals[ikt] = 0;
                R_out_errup[ikt] = 0;
                R_out_errdn[ikt] = 0;
                R_side_vals[ikt] = 0;
                R_side_errup[ikt] = 0;
                R_side_errdn[ikt] = 0;
                R_long_vals[ikt] = 0;
                R_long_errup[ikt] = 0;
                R_long_errdn[ikt] = 0;
            }
        }

        // NEVT_AVG
        for(int in=0; in<3; in++)
        {
            // Finding the correct index for the NEVT_AVG variation to read the correct file; this is a bit convoluted due to the way the variations are set up, but we want to make sure we are reading the intended variation files for each energy, especially since the high-stat energies have a different set of variations than the lower-stat ones
            int index = 0;
            {
                // Default
                index = NEVT_AVG_DEFAULT[ienergy];
                // too low || too high || too high for lowstat (as there the highest is 10000 FIXME if that changes)
                if(index<=0 || index > NEVTAVGS-1 || (ienergy > NENERGIES_highstat && NEVT_AVGsyst[index] > 5000))
                {
                    printf("Error: NEVT_AVG_DEFAULT index %d for energy %s is out of range (cannot choose systematic variation) for low-stat energies (min. index 1). Check your configuration in the common header!\n", index, energies[ienergy]);
                    continue;
                }
                // Smaller
                if(in==1)
                {
                    index=NEVT_AVG_DEFAULT[ienergy] - 1;
                }
                else if(in==2) // in==2 -> Larger
                {
                    index = NEVT_AVG_DEFAULT[ienergy] + 1;
                }
            } // Index selected; now read the file...
            cout << "Selected NEVT_AVG variation index for energy " << energies[ienergy] << ": index " << index << endl;
            
            TFile* f_nevtavg;
            f_nevtavg = new TFile(Form("levyfit/results/UrQMD_onedfitresults_lcms_cent0-10_%s_AVG%d.root", 
                                            energies[ienergy], NEVT_AVGsyst[index]), "READ");
            if(!f_nevtavg || f_nevtavg->IsZombie())
            {
                printf("Could not open nevt_avg file for energy %s: %s\n", energies[ienergy], f_nevtavg ? f_nevtavg->GetName() : "NULL");
                continue;
            }
            
            for(int ikt=0; ikt<NKT; ikt++)
            {
                TH2F* alphaR = (TH2F*)f_nevtavg->Get(Form("alpha_vs_R_ikt%i", ikt));
                if(!alphaR)
                {
                    printf("Warning: [nevt_avg %s] missing alpha_vs_R_ikt%i for energy %s. Skipping this point.\n", (in==0 ? "default" : (in==1 ? "smaller" : "larger")), ikt, energies[ienergy]);
                    continue;
                }
                if(alphaR->GetEntries() < 1e-5)
                {
                    printf("Warning: [nevt_avg %s] alpha_vs_R_ikt%i histogram for energy %s has no entries. Skipping this point.\n", (in==0 ? "default" : (in==1 ? "smaller" : "larger")), ikt, energies[ienergy]);
                    continue;
                }
                alpha_vals[ikt] = alphaR->GetMean(1);
                alpha_errdn[ikt] = alphaR->GetStdDev(1);
                alpha_errup[ikt] = alphaR->GetStdDev(1);
                R_vals[ikt] = alphaR->GetMean(2);
                R_errdn[ikt] = alphaR->GetStdDev(2);
                R_errup[ikt] = alphaR->GetStdDev(2);
                TH1F* hN = (TH1F*)f_nevtavg->Get(Form("Nhist_ikt%i", ikt));
                if(!hN)
                {
                    printf("Warning: [nevt_avg %s] missing Nhist_ikt%i for energy %s. Skipping this point.\n", (in==0 ? "default" : (in==1 ? "smaller" : "larger")), ikt, energies[ienergy]);
                    continue;
                }
                N_vals[ikt] = hN->GetMean();
                N_errdn[ikt] = hN->GetStdDev();
                N_errup[ikt] = hN->GetStdDev();
                fill_3d_radii_from_file(f_nevtavg, ikt,
                                        R_out_vals[ikt], R_out_errup[ikt],
                                        R_side_vals[ikt], R_side_errup[ikt],
                                        R_long_vals[ikt], R_long_errup[ikt]);
                R_out_errdn[ikt] = R_out_errup[ikt];
                R_side_errdn[ikt] = R_side_errup[ikt];
                R_long_errdn[ikt] = R_long_errup[ikt];
                if(debug)
                {
                    cout << "nevt_avg variation: " << (in==0 ? "default" : (in==1 ? "smaller" : "larger")) << endl;
                    cout << Form("ikt %i: alpha = %g +%g -%g; R = %g +%g -%g; N = %g +%g -%g", ikt, alpha_vals[ikt], alpha_errup[ikt], alpha_errdn[ikt], R_vals[ikt], R_errup[ikt], R_errdn[ikt], N_vals[ikt], N_errup[ikt], N_errdn[ikt]) << endl;
                    cout << Form("       R_out = %g +%g -%g; R_side = %g +%g -%g; R_long = %g +%g -%g", 
                                    R_out_vals[ikt], R_out_errup[ikt], R_out_errdn[ikt],
                                    R_side_vals[ikt], R_side_errup[ikt], R_side_errdn[ikt],
                                    R_long_vals[ikt], R_long_errup[ikt], R_long_errdn[ikt]) << endl;
                }
            }
            alpha_nevtavg[ienergy][in] = new TGraphAsymmErrors(NKT, mtbin_centers, alpha_vals, xerr_low, xerr_high, alpha_errdn, alpha_errup);
            R_nevtavg[ienergy][in] = new TGraphAsymmErrors(NKT, mtbin_centers, R_vals, xerr_low, xerr_high, R_errdn, R_errup);
            N_nevtavg[ienergy][in] = new TGraphAsymmErrors(NKT, mtbin_centers, N_vals, xerr_low, xerr_high, N_errdn, N_errup);
            R_out_nevtavg[ienergy][in] = make_graph_if_populated(R_out_vals, R_out_errdn, R_out_errup);
            R_side_nevtavg[ienergy][in] = make_graph_if_populated(R_side_vals, R_side_errdn, R_side_errup);
            R_long_nevtavg[ienergy][in] = make_graph_if_populated(R_long_vals, R_long_errdn, R_long_errup);
            f_nevtavg->Close();
            // Reset arrays for next use
            for(int ikt=0; ikt<NKT; ikt++)
            {
                alpha_vals[ikt] = 0;
                alpha_errdn[ikt] = 0;
                alpha_errup[ikt] = 0;
                R_vals[ikt] = 0;
                R_errup[ikt] = 0;
                R_errdn[ikt] = 0;
                N_vals[ikt] = 0;
                N_errup[ikt] = 0;
                N_errdn[ikt] = 0;
                R_out_vals[ikt] = 0;
                R_out_errup[ikt] = 0;
                R_out_errdn[ikt] = 0;
                R_side_vals[ikt] = 0;
                R_side_errup[ikt] = 0;
                R_side_errdn[ikt] = 0;
                R_long_vals[ikt] = 0;
                R_long_errup[ikt] = 0;
                R_long_errdn[ikt] = 0;
            }
        }

        // Now, having all systematics collected for this energy, calculate the total systematic uncertainties for each mT bin and Lévy parameter
        double alpha_syserr_up[NKT] = {0};
        double alpha_syserr_dn[NKT] = {0};
        double R_syserr_up[NKT] = {0};
        double R_syserr_dn[NKT] = {0};
        double N_syserr_up[NKT] = {0};
        double N_syserr_dn[NKT] = {0};
        double R_out_syserr_up[NKT] = {0};
        double R_out_syserr_dn[NKT] = {0};
        double R_side_syserr_up[NKT] = {0};
        double R_side_syserr_dn[NKT] = {0};
        double R_long_syserr_up[NKT] = {0};
        double R_long_syserr_dn[NKT] = {0};
        // For calculating individual contributions, later in percent
        double alpha_qlcms_up[NKT] = {0}, alpha_qlcms_dn[NKT] = {0};
        double alpha_rhofitmax_up[NKT] = {0}, alpha_rhofitmax_dn[NKT] = {0};
        double alpha_nevtavg_up[NKT] = {0}, alpha_nevtavg_dn[NKT] = {0};
        double R_qlcms_up[NKT] = {0}, R_qlcms_dn[NKT] = {0};
        double R_rhofitmax_up[NKT] = {0}, R_rhofitmax_dn[NKT] = {0};
        double R_nevtavg_up[NKT] = {0}, R_nevtavg_dn[NKT] = {0};
        double N_qlcms_up[NKT] = {0}, N_qlcms_dn[NKT] = {0};
        double N_rhofitmax_up[NKT] = {0}, N_rhofitmax_dn[NKT] = {0};
        double N_nevtavg_up[NKT] = {0}, N_nevtavg_dn[NKT] = {0};
        double R_out_qlcms_up[NKT] = {0}, R_out_qlcms_dn[NKT] = {0};
        double R_out_rhofitmax_up[NKT] = {0}, R_out_rhofitmax_dn[NKT] = {0};
        double R_out_nevtavg_up[NKT] = {0}, R_out_nevtavg_dn[NKT] = {0};
        double R_side_qlcms_up[NKT] = {0}, R_side_qlcms_dn[NKT] = {0};
        double R_side_rhofitmax_up[NKT] = {0}, R_side_rhofitmax_dn[NKT] = {0};
        double R_side_nevtavg_up[NKT] = {0}, R_side_nevtavg_dn[NKT] = {0};
        double R_long_qlcms_up[NKT] = {0}, R_long_qlcms_dn[NKT] = {0};
        double R_long_rhofitmax_up[NKT] = {0}, R_long_rhofitmax_dn[NKT] = {0};
        double R_long_nevtavg_up[NKT] = {0}, R_long_nevtavg_dn[NKT] = {0};

        auto accumulate_syserr_if_ready = [&](double* up, double* dn,
                                              TGraphAsymmErrors* syst1,
                                              TGraphAsymmErrors* syst2,
                                              TGraphAsymmErrors* def,
                                              int ikt)
        {
            if(!syst1 || !syst2 || !def) return;
            correct_syserr_direction(up, dn, syst1, syst2, def, ikt);
        };
        for(int ikt=0; ikt<NKT; ikt++)
        {
            // Collect deviations from default for each systematic source - check the general direction of deviations (up/down)
            accumulate_syserr_if_ready(alpha_syserr_up, alpha_syserr_dn, alpha_qlcms[ienergy][1], alpha_qlcms[ienergy][2], alpha_default[ienergy], ikt); // loose, strict qLCMS
            // TODO deal with saving individual syst. unc. contributions as well
            accumulate_syserr_if_ready(alpha_syserr_up, alpha_syserr_dn, alpha_rhofitmax[ienergy][1], alpha_rhofitmax[ienergy][2], alpha_default[ienergy], ikt); // loose, strict rhofitmax
            accumulate_syserr_if_ready(alpha_syserr_up, alpha_syserr_dn, alpha_nevtavg[ienergy][1], alpha_nevtavg[ienergy][2], alpha_default[ienergy], ikt); // smaller, larger nevt_avg
            alpha_syserr_up[ikt] = sqrt(alpha_syserr_up[ikt]);
            alpha_syserr_dn[ikt] = sqrt(alpha_syserr_dn[ikt]);
            // cerr << "alpha_syserr_up: " << alpha_syserr_up[ikt] << ", alpha_syserr_dn: " << alpha_syserr_dn[ikt] << endl;

            accumulate_syserr_if_ready(R_syserr_up, R_syserr_dn, R_qlcms[ienergy][1], R_qlcms[ienergy][2], R_default[ienergy], ikt); // strict, loose qLCMS
            accumulate_syserr_if_ready(R_syserr_up, R_syserr_dn, R_rhofitmax[ienergy][1], R_rhofitmax[ienergy][2], R_default[ienergy], ikt); // strict, loose rhofitmax
            accumulate_syserr_if_ready(R_syserr_up, R_syserr_dn, R_nevtavg[ienergy][1], R_nevtavg[ienergy][2], R_default[ienergy], ikt); // smaller, larger nevt_avg
            R_syserr_up[ikt] = sqrt(R_syserr_up[ikt]);
            R_syserr_dn[ikt] = sqrt(R_syserr_dn[ikt]);
            // cerr << "R_syserr_up: " << R_syserr_up[ikt] << ", R_syserr_dn: " << R_syserr_dn[ikt] << endl;

            accumulate_syserr_if_ready(N_syserr_up, N_syserr_dn, N_qlcms[ienergy][1], N_qlcms[ienergy][2], N_default[ienergy], ikt); // strict, loose qLCMS
            accumulate_syserr_if_ready(N_syserr_up, N_syserr_dn, N_rhofitmax[ienergy][1], N_rhofitmax[ienergy][2], N_default[ienergy], ikt); // strict, loose rhofitmax
            accumulate_syserr_if_ready(N_syserr_up, N_syserr_dn, N_nevtavg[ienergy][1], N_nevtavg[ienergy][2], N_default[ienergy], ikt); // smaller, larger nevt_avg
            N_syserr_up[ikt] = sqrt(N_syserr_up[ikt]);
            N_syserr_dn[ikt] = sqrt(N_syserr_dn[ikt]);
            // cerr << "N_syserr_up: " << N_syserr_up[ikt] << ", N_syserr_dn: " << N_syserr_dn[ikt] << endl;

            accumulate_syserr_if_ready(R_out_syserr_up, R_out_syserr_dn, R_out_qlcms[ienergy][1], R_out_qlcms[ienergy][2], R_out_default[ienergy], ikt);
            accumulate_syserr_if_ready(R_out_syserr_up, R_out_syserr_dn, R_out_rhofitmax[ienergy][1], R_out_rhofitmax[ienergy][2], R_out_default[ienergy], ikt);
            accumulate_syserr_if_ready(R_out_syserr_up, R_out_syserr_dn, R_out_nevtavg[ienergy][1], R_out_nevtavg[ienergy][2], R_out_default[ienergy], ikt);
            R_out_syserr_up[ikt] = sqrt(R_out_syserr_up[ikt]);
            R_out_syserr_dn[ikt] = sqrt(R_out_syserr_dn[ikt]);

            accumulate_syserr_if_ready(R_side_syserr_up, R_side_syserr_dn, R_side_qlcms[ienergy][1], R_side_qlcms[ienergy][2], R_side_default[ienergy], ikt);
            accumulate_syserr_if_ready(R_side_syserr_up, R_side_syserr_dn, R_side_rhofitmax[ienergy][1], R_side_rhofitmax[ienergy][2], R_side_default[ienergy], ikt);
            accumulate_syserr_if_ready(R_side_syserr_up, R_side_syserr_dn, R_side_nevtavg[ienergy][1], R_side_nevtavg[ienergy][2], R_side_default[ienergy], ikt);
            R_side_syserr_up[ikt] = sqrt(R_side_syserr_up[ikt]);
            R_side_syserr_dn[ikt] = sqrt(R_side_syserr_dn[ikt]);

            accumulate_syserr_if_ready(R_long_syserr_up, R_long_syserr_dn, R_long_qlcms[ienergy][1], R_long_qlcms[ienergy][2], R_long_default[ienergy], ikt);
            accumulate_syserr_if_ready(R_long_syserr_up, R_long_syserr_dn, R_long_rhofitmax[ienergy][1], R_long_rhofitmax[ienergy][2], R_long_default[ienergy], ikt);
            accumulate_syserr_if_ready(R_long_syserr_up, R_long_syserr_dn, R_long_nevtavg[ienergy][1], R_long_nevtavg[ienergy][2], R_long_default[ienergy], ikt);
            R_long_syserr_up[ikt] = sqrt(R_long_syserr_up[ikt]);
            R_long_syserr_dn[ikt] = sqrt(R_long_syserr_dn[ikt]);
        }
        // Create TGraphAsymmErrors for total systematic uncertainties
        alpha_syserr[ienergy] = new TGraphAsymmErrors(NKT, mtbin_centers, alpha_default[ienergy]->GetY(), xerr_low, xerr_high, alpha_syserr_dn, alpha_syserr_up);
        R_syserr[ienergy] = new TGraphAsymmErrors(NKT, mtbin_centers, R_default[ienergy]->GetY(), xerr_low, xerr_high, R_syserr_dn, R_syserr_up);
        N_syserr[ienergy] = new TGraphAsymmErrors(NKT, mtbin_centers, N_default[ienergy]->GetY(), xerr_low, xerr_high, N_syserr_dn, N_syserr_up);
        R_out_syserr[ienergy] = R_out_default[ienergy] ? new TGraphAsymmErrors(NKT, mtbin_centers, R_out_default[ienergy]->GetY(), xerr_low, xerr_high, R_out_syserr_dn, R_out_syserr_up) : nullptr;
        R_side_syserr[ienergy] = R_side_default[ienergy] ? new TGraphAsymmErrors(NKT, mtbin_centers, R_side_default[ienergy]->GetY(), xerr_low, xerr_high, R_side_syserr_dn, R_side_syserr_up) : nullptr;
        R_long_syserr[ienergy] = R_long_default[ienergy] ? new TGraphAsymmErrors(NKT, mtbin_centers, R_long_default[ienergy]->GetY(), xerr_low, xerr_high, R_long_syserr_dn, R_long_syserr_up) : nullptr;

        // Give one-one number about the systematic uncertainty contributions from each source

        // Firstly, extract the individual contributions as well
        for(int ikt=0; ikt<NKT; ikt++)
        {
            accumulate_syserr_if_ready(alpha_qlcms_up, alpha_qlcms_dn, alpha_qlcms[ienergy][1], alpha_qlcms[ienergy][2], alpha_default[ienergy], ikt);
            accumulate_syserr_if_ready(alpha_rhofitmax_up, alpha_rhofitmax_dn, alpha_rhofitmax[ienergy][1], alpha_rhofitmax[ienergy][2], alpha_default[ienergy], ikt);
            accumulate_syserr_if_ready(alpha_nevtavg_up, alpha_nevtavg_dn, alpha_nevtavg[ienergy][1], alpha_nevtavg[ienergy][2], alpha_default[ienergy], ikt);
            accumulate_syserr_if_ready(R_qlcms_up, R_qlcms_dn, R_qlcms[ienergy][1], R_qlcms[ienergy][2], R_default[ienergy], ikt);
            accumulate_syserr_if_ready(R_rhofitmax_up, R_rhofitmax_dn, R_rhofitmax[ienergy][1], R_rhofitmax[ienergy][2], R_default[ienergy], ikt);
            accumulate_syserr_if_ready(R_nevtavg_up, R_nevtavg_dn, R_nevtavg[ienergy][1], R_nevtavg[ienergy][2], R_default[ienergy], ikt);
            accumulate_syserr_if_ready(N_qlcms_up, N_qlcms_dn, N_qlcms[ienergy][1], N_qlcms[ienergy][2], N_default[ienergy], ikt);
            accumulate_syserr_if_ready(N_rhofitmax_up, N_rhofitmax_dn, N_rhofitmax[ienergy][1], N_rhofitmax[ienergy][2], N_default[ienergy], ikt);
            accumulate_syserr_if_ready(N_nevtavg_up, N_nevtavg_dn, N_nevtavg[ienergy][1], N_nevtavg[ienergy][2], N_default[ienergy], ikt);
            accumulate_syserr_if_ready(R_out_qlcms_up, R_out_qlcms_dn, R_out_qlcms[ienergy][1], R_out_qlcms[ienergy][2], R_out_default[ienergy], ikt);
            accumulate_syserr_if_ready(R_out_rhofitmax_up, R_out_rhofitmax_dn, R_out_rhofitmax[ienergy][1], R_out_rhofitmax[ienergy][2], R_out_default[ienergy], ikt);
            accumulate_syserr_if_ready(R_out_nevtavg_up, R_out_nevtavg_dn, R_out_nevtavg[ienergy][1], R_out_nevtavg[ienergy][2], R_out_default[ienergy], ikt);
            accumulate_syserr_if_ready(R_side_qlcms_up, R_side_qlcms_dn, R_side_qlcms[ienergy][1], R_side_qlcms[ienergy][2], R_side_default[ienergy], ikt);
            accumulate_syserr_if_ready(R_side_rhofitmax_up, R_side_rhofitmax_dn, R_side_rhofitmax[ienergy][1], R_side_rhofitmax[ienergy][2], R_side_default[ienergy], ikt);
            accumulate_syserr_if_ready(R_side_nevtavg_up, R_side_nevtavg_dn, R_side_nevtavg[ienergy][1], R_side_nevtavg[ienergy][2], R_side_default[ienergy], ikt);
            accumulate_syserr_if_ready(R_long_qlcms_up, R_long_qlcms_dn, R_long_qlcms[ienergy][1], R_long_qlcms[ienergy][2], R_long_default[ienergy], ikt);
            accumulate_syserr_if_ready(R_long_rhofitmax_up, R_long_rhofitmax_dn, R_long_rhofitmax[ienergy][1], R_long_rhofitmax[ienergy][2], R_long_default[ienergy], ikt);
            accumulate_syserr_if_ready(R_long_nevtavg_up, R_long_nevtavg_dn, R_long_nevtavg[ienergy][1], R_long_nevtavg[ienergy][2], R_long_default[ienergy], ikt);
        }

        // Individual source arrays hold summed squares; convert to magnitudes before percent conversion.
        for(int ikt=0; ikt<NKT; ikt++)
        {
            alpha_qlcms_up[ikt] = sqrt(alpha_qlcms_up[ikt]);
            alpha_qlcms_dn[ikt] = sqrt(alpha_qlcms_dn[ikt]);
            alpha_rhofitmax_up[ikt] = sqrt(alpha_rhofitmax_up[ikt]);
            alpha_rhofitmax_dn[ikt] = sqrt(alpha_rhofitmax_dn[ikt]);
            alpha_nevtavg_up[ikt] = sqrt(alpha_nevtavg_up[ikt]);
            alpha_nevtavg_dn[ikt] = sqrt(alpha_nevtavg_dn[ikt]);

            R_qlcms_up[ikt] = sqrt(R_qlcms_up[ikt]);
            R_qlcms_dn[ikt] = sqrt(R_qlcms_dn[ikt]);
            R_rhofitmax_up[ikt] = sqrt(R_rhofitmax_up[ikt]);
            R_rhofitmax_dn[ikt] = sqrt(R_rhofitmax_dn[ikt]);
            R_nevtavg_up[ikt] = sqrt(R_nevtavg_up[ikt]);
            R_nevtavg_dn[ikt] = sqrt(R_nevtavg_dn[ikt]);

            N_qlcms_up[ikt] = sqrt(N_qlcms_up[ikt]);
            N_qlcms_dn[ikt] = sqrt(N_qlcms_dn[ikt]);
            N_rhofitmax_up[ikt] = sqrt(N_rhofitmax_up[ikt]);
            N_rhofitmax_dn[ikt] = sqrt(N_rhofitmax_dn[ikt]);
            N_nevtavg_up[ikt] = sqrt(N_nevtavg_up[ikt]);
            N_nevtavg_dn[ikt] = sqrt(N_nevtavg_dn[ikt]);

            R_out_qlcms_up[ikt] = sqrt(R_out_qlcms_up[ikt]);
            R_out_qlcms_dn[ikt] = sqrt(R_out_qlcms_dn[ikt]);
            R_out_rhofitmax_up[ikt] = sqrt(R_out_rhofitmax_up[ikt]);
            R_out_rhofitmax_dn[ikt] = sqrt(R_out_rhofitmax_dn[ikt]);
            R_out_nevtavg_up[ikt] = sqrt(R_out_nevtavg_up[ikt]);
            R_out_nevtavg_dn[ikt] = sqrt(R_out_nevtavg_dn[ikt]);

            R_side_qlcms_up[ikt] = sqrt(R_side_qlcms_up[ikt]);
            R_side_qlcms_dn[ikt] = sqrt(R_side_qlcms_dn[ikt]);
            R_side_rhofitmax_up[ikt] = sqrt(R_side_rhofitmax_up[ikt]);
            R_side_rhofitmax_dn[ikt] = sqrt(R_side_rhofitmax_dn[ikt]);
            R_side_nevtavg_up[ikt] = sqrt(R_side_nevtavg_up[ikt]);
            R_side_nevtavg_dn[ikt] = sqrt(R_side_nevtavg_dn[ikt]);

            R_long_qlcms_up[ikt] = sqrt(R_long_qlcms_up[ikt]);
            R_long_qlcms_dn[ikt] = sqrt(R_long_qlcms_dn[ikt]);
            R_long_rhofitmax_up[ikt] = sqrt(R_long_rhofitmax_up[ikt]);
            R_long_rhofitmax_dn[ikt] = sqrt(R_long_rhofitmax_dn[ikt]);
            R_long_nevtavg_up[ikt] = sqrt(R_long_nevtavg_up[ikt]);
            R_long_nevtavg_dn[ikt] = sqrt(R_long_nevtavg_dn[ikt]);
        }
        
        // Average with calc_syst_contribution
        // Median with calc_syst_contrib_median, might make more sense?
        double alpha_qlcms_percent_up = calc_syst_contribution(alpha_default[ienergy], alpha_qlcms_up, alpha_qlcms_dn);
        double alpha_rhofitmax_percent_up = calc_syst_contribution(alpha_default[ienergy], alpha_rhofitmax_up, alpha_rhofitmax_dn);
        double alpha_nevtavg_percent_up = calc_syst_contribution(alpha_default[ienergy], alpha_nevtavg_up, alpha_nevtavg_dn);
        double R_qlcms_percent_up = calc_syst_contribution(R_default[ienergy], R_qlcms_up, R_qlcms_dn);
        double R_rhofitmax_percent_up = calc_syst_contribution(R_default[ienergy], R_rhofitmax_up, R_rhofitmax_dn);
        double R_nevtavg_percent_up = calc_syst_contribution(R_default[ienergy], R_nevtavg_up, R_nevtavg_dn);
        double N_qlcms_percent_up = calc_syst_contribution(N_default[ienergy], N_qlcms_up, N_qlcms_dn);
        double N_rhofitmax_percent_up = calc_syst_contribution(N_default[ienergy], N_rhofitmax_up, N_rhofitmax_dn);
        double N_nevtavg_percent_up = calc_syst_contribution(N_default[ienergy], N_nevtavg_up, N_nevtavg_dn);
        double R_out_qlcms_percent_up = calc_syst_contribution(R_out_default[ienergy], R_out_qlcms_up, R_out_qlcms_dn);
        double R_out_rhofitmax_percent_up = calc_syst_contribution(R_out_default[ienergy], R_out_rhofitmax_up, R_out_rhofitmax_dn);
        double R_out_nevtavg_percent_up = calc_syst_contribution(R_out_default[ienergy], R_out_nevtavg_up, R_out_nevtavg_dn);
        double R_side_qlcms_percent_up = calc_syst_contribution(R_side_default[ienergy], R_side_qlcms_up, R_side_qlcms_dn);
        double R_side_rhofitmax_percent_up = calc_syst_contribution(R_side_default[ienergy], R_side_rhofitmax_up, R_side_rhofitmax_dn);
        double R_side_nevtavg_percent_up = calc_syst_contribution(R_side_default[ienergy], R_side_nevtavg_up, R_side_nevtavg_dn);
        double R_long_qlcms_percent_up = calc_syst_contribution(R_long_default[ienergy], R_long_qlcms_up, R_long_qlcms_dn);
        double R_long_rhofitmax_percent_up = calc_syst_contribution(R_long_default[ienergy], R_long_rhofitmax_up, R_long_rhofitmax_dn);
        double R_long_nevtavg_percent_up = calc_syst_contribution(R_long_default[ienergy], R_long_nevtavg_up, R_long_nevtavg_dn);
        double alpha_qlcms_percent_dn = calc_syst_contribution(alpha_default[ienergy], alpha_qlcms_up, alpha_qlcms_dn, 0);
        double alpha_rhofitmax_percent_dn = calc_syst_contribution(alpha_default[ienergy], alpha_rhofitmax_up, alpha_rhofitmax_dn, 0);
        double alpha_nevtavg_percent_dn = calc_syst_contribution(alpha_default[ienergy], alpha_nevtavg_up, alpha_nevtavg_dn, 0);
        double R_qlcms_percent_dn = calc_syst_contribution(R_default[ienergy], R_qlcms_up, R_qlcms_dn, 0);
        double R_rhofitmax_percent_dn = calc_syst_contribution(R_default[ienergy], R_rhofitmax_up, R_rhofitmax_dn, 0);
        double R_nevtavg_percent_dn = calc_syst_contribution(R_default[ienergy], R_nevtavg_up, R_nevtavg_dn, 0);
        double N_qlcms_percent_dn = calc_syst_contribution(N_default[ienergy], N_qlcms_up, N_qlcms_dn, 0);
        double N_rhofitmax_percent_dn = calc_syst_contribution(N_default[ienergy], N_rhofitmax_up, N_rhofitmax_dn, 0);
        double N_nevtavg_percent_dn = calc_syst_contribution(N_default[ienergy], N_nevtavg_up, N_nevtavg_dn, 0);
        double R_out_qlcms_percent_dn = calc_syst_contribution(R_out_default[ienergy], R_out_qlcms_up, R_out_qlcms_dn, 0);
        double R_out_rhofitmax_percent_dn = calc_syst_contribution(R_out_default[ienergy], R_out_rhofitmax_up, R_out_rhofitmax_dn, 0);
        double R_out_nevtavg_percent_dn = calc_syst_contribution(R_out_default[ienergy], R_out_nevtavg_up, R_out_nevtavg_dn, 0);
        double R_side_qlcms_percent_dn = calc_syst_contribution(R_side_default[ienergy], R_side_qlcms_up, R_side_qlcms_dn, 0);
        double R_side_rhofitmax_percent_dn = calc_syst_contribution(R_side_default[ienergy], R_side_rhofitmax_up, R_side_rhofitmax_dn, 0);
        double R_side_nevtavg_percent_dn = calc_syst_contribution(R_side_default[ienergy], R_side_nevtavg_up, R_side_nevtavg_dn, 0);
        double R_long_qlcms_percent_dn = calc_syst_contribution(R_long_default[ienergy], R_long_qlcms_up, R_long_qlcms_dn, 0);
        double R_long_rhofitmax_percent_dn = calc_syst_contribution(R_long_default[ienergy], R_long_rhofitmax_up, R_long_rhofitmax_dn, 0);
        double R_long_nevtavg_percent_dn = calc_syst_contribution(R_long_default[ienergy], R_long_nevtavg_up, R_long_nevtavg_dn, 0);

        double alpha_total_percent_up = calc_syst_contribution(alpha_default[ienergy], alpha_syserr_up, alpha_syserr_dn);
        double alpha_total_percent_dn = calc_syst_contribution(alpha_default[ienergy], alpha_syserr_up, alpha_syserr_dn, 0);
        double R_total_percent_up = calc_syst_contribution(R_default[ienergy], R_syserr_up, R_syserr_dn);
        double R_total_percent_dn = calc_syst_contribution(R_default[ienergy], R_syserr_up, R_syserr_dn, 0);
        double N_total_percent_up = calc_syst_contribution(N_default[ienergy], N_syserr_up, N_syserr_dn);
        double N_total_percent_dn = calc_syst_contribution(N_default[ienergy], N_syserr_up, N_syserr_dn, 0);
        double R_out_total_percent_up = calc_syst_contribution(R_out_default[ienergy], R_out_syserr_up, R_out_syserr_dn);
        double R_out_total_percent_dn = calc_syst_contribution(R_out_default[ienergy], R_out_syserr_up, R_out_syserr_dn, 0);
        double R_side_total_percent_up = calc_syst_contribution(R_side_default[ienergy], R_side_syserr_up, R_side_syserr_dn);
        double R_side_total_percent_dn = calc_syst_contribution(R_side_default[ienergy], R_side_syserr_up, R_side_syserr_dn, 0);
        double R_long_total_percent_up = calc_syst_contribution(R_long_default[ienergy], R_long_syserr_up, R_long_syserr_dn);
        double R_long_total_percent_dn = calc_syst_contribution(R_long_default[ienergy], R_long_syserr_up, R_long_syserr_dn, 0);

        // Append CSV block lines: energy header, parameter sub-header, then lines for alpha,R,N
        std::ostringstream h1;
        h1 << "[" << energies[ienergy] << " GeV]";
        csv_lines.push_back(h1.str());

        csv_lines.push_back("Levy parameter,qlcms cut [%],rhofitmax [%],nevtavg [%],total [%]");
        std::ostringstream arow;
        arow << "alpha," << std::fixed << std::setprecision(6)
            << " +" << alpha_qlcms_percent_up << " / -" << alpha_qlcms_percent_dn << ","
            << " +" << alpha_rhofitmax_percent_up << " / -" << alpha_rhofitmax_percent_dn << ","
            << " +" << alpha_nevtavg_percent_up << " / -" << alpha_nevtavg_percent_dn << ","
            << " +" << alpha_total_percent_up << " / -" << alpha_total_percent_dn;
        csv_lines.push_back(arow.str());
        std::ostringstream rrow;
        rrow << "R," << std::fixed << std::setprecision(6)
            << " +" << R_qlcms_percent_up << " / -" << R_qlcms_percent_dn << ","
            << " +" << R_rhofitmax_percent_up << " / -" << R_rhofitmax_percent_dn << ","
            << " +" << R_nevtavg_percent_up << " / -" << R_nevtavg_percent_dn << ","
            << " +" << R_total_percent_up << " / -" << R_total_percent_dn;
        csv_lines.push_back(rrow.str());
        std::ostringstream nrow;
        nrow << "N," << std::fixed << std::setprecision(6)
            << " +" << N_qlcms_percent_up << " / -" << N_qlcms_percent_dn << ","
            << " +" << N_rhofitmax_percent_up << " / -" << N_rhofitmax_percent_dn << ","
            << " +" << N_nevtavg_percent_up << " / -" << N_nevtavg_percent_dn << ","
            << " +" << N_total_percent_up << " / -" << N_total_percent_dn;
        csv_lines.push_back(nrow.str());
        std::ostringstream routrow;
        routrow << "R_out," << std::fixed << std::setprecision(6)
            << " +" << R_out_qlcms_percent_up << " / -" << R_out_qlcms_percent_dn << ","
            << " +" << R_out_rhofitmax_percent_up << " / -" << R_out_rhofitmax_percent_dn << ","
            << " +" << R_out_nevtavg_percent_up << " / -" << R_out_nevtavg_percent_dn << ","
            << " +" << R_out_total_percent_up << " / -" << R_out_total_percent_dn;
        csv_lines.push_back(routrow.str());
        std::ostringstream rsiderow;
        rsiderow << "R_side," << std::fixed << std::setprecision(6)
             << " +" << R_side_qlcms_percent_up << " / -" << R_side_qlcms_percent_dn << ","
             << " +" << R_side_rhofitmax_percent_up << " / -" << R_side_rhofitmax_percent_dn << ","
             << " +" << R_side_nevtavg_percent_up << " / -" << R_side_nevtavg_percent_dn << ","
             << " +" << R_side_total_percent_up << " / -" << R_side_total_percent_dn;
        csv_lines.push_back(rsiderow.str());
        std::ostringstream rlongrow;
        rlongrow << "R_long," << std::fixed << std::setprecision(6)
             << " +" << R_long_qlcms_percent_up << " / -" << R_long_qlcms_percent_dn << ","
             << " +" << R_long_rhofitmax_percent_up << " / -" << R_long_rhofitmax_percent_dn << ","
             << " +" << R_long_nevtavg_percent_up << " / -" << R_long_nevtavg_percent_dn << ","
             << " +" << R_long_total_percent_up << " / -" << R_long_total_percent_dn;
        csv_lines.push_back(rlongrow.str());
        
        // ---- m_T choice systematic contribution (percent) ----
        // Add separate energy header for m_T choice block
        std::ostringstream h_mt;
        h_mt << "[" << energies[ienergy] << " GeV]";
        csv_lines.push_back(h_mt.str());
        
        // compute mT center and neighbor differences in MeV
        double default_mt_mev = mtbin_centers[2] * 1000.0;
        double upper_diff_mev = (mtbin_centers[3] - mtbin_centers[2]) * 1000.0;
        double lower_diff_mev = (mtbin_centers[2] - mtbin_centers[1]) * 1000.0;
        double mt_pct_up = 0.0, mt_pct_dn = 0.0;
        TGraphAsymmErrors* gdef_mt = alpha_default[ienergy];
        if(gdef_mt && gdef_mt->GetN() > 2){
            double center = gdef_mt->GetY()[2];
            double neigh_up = 0., neigh_dn = 0.;
            calc_neighbor_uncert(gdef_mt, 2, neigh_up, neigh_dn);
            if(std::isfinite(center) && center != 0.){
                // Leaves mT choice uncertainty zero if turned off via mtchoice_syst=false
                if(mtchoice_syst) mt_pct_up = 100.0 * neigh_up / fabs(center);
                if(mtchoice_syst) mt_pct_dn = 100.0 * neigh_dn / fabs(center);
            }
        }
        
        // Subheader with m_T details
        std::ostringstream mtsub;
        mtsub << "m_T choice (m_T=" << std::fixed << std::setprecision(0) << default_mt_mev
              << " MeV +" << std::fixed << std::setprecision(0) << upper_diff_mev
              << "/-" << std::fixed << std::setprecision(0) << lower_diff_mev << " MeV)";
        csv_lines.push_back(mtsub.str());
        
        // Data line with energy and percentages
        std::ostringstream mtline;
        mtline << energies[ienergy] << "," << std::fixed << std::setprecision(6)
               << " +" << mt_pct_up << " / -" << mt_pct_dn;
        csv_lines.push_back(mtline.str());
        
        // Store per-energy percent averages/medians for later sqrt(sNN) plotting
        alpha_qlcms_pct_up_arr[ienergy] = alpha_qlcms_percent_up;
        alpha_qlcms_pct_dn_arr[ienergy] = alpha_qlcms_percent_dn;
        alpha_rhofit_pct_up_arr[ienergy] = alpha_rhofitmax_percent_up;
        alpha_rhofit_pct_dn_arr[ienergy] = alpha_rhofitmax_percent_dn;
        alpha_nevt_pct_up_arr[ienergy] = alpha_nevtavg_percent_up;
        alpha_nevt_pct_dn_arr[ienergy] = alpha_nevtavg_percent_dn;

        R_qlcms_pct_up_arr[ienergy] = R_qlcms_percent_up;
        R_qlcms_pct_dn_arr[ienergy] = R_qlcms_percent_dn;
        R_rhofit_pct_up_arr[ienergy] = R_rhofitmax_percent_up;
        R_rhofit_pct_dn_arr[ienergy] = R_rhofitmax_percent_dn;
        R_nevt_pct_up_arr[ienergy] = R_nevtavg_percent_up;
        R_nevt_pct_dn_arr[ienergy] = R_nevtavg_percent_dn;

        N_qlcms_pct_up_arr[ienergy] = N_qlcms_percent_up;
        N_qlcms_pct_dn_arr[ienergy] = N_qlcms_percent_dn;
        N_rhofit_pct_up_arr[ienergy] = N_rhofitmax_percent_up;
        N_rhofit_pct_dn_arr[ienergy] = N_rhofitmax_percent_dn;
        N_nevt_pct_up_arr[ienergy] = N_nevtavg_percent_up;
        N_nevt_pct_dn_arr[ienergy] = N_nevtavg_percent_dn;

        R_out_qlcms_pct_up_arr[ienergy] = R_out_qlcms_percent_up;
        R_out_qlcms_pct_dn_arr[ienergy] = R_out_qlcms_percent_dn;
        R_out_rhofit_pct_up_arr[ienergy] = R_out_rhofitmax_percent_up;
        R_out_rhofit_pct_dn_arr[ienergy] = R_out_rhofitmax_percent_dn;
        R_out_nevt_pct_up_arr[ienergy] = R_out_nevtavg_percent_up;
        R_out_nevt_pct_dn_arr[ienergy] = R_out_nevtavg_percent_dn;

        R_side_qlcms_pct_up_arr[ienergy] = R_side_qlcms_percent_up;
        R_side_qlcms_pct_dn_arr[ienergy] = R_side_qlcms_percent_dn;
        R_side_rhofit_pct_up_arr[ienergy] = R_side_rhofitmax_percent_up;
        R_side_rhofit_pct_dn_arr[ienergy] = R_side_rhofitmax_percent_dn;
        R_side_nevt_pct_up_arr[ienergy] = R_side_nevtavg_percent_up;
        R_side_nevt_pct_dn_arr[ienergy] = R_side_nevtavg_percent_dn;

        R_long_qlcms_pct_up_arr[ienergy] = R_long_qlcms_percent_up;
        R_long_qlcms_pct_dn_arr[ienergy] = R_long_qlcms_percent_dn;
        R_long_rhofit_pct_up_arr[ienergy] = R_long_rhofitmax_percent_up;
        R_long_rhofit_pct_dn_arr[ienergy] = R_long_rhofitmax_percent_dn;
        R_long_nevt_pct_up_arr[ienergy] = R_long_nevtavg_percent_up;
        R_long_nevt_pct_dn_arr[ienergy] = R_long_nevtavg_percent_dn;

        alpha_total_pct_up_arr[ienergy] = alpha_total_percent_up;
        alpha_total_pct_dn_arr[ienergy] = alpha_total_percent_dn;
        R_total_pct_up_arr[ienergy] = R_total_percent_up;
        R_total_pct_dn_arr[ienergy] = R_total_percent_dn;
        N_total_pct_up_arr[ienergy] = N_total_percent_up;
        N_total_pct_dn_arr[ienergy] = N_total_percent_dn;
        R_out_total_pct_up_arr[ienergy] = R_out_total_percent_up;
        R_out_total_pct_dn_arr[ienergy] = R_out_total_percent_dn;
        R_side_total_pct_up_arr[ienergy] = R_side_total_percent_up;
        R_side_total_pct_dn_arr[ienergy] = R_side_total_percent_dn;
        R_long_total_pct_up_arr[ienergy] = R_long_total_percent_up;
        R_long_total_pct_dn_arr[ienergy] = R_long_total_percent_dn;

        // mT-averaged defaults and ikt=2 values
        double sum_a=0., sum_R=0., sum_Nv=0., sum_Rout=0., sum_Rside=0., sum_Rlong=0.;
        int cnt=0;
        int cnt3d=0;
        for(int ikt=0; ikt<NKT; ikt++){
            double aya = alpha_default[ienergy]->GetY()[ikt];
            double rya = R_default[ienergy]->GetY()[ikt];
            double nya = N_default[ienergy]->GetY()[ikt];
            if(!std::isfinite(aya) || !std::isfinite(rya) || !std::isfinite(nya)) continue;
            sum_a += aya; sum_R += rya; sum_Nv += nya;
            if(R_out_default[ienergy] && R_side_default[ienergy] && R_long_default[ienergy]){
                double ro = R_out_default[ienergy]->GetY()[ikt];
                double rs = R_side_default[ienergy]->GetY()[ikt];
                double rl = R_long_default[ienergy]->GetY()[ikt];
                if(std::isfinite(ro) && std::isfinite(rs) && std::isfinite(rl)){
                    sum_Rout += ro;
                    sum_Rside += rs;
                    sum_Rlong += rl;
                    cnt3d++;
                }
            }
            cnt++;
        }
        if(cnt>0){
            avg_alpha[ienergy] = sum_a / cnt;
            avg_R[ienergy] = sum_R / cnt;
            avg_N[ienergy] = sum_Nv / cnt;
        }
        if(cnt3d>0){
            avg_R_out[ienergy] = sum_Rout / cnt3d;
            avg_R_side[ienergy] = sum_Rside / cnt3d;
            avg_R_long[ienergy] = sum_Rlong / cnt3d;
        }
        // store ikt=2 (3rd mT bin) default values if available
        if(NKT>2){
            val_alpha_ikt2[ienergy] = alpha_default[ienergy]->GetY()[2];
            val_R_ikt2[ienergy] = R_default[ienergy]->GetY()[2];
            val_N_ikt2[ienergy] = N_default[ienergy]->GetY()[2];
            if(R_out_default[ienergy] && R_side_default[ienergy] && R_long_default[ienergy]){
                val_R_out_ikt2[ienergy] = R_out_default[ienergy]->GetY()[2];
                val_R_side_ikt2[ienergy] = R_side_default[ienergy]->GetY()[2];
                val_R_long_ikt2[ienergy] = R_long_default[ienergy]->GetY()[2];
            }
        }

        // Populate m_T-choice source arrays in active code path for optional overlay usage.
        if(mtchoice_syst && NKT>2)
        {
            calc_neighbor_uncert(alpha_default[ienergy], 2, alpha_mTchoice_syst_up[ienergy], alpha_mTchoice_syst_dn[ienergy]);
            calc_neighbor_uncert(R_default[ienergy], 2, R_mTchoice_syst_up[ienergy], R_mTchoice_syst_dn[ienergy]);
            calc_neighbor_uncert(N_default[ienergy], 2, N_mTchoice_syst_up[ienergy], N_mTchoice_syst_dn[ienergy]);
            if(R_out_default[ienergy]) calc_neighbor_uncert(R_out_default[ienergy], 2, R_out_mTchoice_syst_up[ienergy], R_out_mTchoice_syst_dn[ienergy]);
            if(R_side_default[ienergy]) calc_neighbor_uncert(R_side_default[ienergy], 2, R_side_mTchoice_syst_up[ienergy], R_side_mTchoice_syst_dn[ienergy]);
            if(R_long_default[ienergy]) calc_neighbor_uncert(R_long_default[ienergy], 2, R_long_mTchoice_syst_up[ienergy], R_long_mTchoice_syst_dn[ienergy]);
        }
        
         
    } // end energy loop

    // NOTE: From now on, even though collected, I won't deal with "statistical" uncertainties here (StdDev)
    
    ////////////////////////////////////////////////////
    /////////////  PLOTTING   /////////////////////////
    ////////////////////////////////////////////////////
    
    
    ////////////////////////////////////////////////////
    ////////  PLOTTING mT vs parameter with sys bands (single panel)  
    ////////////////////////////////////////////////////
    // One combined panel per parameter across all energies.

    // Parameter loop: 0=alpha,1=R,2=N
    // ensure output dir exists
    gSystem->mkdir("figs/syserr", kTRUE);

    // Publication-oriented plot typography and output helpers.
    const double kAxisTitleSize = 0.060;
    const double kAxisLabelSize = 0.048;
    const double kLegendTextSizeMain = 0.034;
    const double kLegendTextSizeOverlay = 0.042;
    const double kTitleTextSize = 0.045;
    const double kAnnotTextSize = 0.043;

    auto style_frame_axes = [&](TH1* frame)
    {
        if(!frame) return;
        frame->GetXaxis()->SetTitleSize(kAxisTitleSize);
        frame->GetYaxis()->SetTitleSize(kAxisTitleSize);
        frame->GetXaxis()->SetLabelSize(kAxisLabelSize);
        frame->GetYaxis()->SetLabelSize(kAxisLabelSize);
        frame->GetXaxis()->SetTitleOffset(1.20);
        frame->GetYaxis()->SetTitleOffset(1.25);
    };

    // Convert coordinates defined in the inner frame [0,1]x[0,1] to pad NDC.
    auto frame_to_ndc_x = [&] (double x_rel)
    {
        return gPad->GetLeftMargin() + x_rel * (1.0 - gPad->GetLeftMargin() - gPad->GetRightMargin());
    };
    auto frame_to_ndc_y = [&] (double y_rel)
    {
        return gPad->GetBottomMargin() + y_rel * (1.0 - gPad->GetBottomMargin() - gPad->GetTopMargin());
    };

    auto save_canvas_publication = [&](TCanvas* can, const std::string& basepath_noext)
    {
        if(!can) return;
        can->SaveAs((basepath_noext + ".pdf").c_str());
        can->SaveAs((basepath_noext + ".png").c_str());
    };

    std::ofstream graph_txt("syserr_graphs.txt");
    if(!graph_txt.is_open())
    {
        cerr << "Warning: could not open syserr_graphs.txt for writing" << endl;
    }

    auto write_graph_table = [&](std::ofstream& out, TGraphAsymmErrors* g, const std::string& label, const std::string& x_label)
    {
        if(!g || !out.is_open()) return;
        out << label << "\n";
        out << x_label << " y err_up err_dn\n";
        int npts = g->GetN();
        for(int ipt=0; ipt<npts; ipt++)
        {
            double x = 0.0, y = 0.0;
            g->GetPoint(ipt, x, y);
            double err_up = g->GetEYhigh()[ipt];
            double err_dn = g->GetEYlow()[ipt];
            out << std::fixed << std::setprecision(6)
                << x << " " << y << " " << err_up << " " << err_dn << "\n";
        }
        out << "\n";
    };

    // Manual legend placement per panel to avoid covering error bars.
    const double leg_x1_fixed[3] = {0.28, 0.36, 0.26}; // alpha: BR, R: UR, N(lambda): bottom-middle (frame-relative)
    const double leg_y1_fixed[3] = {0.06, 0.65, 0.06};
    const double leg_x2_fixed[3] = {0.99, 0.99, 0.92};

    for(int iparam=0; iparam<3; iparam++)
    {
        TCanvas* can = new TCanvas(Form("can_param_%d", iparam), "", 1600, 1200);
        gPad->SetLeftMargin(0.18);
        gPad->SetBottomMargin(0.18);
        gPad->SetRightMargin(0.04);
        gPad->SetTopMargin(0.07);

        double ymins[] = {1.3, 3, 0.5};
        double ymaxs[] = {2.1, 9., 1.2};
        double ymin = ymins[iparam];
        double ymax = ymaxs[iparam];

        // Draw an empty frame with axis titles
        double xmin = mtbin_centers[0]-0.05;
        double xmax = mtbin_centers[NKT-1]+0.05;
        TH1F* frame = gPad->DrawFrame(xmin, ymin, xmax, ymax);
        const char* ytitle = (iparam==0)?"#alpha":(iparam==1)?"#LT#kern[-0.1]{R}#GT [fm]":"#lambda*"; // N
        frame->GetXaxis()->SetTitle("m_{T} [GeV/c^{2}]");
        frame->GetYaxis()->SetTitle(ytitle);
        style_frame_axes(frame);

        int nEntries = 0;
        for(int ie=0; ie<NENERGIES; ie++){
            if(energy_to_plot>-1 && ie!=energy_to_plot) continue;
            TGraphAsymmErrors* gdef = (iparam==0)? alpha_default[ie] : (iparam==1)? R_default[ie] : N_default[ie];
            TGraphAsymmErrors* gsys = (iparam==0)? alpha_syserr[ie] : (iparam==1)? R_syserr[ie] : N_syserr[ie];
            if(!gdef || !gsys) continue;
            nEntries++;
        }
        int nCols = (nEntries>14)? 3 : (nEntries>7)? 2 : 1;
        int nRows = (nEntries==0)? 0 : ( (nEntries + nCols - 1) / nCols );
        double lineH = 0.038;
        double vPad = 0.012;
        double legH = std::max(0.06, nRows*lineH + 2*vPad);
        double x1_rel = leg_x1_fixed[iparam] - ((nCols==1)? 0.0 : (nCols==2)? 0.05 : 0.10);
        double x2_rel = leg_x2_fixed[iparam];
        double y1_rel = leg_y1_fixed[iparam];
        x1_rel = std::max(0.02, x1_rel);
        x2_rel = std::min(0.98, x2_rel);

        double x1 = frame_to_ndc_x(x1_rel);
        double x2 = frame_to_ndc_x(x2_rel);
        double y1 = frame_to_ndc_y(y1_rel);
        double y2 = y1 + legH;
        double maxY2 = frame_to_ndc_y(0.98);
        if(y2 > maxY2)
        {
            y2 = maxY2;
            y1 = std::max(frame_to_ndc_y(0.02), y2 - legH);
        }
        TLegend* leg = new TLegend(x1, y1, x2, y2);
        leg->SetNColumns(nCols);
        leg->SetFillColor(0);
        leg->SetBorderSize(1);
        leg->SetTextSize(kLegendTextSizeMain-0.005);
        leg->SetMargin(0.20);
        leg->SetEntrySeparation(0.005);
        leg->SetColumnSeparation(0.02);

        // Colors
        int colors[] = {kBlack,kBlue+2,kRed+1,kGreen+2,kMagenta+2,kOrange+7,kViolet+1,kCyan+1};

        int colorIdx=0;
        for(int ie=0; ie<NENERGIES; ie++)
        {
            if(energy_to_plot>-1 && ie!=energy_to_plot) continue; // plot only selected energy if specified
            TGraphAsymmErrors* gdef = (iparam==0)? alpha_default[ie] : (iparam==1)? R_default[ie] : N_default[ie];
            TGraphAsymmErrors* gsys = (iparam==0)? alpha_syserr[ie] : (iparam==1)? R_syserr[ie] : N_syserr[ie];
            if(!gdef || !gsys) continue;
            int col = colors[colorIdx % (sizeof(colors)/sizeof(int))];
            // filled band: use sys graph draw option 3
            int fillcol = TColor::GetColorTransparent(col, 0.25); // 25-35% transparent fill of line color
            gsys->SetFillColor(fillcol);
            gsys->SetFillStyle(1001);
            gsys->SetLineColor(col);
            gsys->SetLineWidth(1);
            gsys->Draw("3 same");

            gsys->SetLineColor(col);
            gsys->SetLineWidth(2);
            gsys->SetMarkerStyle(20+colorIdx); // differentiate points
            gsys->SetMarkerSize(2.0); // instead of 0.0, points needed for better visibility in case of mT vs param
            gsys->SetMarkerColor(col);
            gsys->Draw("LPX same");

            write_graph_table(graph_txt, gsys, Form("mT_vs_%s_syserr_%s", levy_params[iparam], energies[ie]), "mT_GeV");

            std::stringstream label;
            label << std::fixed << std::setprecision(1) << energydouble[ie];
            leg->AddEntry(gsys, Form("#sqrt{s_{NN}} = %s GeV",label.str().c_str()), "lpx"); // changed from "l", let the marker show in legend
            colorIdx++;
        }
        leg->Draw();
        {
            // Draw title
            TLatex title;
            title.SetNDC();
            title.SetTextSize(kTitleTextSize);
            title.DrawLatex(frame_to_ndc_x(0.03), frame_to_ndc_y(0.95), "UrQMD Au+Au 0#minus10%, #kern[-0.3]{#pi^{#pm}#pi^{#pm}}, L#acute{e}vy fit");
        }

        // Save canvas (single-panel figure)
        const char* energyplotted = (energy_to_plot>-1)? Form("_energy_%s", energies[energy_to_plot]) : "_allenergies";
        can->SetTitle(Form("m_{T} vs %s systematic uncertainties%s",
                            (iparam==0)?"#alpha":(iparam==1)?"#LT#kern[-0.1]{R}#GT":"#lambda*", energyplotted));
        save_canvas_publication(can, Form("figs/syserr/mT_vs_param_%s%s", levy_params[iparam], energyplotted));
        delete can;
    }

    // Additional 3D-radius panels in one 3-panel canvas: R_out, R_side, R_long
    TCanvas* can3d = new TCanvas("can_param_Rosl", "", 2400, 800);
    TLegend* leg3d = nullptr;
    TPad* pads3d[3] = {nullptr, nullptr, nullptr};
    const double pad_x1[3] = {0.0, 1.0/3.0, 2.0/3.0};
    const double pad_x2[3] = {1.0/3.0, 2.0/3.0, 1.0};
    const double pad_left_margins[3] = {0.18, 0.0, 0.0};
    const double pad_right_margins[3] = {0.0, 0.0, 0.06};
    const double pad_bottom_margin = 0.18;
    const double pad_top_margin = 0.07;

    for(int ip=0; ip<3; ip++)
    {
        pads3d[ip] = new TPad(Form("pad_Rosl_%d", ip), "", pad_x1[ip], 0.0, pad_x2[ip], 1.0);
        pads3d[ip]->SetLeftMargin(pad_left_margins[ip]);
        pads3d[ip]->SetRightMargin(pad_right_margins[ip]);
        pads3d[ip]->SetBottomMargin(pad_bottom_margin);
        pads3d[ip]->SetTopMargin(pad_top_margin);
        pads3d[ip]->Draw();
    }

    int nEntries3d = 0;
    for(int ie=0; ie<NENERGIES; ie++)
    {
        if(energy_to_plot>-1 && ie!=energy_to_plot) continue;
        // Use middle panel availability to define shared legend content.
        TGraphAsymmErrors* gdef = R_side_default[ie];
        TGraphAsymmErrors* gsys = R_side_syserr[ie];
        if(!gdef || !gsys) continue;
        nEntries3d++;
    }
    int nCols3d = (nEntries3d>14)? 3 : (nEntries3d>7)? 2 : 1;

    for(int iparam=0; iparam<3; iparam++)
    {
        pads3d[iparam]->cd();

        double ymins[] = {2.0, 2.0, 2.0};
        double ymaxs[] = {11.0, 11.0, 11.0};
        double ymin = ymins[iparam];
        double ymax = ymaxs[iparam];

        double xmin = mtbin_centers[0]-0.05;
        double xmax = mtbin_centers[NKT-1]+0.05;
        TH1F* frame = gPad->DrawFrame(xmin, ymin, xmax, ymax);
        frame->GetXaxis()->SetTitle("m_{T} [GeV/c^{2}]");
        frame->GetYaxis()->SetTitle((iparam==0) ? "R_{out,side,long} [fm]" : "");
        style_frame_axes(frame);
        if(iparam>0)
        {
            frame->GetYaxis()->SetLabelSize(0.0);
            frame->GetYaxis()->SetTitleSize(0.0);
        }

        if(iparam==0)
        {
            // Title on the first panel only
            TLatex title;
            title.SetNDC();
            title.SetTextSize(kTitleTextSize);
            title.DrawLatex(frame_to_ndc_x(0.03), frame_to_ndc_y(0.95), "UrQMD Au+Au 0#minus10%, #kern[-0.3]{#pi^{#pm}#pi^{#pm}}, L#acute{e}vy fit");
        }
        if(iparam==1)
        {
            // Shared legend on middle panel; make it wider and place slightly lower.
            double lx1 = frame_to_ndc_x(0.08);
            double ly1 = frame_to_ndc_y(0.46);
            double lx2 = frame_to_ndc_x(0.98);
            double ly2 = frame_to_ndc_y(0.94);
            leg3d = new TLegend(lx1, ly1, lx2, ly2);
            leg3d->SetNColumns(nCols3d);
            leg3d->SetFillColor(0);
            leg3d->SetBorderSize(1);
            leg3d->SetTextSize(kLegendTextSizeMain);
            leg3d->SetMargin(0.20);
            leg3d->SetEntrySeparation(0.005);
            leg3d->SetColumnSeparation(0.02);
        }

        int colors[] = {kBlack,kBlue+2,kRed+1,kGreen+2,kMagenta+2,kOrange+7,kViolet+1,kCyan+1};
        int colorIdx=0;
        for(int ie=0; ie<NENERGIES; ie++)
        {
            if(energy_to_plot>-1 && ie!=energy_to_plot) continue;
            TGraphAsymmErrors* gdef = (iparam==0)? R_out_default[ie] : (iparam==1)? R_side_default[ie] : R_long_default[ie];
            TGraphAsymmErrors* gsys = (iparam==0)? R_out_syserr[ie] : (iparam==1)? R_side_syserr[ie] : R_long_syserr[ie];
            if(!gdef || !gsys) continue;
            int col = colors[colorIdx % (sizeof(colors)/sizeof(int))];
            int fillcol = TColor::GetColorTransparent(col, 0.25);
            gsys->SetFillColor(fillcol);
            gsys->SetFillStyle(1001);
            gsys->SetLineColor(col);
            gsys->SetLineWidth(1);
            gsys->Draw("3 same");

            gsys->SetLineColor(col);
            gsys->SetLineWidth(1);
            gsys->SetMarkerStyle(20+colorIdx);
            gsys->SetMarkerSize(2.0);
            gsys->SetMarkerColor(col);
            gsys->Draw("LPX same");

            write_graph_table(graph_txt, gsys, Form("mT_vs_%s_syserr_%s", levy_params_3d[iparam], energies[ie]), "mT_GeV");

            std::stringstream label;
            label << std::fixed << std::setprecision(1) << energydouble[ie];
            if(iparam==1 && leg3d)
            {
                leg3d->AddEntry(gsys, Form("#sqrt{s_{NN}} = %s GeV",label.str().c_str()), "lpx");
            }
            colorIdx++;
        }

        if(iparam==1 && leg3d)
        {
            leg3d->Draw();
        }

        {
            TLatex panel_label;
            panel_label.SetNDC();
            panel_label.SetTextSize(kAnnotTextSize);
            const char* label = (iparam==0)?"out":(iparam==1)?"side":"long";
            panel_label.DrawLatex(frame_to_ndc_x(0.06), frame_to_ndc_y(0.10), label);
        }
    }
    const char* energyplotted_3d = (energy_to_plot>-1)? Form("_energy_%s", energies[energy_to_plot]) : "_allenergies";
    can3d->SetTitle(Form("m_{T} vs 3D radii with systematic uncertainties%s", energyplotted_3d));
    save_canvas_publication(can3d, Form("figs/syserr/mT_vs_param_R_osl%s", energyplotted_3d));
    delete can3d;

    
    /////////////////////////////////////////////////////
    // SYSTEMATICS TO sqrt(s_NN) vs Lévy parameter /////
    ///////////////////////////////////////////////////
    
    // Prepare x-axis (collision energies as numbers)
    double x_energy[NENERGIES];
    for(int ie=0; ie<NENERGIES; ie++) x_energy[ie] = energydouble[ie];
    double xerr_low_s[NENERGIES] = {0};
    double xerr_high_s[NENERGIES] = {0};

    

    ////////////////////////////////////////////////////////////
    // Single-panel overlays: each parameter                ////
    // showing optional mT-averaged, and selected mT bins   ////
    ////////////////////////////////////////////////////////////
    
    for(int iparam=0; iparam<3; iparam++){
        TCanvas* can_overlay = new TCanvas(Form("can_overlay_param_%d", iparam), "", 1200, 800);
        can_overlay->SetLogx(1); // !!! FIXME if not needed to compare with STAR data
        gPad->SetLeftMargin(0.18);
        gPad->SetBottomMargin(0.18);
        gPad->SetRightMargin(0.04);
        gPad->SetTopMargin(0.07);

        
        // Precompute y-range for both datasets
        double ymin_both = 1e9, ymax_both = -1e9;
        
        // mT-averaged range
        for(int ie=0; ie<NENERGIES; ie++){
            double avg = (iparam==0)? avg_alpha[ie] : (iparam==1)? avg_R[ie] : avg_N[ie];
            double up_sq = 0., dn_sq = 0.;
            if(iparam==0){
                up_sq += pow(alpha_qlcms_pct_up_arr[ie]/100.0 * avg, 2);
                up_sq += pow(alpha_rhofit_pct_up_arr[ie]/100.0 * avg, 2);
                up_sq += pow(alpha_nevt_pct_up_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(alpha_qlcms_pct_dn_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(alpha_rhofit_pct_dn_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(alpha_nevt_pct_dn_arr[ie]/100.0 * avg, 2);
                if(mtchoice_syst){
                    up_sq += pow(alpha_mTchoice_syst_up[ie], 2);
                    dn_sq += pow(alpha_mTchoice_syst_dn[ie], 2);
                }
            } else if(iparam==1){
                up_sq += pow(R_qlcms_pct_up_arr[ie]/100.0 * avg, 2);
                up_sq += pow(R_rhofit_pct_up_arr[ie]/100.0 * avg, 2);
                up_sq += pow(R_nevt_pct_up_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(R_qlcms_pct_dn_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(R_rhofit_pct_dn_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(R_nevt_pct_dn_arr[ie]/100.0 * avg, 2);
                if(mtchoice_syst){
                    up_sq += pow(R_mTchoice_syst_up[ie], 2);
                    dn_sq += pow(R_mTchoice_syst_dn[ie], 2);
                }
            } else {
                up_sq += pow(N_qlcms_pct_up_arr[ie]/100.0 * avg, 2);
                up_sq += pow(N_rhofit_pct_up_arr[ie]/100.0 * avg, 2);
                up_sq += pow(N_nevt_pct_up_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(N_qlcms_pct_dn_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(N_rhofit_pct_dn_arr[ie]/100.0 * avg, 2);
                dn_sq += pow(N_nevt_pct_dn_arr[ie]/100.0 * avg, 2);
                if(mtchoice_syst){
                    up_sq += pow(N_mTchoice_syst_up[ie], 2);
                    dn_sq += pow(N_mTchoice_syst_dn[ie], 2);
                }
            }
            double errup = sqrt(up_sq);
            double errdn = sqrt(dn_sq);
            if(std::isfinite(avg)){
                ymin_both = std::min(ymin_both, avg - errdn);
                ymax_both = std::max(ymax_both, avg + errup);
            }
        }
        
        // selected mT-bin ranges
        std::vector<int> valid_ikts_overlay;
        for(int iktsel=0; iktsel<NIKT_SNN_PLOTS; iktsel++){
            int ikt = ikt_indices_to_plot[iktsel];
            if(ikt < 0 || ikt >= NKT) continue;
            valid_ikts_overlay.push_back(ikt);
            for(int ie=0; ie<NENERGIES; ie++){
                TGraphAsymmErrors* gdef = (iparam==0)? alpha_default[ie] : (iparam==1)? R_default[ie] : N_default[ie];
                TGraphAsymmErrors* gsys = (iparam==0)? alpha_syserr[ie] : (iparam==1)? R_syserr[ie] : N_syserr[ie];
                if(!gdef || !gsys || gdef->GetN()<=ikt || gsys->GetN()<=ikt) continue;

                double val = gdef->GetY()[ikt];
                double base_up = gsys->GetEYhigh()[ikt];
                double base_dn = gsys->GetEYlow()[ikt];
                double neigh_up_sq = 0., neigh_dn_sq = 0.;
                double center = gdef->GetY()[ikt];
                for(int j : {ikt-1, ikt+1}){
                    if(j<0 || j>=gdef->GetN()) continue;
                    double nb = gdef->GetY()[j];
                    double diff = nb - center;
                    if(diff>=0) neigh_up_sq += diff*diff; else neigh_dn_sq += diff*diff;
                }
                double neigh_up = sqrt(neigh_up_sq);
                double neigh_dn = sqrt(neigh_dn_sq);
                double extra_up = mtchoice_syst ? neigh_up : 0.;
                double extra_dn = mtchoice_syst ? neigh_dn : 0.;
                double errup2 = sqrt(base_up*base_up + extra_up*extra_up);
                double errdn2 = sqrt(base_dn*base_dn + extra_dn*extra_dn);
                if(std::isfinite(val)){
                    ymin_both = std::min(ymin_both, val - errdn2);
                    ymax_both = std::max(ymax_both, val + errup2);
                }
            }
        }
        
        // For alpha, also account for analytic curve range
        if(iparam==0){
            double fxmin = x_energy[0];
            double fxmax = x_energy[NENERGIES-1];
            int Nsample = 500;
            for(int is=0; is<=Nsample; ++is){
                double x = fxmin + (fxmax - fxmin) * (double)is / (double)Nsample;
                double y = 0.85 + pow(x, -0.14);
                if(std::isfinite(y)){
                    ymin_both = std::min(ymin_both, y);
                    ymax_both = std::max(ymax_both, y);
                }
            }
        }
        
        if(ymin_both > ymax_both){ ymin_both = 0.; ymax_both = 1.; }
        double ypadmargin = 0.06 * (ymax_both - ymin_both);
        ymin_both -= ypadmargin;
        ymax_both += ypadmargin;
        
        // Draw frame
        double xmin = x_energy[0] - 1.0;
        double xmax = x_energy[NENERGIES-1] + 1.0;
        TH1F* frame_overlay = gPad->DrawFrame(xmin, ymin_both, xmax, ymax_both);
        if(iparam==2)
        {
            frame_overlay = gPad->DrawFrame(xmin, 0.7, xmax, 1.08);
        }
        else if(iparam==0)
        {
            //frame_overlay = gPad->DrawFrame(1.0, 0.1, 10000., 2.1); // FIXME if not needed to compare with STAR data
        }
        frame_overlay->GetXaxis()->SetTitle("#sqrt{s_{NN}} [GeV]");
        const char* ytitle_overlay = (iparam==0)?"#alpha":(iparam==1)?"#LT#kern[-0.1]{R}#GT [fm]":"#lambda*"; // =N
        frame_overlay->GetYaxis()->SetTitle(ytitle_overlay);
        style_frame_axes(frame_overlay);
        
        // Draw mT-averaged data
        TGraphAsymmErrors* g_avg_overlay;
        {
            double yavg_overlay[NENERGIES];
            double yavg_errup_overlay[NENERGIES];
            double yavg_errdn_overlay[NENERGIES];
            for(int ie=0; ie<NENERGIES; ie++){
                double avg = (iparam==0)? avg_alpha[ie] : (iparam==1)? avg_R[ie] : avg_N[ie];
                yavg_overlay[ie] = avg;
                double up_sq = 0., dn_sq = 0.;
                if(iparam==0){
                    up_sq += pow(alpha_qlcms_pct_up_arr[ie]/100.0 * avg, 2);
                    up_sq += pow(alpha_rhofit_pct_up_arr[ie]/100.0 * avg, 2);
                    up_sq += pow(alpha_nevt_pct_up_arr[ie]/100.0 * avg, 2);
                    dn_sq += pow(alpha_qlcms_pct_dn_arr[ie]/100.0 * avg, 2);
                    dn_sq += pow(alpha_rhofit_pct_dn_arr[ie]/100.0 * avg, 2);
                    dn_sq += pow(alpha_nevt_pct_dn_arr[ie]/100.0 * avg, 2);
                    // also add m_T choice sys uncert
                    if(mtchoice_syst){
                        up_sq += pow(alpha_mTchoice_syst_up[ie], 2);
                        dn_sq += pow(alpha_mTchoice_syst_dn[ie], 2);
                    }
                } else if(iparam==1){
                    up_sq += pow(R_qlcms_pct_up_arr[ie]/100.0 * avg, 2);
                    up_sq += pow(R_rhofit_pct_up_arr[ie]/100.0 * avg, 2);
                    up_sq += pow(R_nevt_pct_up_arr[ie]/100.0 * avg, 2);
                    dn_sq += pow(R_qlcms_pct_dn_arr[ie]/100.0 * avg, 2);
                    dn_sq += pow(R_rhofit_pct_dn_arr[ie]/100.0 * avg, 2);
                    dn_sq += pow(R_nevt_pct_dn_arr[ie]/100.0 * avg, 2);
                    if(mtchoice_syst){
                        up_sq += pow(R_mTchoice_syst_up[ie], 2);
                        dn_sq += pow(R_mTchoice_syst_dn[ie], 2);
                    }
                } else {
                    up_sq += pow(N_qlcms_pct_up_arr[ie]/100.0 * avg, 2);
                    up_sq += pow(N_rhofit_pct_up_arr[ie]/100.0 * avg, 2);
                    up_sq += pow(N_nevt_pct_up_arr[ie]/100.0 * avg, 2);
                    dn_sq += pow(N_qlcms_pct_dn_arr[ie]/100.0 * avg, 2);
                    dn_sq += pow(N_rhofit_pct_dn_arr[ie]/100.0 * avg, 2);
                    dn_sq += pow(N_nevt_pct_dn_arr[ie]/100.0 * avg, 2);
                    if(mtchoice_syst){
                        up_sq += pow(N_mTchoice_syst_up[ie], 2);
                        dn_sq += pow(N_mTchoice_syst_dn[ie], 2);
                    }
                }
                yavg_errup_overlay[ie] = sqrt(up_sq);
                yavg_errdn_overlay[ie] = sqrt(dn_sq);
            }
            g_avg_overlay = new TGraphAsymmErrors(NENERGIES, x_energy, yavg_overlay, xerr_low_s, xerr_high_s, yavg_errdn_overlay, yavg_errup_overlay);
            int col_avg = kBlue+2;
            int fillcol_avg = TColor::GetColorTransparent(col_avg, 0.35);
            g_avg_overlay->SetFillColor(fillcol_avg);
            g_avg_overlay->SetFillStyle(1001);
            g_avg_overlay->SetLineColor(col_avg);
            g_avg_overlay->SetLineWidth(3);
            if(plot_mtavg) g_avg_overlay->Draw("3 same");
            g_avg_overlay->SetMarkerStyle(21);
            g_avg_overlay->SetMarkerSize(0.0);
            if(iparam==2)
            {
                //g_avg_overlay->GetYaxis()->SetRangeUser(0.9,1.1);
            }
            if(plot_mtavg) g_avg_overlay->Draw("LPX same"); // FIXME uncomment if averaged whenever needed
        }
        
        // Draw selected mT-bin data
        std::vector<TGraphAsymmErrors*> g_ikt_overlay;
        std::vector<std::string> g_ikt_overlay_labels;
        {
            int ikt_colors[] = {kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+7};
            int ikt_markers[] = {21, 22, 23, 33, 34};
            for(size_t ik=0; ik<valid_ikts_overlay.size(); ik++){
                int ikt = valid_ikts_overlay[ik];
                double y_overlay[NENERGIES] = {0};
                double yerrup_overlay[NENERGIES] = {0};
                double yerrdn_overlay[NENERGIES] = {0};
                for(int ie=0; ie<NENERGIES; ie++){
                    // FIXME? the syserr graph should have the same points as the default graph, only the errors differ - gdef could be omitted at this point if stat err not plotted
                    TGraphAsymmErrors* gdef = (iparam==0)? alpha_default[ie] : (iparam==1)? R_default[ie] : N_default[ie];
                    //TGraphAsymmErrors* gdef = (iparam==0)? alpha_syserr[ie] : (iparam==1)? R_syserr[ie] : N_syserr[ie]; // like that...
                    TGraphAsymmErrors* gsys = (iparam==0)? alpha_syserr[ie] : (iparam==1)? R_syserr[ie] : N_syserr[ie]; // still, can be left like this for sanity check
                    if(!gdef || !gsys || gdef->GetN()<=ikt || gsys->GetN()<=ikt) continue;
                    if(debug)
                    {
                        // Print out the values and errors for this point for debugging
                        double val = gdef->GetY()[ikt];
                        double errup = gsys->GetEYhigh()[ikt];
                        double errdn = gsys->GetEYlow()[ikt];
                        cerr << "Debug gdef/gsys: param=" << iparam << ", ie=" << ie << ", ikt=" << ikt << ", val=" << val << ", errup=" << errup << ", errdn=" << errdn << endl;
                    }

                    y_overlay[ie] = gdef->GetY()[ikt];
                    double base_up = gsys->GetEYhigh()[ikt];
                    double base_dn = gsys->GetEYlow()[ikt];
                    double neigh_up_sq = 0., neigh_dn_sq = 0.;
                    double center = gdef->GetY()[ikt];
                    for(int j : {ikt-1, ikt+1}){
                        if(j<0 || j>=gdef->GetN()) continue;
                        double nb = gdef->GetY()[j];
                        double diff = nb - center;
                        if(diff>=0) neigh_up_sq += diff*diff; else neigh_dn_sq += diff*diff;
                    }
                    double neigh_up = sqrt(neigh_up_sq);
                    double neigh_dn = sqrt(neigh_dn_sq);
                    double extra_up = mtchoice_syst ? neigh_up : 0.;
                    double extra_dn = mtchoice_syst ? neigh_dn : 0.;
                    yerrup_overlay[ie] = sqrt(base_up*base_up + extra_up*extra_up);
                    yerrdn_overlay[ie] = sqrt(base_dn*base_dn + extra_dn*extra_dn);
                }

                TGraphAsymmErrors* g_ikt = new TGraphAsymmErrors(NENERGIES, x_energy, y_overlay, xerr_low_s, xerr_high_s, yerrdn_overlay, yerrup_overlay);
                int col_ikt = ikt_colors[ik % (sizeof(ikt_colors)/sizeof(int))];
                int fillcol_ikt = TColor::GetColorTransparent(col_ikt, 0.35);
                g_ikt->SetFillColor(fillcol_ikt);
                g_ikt->SetFillStyle(1001);
                g_ikt->SetLineColor(col_ikt);
                g_ikt->SetLineWidth(3);
                g_ikt->Draw("3 same");
                g_ikt->SetMarkerStyle(ikt_markers[ik % (sizeof(ikt_markers)/sizeof(int))]);
                g_ikt->SetMarkerSize(0.0);
                g_ikt->SetMarkerColor(col_ikt);
                g_ikt->Draw("LPX same");

                g_ikt_overlay.push_back(g_ikt);
                double mt_edge_dn = sqrt(ktbins[ikt]*ktbins[ikt] + Mass2_pi); double mt_edge_up = sqrt(ktbins[ikt+1]*ktbins[ikt+1] + Mass2_pi);
                g_ikt_overlay_labels.push_back(Form("m_{T} = %.0f#minus%.0f MeV", mt_edge_dn*1000.0, mt_edge_up*1000.0));
                write_graph_table(graph_txt, g_ikt, Form("sqrtS_overlay_%s_ikt%d_mT=%.0f-%.0fMeV", levy_params[iparam], ikt, mt_edge_dn*1000.0, mt_edge_up*1000.0), "sqrtSNN_GeV");
            }
        }
        
        // For alpha, draw analytic curve
        TF1* f_analytic = nullptr;
        if(iparam==0){
            double fxmin = x_energy[0];
            double fxmax = x_energy[NENERGIES-1];
            f_analytic = new TF1("alpha_curve_overlay", "0.85 + pow(x, -0.14)", fxmin, fxmax);
            f_analytic->SetLineColor(kBlack);
            f_analytic->SetLineStyle(1);
            f_analytic->SetLineWidth(4);
            f_analytic->SetNpx(500);
            if(comparewithSTAR) f_analytic->Draw("L same");
        }
        
        // Legend
        double leg_x1, leg_x2, leg_y1, leg_y2;
        if(iparam==0){
            // bottom left
            leg_x1 = 0.03; leg_x2 = 0.56; leg_y1 = 0.06; leg_y2 = 0.42;
        } else if(iparam==1){
            leg_x1 = 0.12; leg_x2 = 0.46; leg_y1 = 0.60; leg_y2 = 0.90;
        } else {
            leg_x1 = 0.12; leg_x2 = 0.46; leg_y1 = 0.06; leg_y2 = 0.36;
        }
        TLegend* leg_overlay = new TLegend(frame_to_ndc_x(leg_x1), frame_to_ndc_y(leg_y1), frame_to_ndc_x(leg_x2), frame_to_ndc_y(leg_y2));
        leg_overlay->SetBorderSize(0);
        leg_overlay->SetFillStyle(0);
        leg_overlay->SetTextSize(kLegendTextSizeOverlay);
        if(plot_mtavg) leg_overlay->AddEntry(g_avg_overlay, mtchoice_syst ? "UrQMD <m_{T}>, incl. m_{T} choice sys.unc." : "UrQMD <m_{T}>", "f");
        for(size_t ig=0; ig<g_ikt_overlay.size(); ig++){
            leg_overlay->AddEntry(g_ikt_overlay[ig], g_ikt_overlay_labels[ig].c_str(), "f");
        }
        if(iparam==0){
            if(comparewithSTAR) leg_overlay->AddEntry(f_analytic, "#alpha=0.85 + #sqrt{s_{NN}}^{-0.14} (trend of STAR data)", "l");
        }
        leg_overlay->Draw("L same");

        // Plot title
        {
            TLatex title;
            title.SetNDC();
            title.SetTextSize(kTitleTextSize);
            title.DrawLatex(frame_to_ndc_x(0.03), frame_to_ndc_y(0.95), Form("UrQMD Au+Au 0#minus10%%, #kern[-0.3]{#pi^{#pm}#pi^{#pm}}, L#acute{e}vy fit%s", (energy_to_plot>-1)? Form(", #sqrt{s_{NN}} = %s GeV", energies[energy_to_plot]) : ""));
        }
        
        save_canvas_publication(can_overlay, Form("figs/syserr/sqrtS_overlay_%s", levy_params[iparam]));
        delete can_overlay;
    }

    // Transposed sqrt(sNN) overlays:
    // one panel per selected mT bin, each showing
    // R (from existing R_default/R_syserr at that bin), R_out, R_side, R_long.

    auto get_ikt_point_from_graphs = [&](TGraphAsymmErrors* gdef, TGraphAsymmErrors* gsys, int ikt,
                                         double& val, double& errup, double& errdn) -> bool
    {
        if(!gdef || !gsys) return false;
        if(gdef->GetN()<=ikt || gsys->GetN()<=ikt) return false;
        val = gdef->GetY()[ikt];
        if(!std::isfinite(val)) return false;

        double base_up = gsys->GetEYhigh()[ikt];
        double base_dn = gsys->GetEYlow()[ikt];
        double neigh_up_sq = 0., neigh_dn_sq = 0.;
        if(mtchoice_syst)
        {
            double center = gdef->GetY()[ikt];
            for(int j : {ikt-1, ikt+1})
            {
                if(j<0 || j>=gdef->GetN()) continue;
                double nb = gdef->GetY()[j];
                double diff = nb - center;
                if(diff>=0.) neigh_up_sq += diff*diff;
                else neigh_dn_sq += diff*diff;
            }
        }
        errup = sqrt(base_up*base_up + neigh_up_sq);
        errdn = sqrt(base_dn*base_dn + neigh_dn_sq);
        return true;
    };

    std::vector<int> valid_ikts_overlay;
    for(int iktsel=0; iktsel<NIKT_SNN_PLOTS; iktsel++)
    {
        int ikt = ikt_indices_to_plot[iktsel];
        if(ikt >= 0 && ikt < NKT) valid_ikts_overlay.push_back(ikt);
    }

    // Use a shared y-range for all transposed panels so they are directly comparable.
    double global_ymin = 1e9, global_ymax = -1e9;
    for(size_t ik=0; ik<valid_ikts_overlay.size(); ik++)
    {
        int ikt = valid_ikts_overlay[ik];
        for(int ie=0; ie<NENERGIES; ie++)
        {
            double val=0., up=0., dn=0.;
            if(get_ikt_point_from_graphs(R_default[ie], R_syserr[ie], ikt, val, up, dn))
            {
                global_ymin = std::min(global_ymin, val - dn);
                global_ymax = std::max(global_ymax, val + up);
            }
            if(get_ikt_point_from_graphs(R_out_default[ie], R_out_syserr[ie], ikt, val, up, dn))
            {
                global_ymin = std::min(global_ymin, val - dn);
                global_ymax = std::max(global_ymax, val + up);
            }
            if(get_ikt_point_from_graphs(R_side_default[ie], R_side_syserr[ie], ikt, val, up, dn))
            {
                global_ymin = std::min(global_ymin, val - dn);
                global_ymax = std::max(global_ymax, val + up);
            }
            if(get_ikt_point_from_graphs(R_long_default[ie], R_long_syserr[ie], ikt, val, up, dn))
            {
                global_ymin = std::min(global_ymin, val - dn);
                global_ymax = std::max(global_ymax, val + up);
            }
        }
    }
    if(global_ymin > global_ymax)
    {
        global_ymin = 0.;
        global_ymax = 1.;
    }
    double global_ypadmargin = 0.06 * (global_ymax - global_ymin);
    global_ymin -= global_ypadmargin;
    global_ymax += global_ypadmargin;

    for(size_t ik=0; ik<valid_ikts_overlay.size(); ik++)
    {
        int ikt = valid_ikts_overlay[ik];
        TCanvas* can_overlay = new TCanvas(Form("can_overlay_R_osl_transposed_ikt_%d", ikt), "", 1200, 800);
        can_overlay->SetLogx(1);
        gPad->SetLeftMargin(0.18);
        gPad->SetBottomMargin(0.18);
        gPad->SetRightMargin(0.04);
        gPad->SetTopMargin(0.07);

        double y_Ravg[NENERGIES] = {0}, y_Ravg_up[NENERGIES] = {0}, y_Ravg_dn[NENERGIES] = {0};
        double y_Rout[NENERGIES] = {0}, y_Rout_up[NENERGIES] = {0}, y_Rout_dn[NENERGIES] = {0};
        double y_Rside[NENERGIES] = {0}, y_Rside_up[NENERGIES] = {0}, y_Rside_dn[NENERGIES] = {0};
        double y_Rlong[NENERGIES] = {0}, y_Rlong_up[NENERGIES] = {0}, y_Rlong_dn[NENERGIES] = {0};

        double ymin_both = global_ymin, ymax_both = global_ymax;
        for(int ie=0; ie<NENERGIES; ie++)
        {
            double val=0., up=0., dn=0.;

            if(get_ikt_point_from_graphs(R_default[ie], R_syserr[ie], ikt, val, up, dn))
            {
                y_Ravg[ie] = val;
                y_Ravg_up[ie] = up;
                y_Ravg_dn[ie] = dn;
                ymin_both = std::min(ymin_both, val - dn);
                ymax_both = std::max(ymax_both, val + up);
            }

            if(get_ikt_point_from_graphs(R_out_default[ie], R_out_syserr[ie], ikt, val, up, dn))
            {
                y_Rout[ie] = val;
                y_Rout_up[ie] = up;
                y_Rout_dn[ie] = dn;
                ymin_both = std::min(ymin_both, val - dn);
                ymax_both = std::max(ymax_both, val + up);
            }

            if(get_ikt_point_from_graphs(R_side_default[ie], R_side_syserr[ie], ikt, val, up, dn))
            {
                y_Rside[ie] = val;
                y_Rside_up[ie] = up;
                y_Rside_dn[ie] = dn;
                ymin_both = std::min(ymin_both, val - dn);
                ymax_both = std::max(ymax_both, val + up);
            }

            if(get_ikt_point_from_graphs(R_long_default[ie], R_long_syserr[ie], ikt, val, up, dn))
            {
                y_Rlong[ie] = val;
                y_Rlong_up[ie] = up;
                y_Rlong_dn[ie] = dn;
                ymin_both = std::min(ymin_both, val - dn);
                ymax_both = std::max(ymax_both, val + up);
            }
        }

        double xmin = x_energy[0] - 1.0;
        double xmax = x_energy[NENERGIES-1] + 1.0;
        TH1F* frame_overlay = gPad->DrawFrame(xmin, ymin_both, xmax, ymax_both);
        frame_overlay->GetXaxis()->SetTitle("#sqrt{s_{NN}} [GeV]");
        frame_overlay->GetYaxis()->SetTitle("R [fm]");
        style_frame_axes(frame_overlay);

        int cols[] = {kBlack, kRed+1, kBlue+1, kGreen+2};
        const char* labels[] = {"#LT#kern[-0.1]{R}#GT = #sqrt{(R_{out}^{2} + R_{side}^{2} + R_{long}^{2})/3}", "R_{out}", "R_{side}", "R_{long}"};
        int markers[] = {20, 21, 22, 23};

        TGraphAsymmErrors* g_Ravg = new TGraphAsymmErrors(NENERGIES, x_energy, y_Ravg, xerr_low_s, xerr_high_s, y_Ravg_dn, y_Ravg_up);
        TGraphAsymmErrors* g_Rout = new TGraphAsymmErrors(NENERGIES, x_energy, y_Rout, xerr_low_s, xerr_high_s, y_Rout_dn, y_Rout_up);
        TGraphAsymmErrors* g_Rside = new TGraphAsymmErrors(NENERGIES, x_energy, y_Rside, xerr_low_s, xerr_high_s, y_Rside_dn, y_Rside_up);
        TGraphAsymmErrors* g_Rlong = new TGraphAsymmErrors(NENERGIES, x_energy, y_Rlong, xerr_low_s, xerr_high_s, y_Rlong_dn, y_Rlong_up);
        TGraphAsymmErrors* graphs[] = {g_Ravg, g_Rout, g_Rside, g_Rlong};

        write_graph_table(graph_txt, g_Ravg, Form("sqrtS_transposed_ikt%d_Ravg", ikt), "sqrtSNN_GeV");
        write_graph_table(graph_txt, g_Rout, Form("sqrtS_transposed_ikt%d_Rout", ikt), "sqrtSNN_GeV");
        write_graph_table(graph_txt, g_Rside, Form("sqrtS_transposed_ikt%d_Rside", ikt), "sqrtSNN_GeV");
        write_graph_table(graph_txt, g_Rlong, Form("sqrtS_transposed_ikt%d_Rlong", ikt), "sqrtSNN_GeV");

        for(int ig=0; ig<4; ig++)
        {
            int fillcol = TColor::GetColorTransparent(cols[ig], 0.28);
            graphs[ig]->SetFillColor(fillcol);
            graphs[ig]->SetFillStyle(1001);
            graphs[ig]->SetLineColor(cols[ig]);
            graphs[ig]->SetLineWidth(3);
            graphs[ig]->Draw("3 same");
            graphs[ig]->SetMarkerStyle(markers[ig]);
            graphs[ig]->SetMarkerSize(0.0);
            graphs[ig]->SetMarkerColor(cols[ig]);
            graphs[ig]->Draw("LPX same");
        }

        TLegend* leg_overlay = new TLegend(frame_to_ndc_x(0.03), frame_to_ndc_y(0.58), frame_to_ndc_x(0.42), frame_to_ndc_y(0.91));
        leg_overlay->SetBorderSize(0);
        leg_overlay->SetFillStyle(0);
        leg_overlay->SetTextSize(kLegendTextSizeOverlay);
        for(int ig=0; ig<4; ig++)
        {
            leg_overlay->AddEntry(graphs[ig], labels[ig], "f");
        }
        leg_overlay->Draw("same");

        // Plot title
        {
            TLatex lat;
            lat.SetNDC();
            lat.SetTextSize(kAnnotTextSize);
            double mt_edge_dn = sqrt(ktbins[ikt]*ktbins[ikt] + Mass2_pi); double mt_edge_up = sqrt(ktbins[ikt+1]*ktbins[ikt+1] + Mass2_pi);
            lat.DrawLatex(frame_to_ndc_x(0.03), frame_to_ndc_y(0.95), Form("UrQMD Au+Au, #kern[-0.3]{#pi^{#pm}#pi^{#pm}}, m_{T} = %.0f#minus%.0f MeV", mt_edge_dn*1000.0, mt_edge_up*1000.0));
        }

        save_canvas_publication(can_overlay, Form("figs/syserr/sqrtS_overlay_R_osl_and_Ravg_ikt%d", ikt));
        delete can_overlay;
    }

    // Derived sqrt(sNN) overlays from 3D radii:
    //   0 -> R_out / R_side
    //   1 -> R_out^{2} - R_side^{2}
    auto build_total_err_for_3davg = [&](int ie, double Rout, double Rside,
                                         double& dRout_up, double& dRout_dn,
                                         double& dRside_up, double& dRside_dn)
    {
        dRout_up = 0.; dRout_dn = 0.;
        dRside_up = 0.; dRside_dn = 0.;

        double up_sq_out = 0., dn_sq_out = 0.;
        up_sq_out += pow(R_out_qlcms_pct_up_arr[ie]/100.0 * Rout, 2);
        up_sq_out += pow(R_out_rhofit_pct_up_arr[ie]/100.0 * Rout, 2);
        up_sq_out += pow(R_out_nevt_pct_up_arr[ie]/100.0 * Rout, 2);
        dn_sq_out += pow(R_out_qlcms_pct_dn_arr[ie]/100.0 * Rout, 2);
        dn_sq_out += pow(R_out_rhofit_pct_dn_arr[ie]/100.0 * Rout, 2);
        dn_sq_out += pow(R_out_nevt_pct_dn_arr[ie]/100.0 * Rout, 2);
        if(mtchoice_syst){
            up_sq_out += pow(R_out_mTchoice_syst_up[ie], 2);
            dn_sq_out += pow(R_out_mTchoice_syst_dn[ie], 2);
        }

        double up_sq_side = 0., dn_sq_side = 0.;
        up_sq_side += pow(R_side_qlcms_pct_up_arr[ie]/100.0 * Rside, 2);
        up_sq_side += pow(R_side_rhofit_pct_up_arr[ie]/100.0 * Rside, 2);
        up_sq_side += pow(R_side_nevt_pct_up_arr[ie]/100.0 * Rside, 2);
        dn_sq_side += pow(R_side_qlcms_pct_dn_arr[ie]/100.0 * Rside, 2);
        dn_sq_side += pow(R_side_rhofit_pct_dn_arr[ie]/100.0 * Rside, 2);
        dn_sq_side += pow(R_side_nevt_pct_dn_arr[ie]/100.0 * Rside, 2);
        if(mtchoice_syst){
            up_sq_side += pow(R_side_mTchoice_syst_up[ie], 2);
            dn_sq_side += pow(R_side_mTchoice_syst_dn[ie], 2);
        }

        dRout_up = sqrt(up_sq_out);
        dRout_dn = sqrt(dn_sq_out);
        dRside_up = sqrt(up_sq_side);
        dRside_dn = sqrt(dn_sq_side);
    };

    auto build_total_err_for_3d_ikt = [&](int ie, int ikt,
                                          double& Rout, double& Rside,
                                          double& dRout_up, double& dRout_dn,
                                          double& dRside_up, double& dRside_dn) -> bool
    {
        TGraphAsymmErrors* gout_def = R_out_default[ie];
        TGraphAsymmErrors* gside_def = R_side_default[ie];
        TGraphAsymmErrors* gout_sys = R_out_syserr[ie];
        TGraphAsymmErrors* gside_sys = R_side_syserr[ie];
        if(!gout_def || !gside_def || !gout_sys || !gside_sys) return false;
        if(gout_def->GetN()<=ikt || gside_def->GetN()<=ikt || gout_sys->GetN()<=ikt || gside_sys->GetN()<=ikt) return false;

        Rout = gout_def->GetY()[ikt];
        Rside = gside_def->GetY()[ikt];

        double base_out_up = gout_sys->GetEYhigh()[ikt];
        double base_out_dn = gout_sys->GetEYlow()[ikt];
        double base_side_up = gside_sys->GetEYhigh()[ikt];
        double base_side_dn = gside_sys->GetEYlow()[ikt];

        double neigh_out_up_sq = 0., neigh_out_dn_sq = 0.;
        double neigh_side_up_sq = 0., neigh_side_dn_sq = 0.;
        if(mtchoice_syst)
        {
            for(int j : {ikt-1, ikt+1})
            {
                if(j>=0 && j<gout_def->GetN())
                {
                    double diff = gout_def->GetY()[j] - Rout;
                    if(diff>=0.) neigh_out_up_sq += diff*diff;
                    else neigh_out_dn_sq += diff*diff;
                }
                if(j>=0 && j<gside_def->GetN())
                {
                    double diff = gside_def->GetY()[j] - Rside;
                    if(diff>=0.) neigh_side_up_sq += diff*diff;
                    else neigh_side_dn_sq += diff*diff;
                }
            }
        }

        dRout_up = sqrt(base_out_up*base_out_up + neigh_out_up_sq);
        dRout_dn = sqrt(base_out_dn*base_out_dn + neigh_out_dn_sq);
        dRside_up = sqrt(base_side_up*base_side_up + neigh_side_up_sq);
        dRside_dn = sqrt(base_side_dn*base_side_dn + neigh_side_dn_sq);
        return std::isfinite(Rout) && std::isfinite(Rside);
    };

    auto propagate_ratio_err = [](double Rout, double Rside, double dRout, double dRside) -> double
    {
        if(!std::isfinite(Rout) || !std::isfinite(Rside) || !std::isfinite(dRout) || !std::isfinite(dRside)) return 0.;
        if(fabs(Rout) < 1e-12 || fabs(Rside) < 1e-12) return 0.;
        double ratio = Rout / Rside;
        return fabs(ratio) * sqrt(pow(dRout / Rout, 2) + pow(dRside / Rside, 2));
    };

    auto propagate_diff2_err = [](double Rout, double Rside, double dRout, double dRside) -> double
    {
        if(!std::isfinite(Rout) || !std::isfinite(Rside) || !std::isfinite(dRout) || !std::isfinite(dRside)) return 0.;
        return sqrt(pow(2.0 * Rout * dRout, 2) + pow(2.0 * Rside * dRside, 2));
    };

    for(int ider=0; ider<2; ider++)
    {
        TCanvas* can_overlay = new TCanvas(Form("can_overlay_derived_3d_%d", ider), "", 1200, 800);
        can_overlay->SetLogx(1);
        gPad->SetLeftMargin(0.18);
        gPad->SetBottomMargin(0.18);
        gPad->SetRightMargin(0.04);
        gPad->SetTopMargin(0.07);

        double ymin_both = 1e9, ymax_both = -1e9;

        // Optional mT-averaged trend
        for(int ie=0; ie<NENERGIES; ie++)
        {
            double Rout = avg_R_out[ie];
            double Rside = avg_R_side[ie];
            if(!std::isfinite(Rout) || !std::isfinite(Rside) || fabs(Rside) < 1e-12) continue;

            double dRout_up=0., dRout_dn=0., dRside_up=0., dRside_dn=0.;
            build_total_err_for_3davg(ie, Rout, Rside, dRout_up, dRout_dn, dRside_up, dRside_dn);

            double y = (ider==0) ? (Rout / Rside) : (Rout*Rout - Rside*Rside);
            double errup = (ider==0) ? propagate_ratio_err(Rout, Rside, dRout_up, dRside_up)
                                     : propagate_diff2_err(Rout, Rside, dRout_up, dRside_up);
            double errdn = (ider==0) ? propagate_ratio_err(Rout, Rside, dRout_dn, dRside_dn)
                                     : propagate_diff2_err(Rout, Rside, dRout_dn, dRside_dn);
            if(std::isfinite(y))
            {
                ymin_both = std::min(ymin_both, y - errdn);
                ymax_both = std::max(ymax_both, y + errup);
            }
        }

        std::vector<int> valid_ikts_overlay;
        for(int iktsel=0; iktsel<NIKT_SNN_PLOTS; iktsel++)
        {
            int ikt = ikt_indices_to_plot[iktsel];
            if(ikt < 0 || ikt >= NKT) continue;
            valid_ikts_overlay.push_back(ikt);

            for(int ie=0; ie<NENERGIES; ie++)
            {
                double Rout=0., Rside=0.;
                double dRout_up=0., dRout_dn=0., dRside_up=0., dRside_dn=0.;
                if(!build_total_err_for_3d_ikt(ie, ikt, Rout, Rside, dRout_up, dRout_dn, dRside_up, dRside_dn)) continue;
                if(fabs(Rside) < 1e-12) continue;

                double y = (ider==0) ? (Rout / Rside) : (Rout*Rout - Rside*Rside);
                double errup = (ider==0) ? propagate_ratio_err(Rout, Rside, dRout_up, dRside_up)
                                         : propagate_diff2_err(Rout, Rside, dRout_up, dRside_up);
                double errdn = (ider==0) ? propagate_ratio_err(Rout, Rside, dRout_dn, dRside_dn)
                                         : propagate_diff2_err(Rout, Rside, dRout_dn, dRside_dn);

                if(std::isfinite(y))
                {
                    ymin_both = std::min(ymin_both, y - errdn);
                    ymax_both = std::max(ymax_both, y + errup);
                }
            }
        }

        if(ymin_both > ymax_both)
        {
            ymin_both = (ider==0) ? 0.8 : -5.0;
            ymax_both = (ider==0) ? 1.8 : 15.0;
        }
        double ypadmargin = 0.06 * (ymax_both - ymin_both);
        ymin_both -= ypadmargin;
        ymax_both += ypadmargin;

        double xmin = x_energy[0] - 1.0;
        double xmax = x_energy[NENERGIES-1] + 1.0;
        TH1F* frame_overlay = gPad->DrawFrame(xmin, ymin_both, xmax, ymax_both);
        frame_overlay->GetXaxis()->SetTitle("#sqrt{s_{NN}} [GeV]");
        const char* ytitle_overlay = (ider==0) ? "R_{out}/R_{side}" : "R_{out}^{2}-R_{side}^{2} [fm^{2}]";
        frame_overlay->GetYaxis()->SetTitle(ytitle_overlay);
        style_frame_axes(frame_overlay);

        TGraphAsymmErrors* g_avg_overlay = nullptr;
        {
            double yavg_overlay[NENERGIES] = {0};
            double yavg_errup_overlay[NENERGIES] = {0};
            double yavg_errdn_overlay[NENERGIES] = {0};
            for(int ie=0; ie<NENERGIES; ie++)
            {
                double Rout = avg_R_out[ie];
                double Rside = avg_R_side[ie];
                if(!std::isfinite(Rout) || !std::isfinite(Rside) || fabs(Rside) < 1e-12) continue;

                double dRout_up=0., dRout_dn=0., dRside_up=0., dRside_dn=0.;
                build_total_err_for_3davg(ie, Rout, Rside, dRout_up, dRout_dn, dRside_up, dRside_dn);

                yavg_overlay[ie] = (ider==0) ? (Rout / Rside) : (Rout*Rout - Rside*Rside);
                yavg_errup_overlay[ie] = (ider==0) ? propagate_ratio_err(Rout, Rside, dRout_up, dRside_up)
                                                   : propagate_diff2_err(Rout, Rside, dRout_up, dRside_up);
                yavg_errdn_overlay[ie] = (ider==0) ? propagate_ratio_err(Rout, Rside, dRout_dn, dRside_dn)
                                                   : propagate_diff2_err(Rout, Rside, dRout_dn, dRside_dn);
            }

            g_avg_overlay = new TGraphAsymmErrors(NENERGIES, x_energy, yavg_overlay, xerr_low_s, xerr_high_s, yavg_errdn_overlay, yavg_errup_overlay);
            int col_avg = kBlue+2;
            int fillcol_avg = TColor::GetColorTransparent(col_avg, 0.35);
            g_avg_overlay->SetFillColor(fillcol_avg);
            g_avg_overlay->SetFillStyle(1001);
            g_avg_overlay->SetLineColor(col_avg);
            g_avg_overlay->SetLineWidth(3);
            if(plot_mtavg) g_avg_overlay->Draw("3 same");
            g_avg_overlay->SetMarkerStyle(21);
            g_avg_overlay->SetMarkerSize(0.0);
            if(plot_mtavg) g_avg_overlay->Draw("LPX same");
        }

        std::vector<TGraphAsymmErrors*> g_ikt_overlay;
        std::vector<std::string> g_ikt_overlay_labels;
        {
            int ikt_colors[] = {kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+7};
            int ikt_markers[] = {21, 22, 23, 33, 34};
            for(size_t ik=0; ik<valid_ikts_overlay.size(); ik++)
            {
                int ikt = valid_ikts_overlay[ik];
                double y_overlay[NENERGIES] = {0};
                double yerrup_overlay[NENERGIES] = {0};
                double yerrdn_overlay[NENERGIES] = {0};

                for(int ie=0; ie<NENERGIES; ie++)
                {
                    double Rout=0., Rside=0.;
                    double dRout_up=0., dRout_dn=0., dRside_up=0., dRside_dn=0.;
                    if(!build_total_err_for_3d_ikt(ie, ikt, Rout, Rside, dRout_up, dRout_dn, dRside_up, dRside_dn)) continue;
                    if(fabs(Rside) < 1e-12) continue;

                    y_overlay[ie] = (ider==0) ? (Rout / Rside) : (Rout*Rout - Rside*Rside);
                    yerrup_overlay[ie] = (ider==0) ? propagate_ratio_err(Rout, Rside, dRout_up, dRside_up)
                                                   : propagate_diff2_err(Rout, Rside, dRout_up, dRside_up);
                    yerrdn_overlay[ie] = (ider==0) ? propagate_ratio_err(Rout, Rside, dRout_dn, dRside_dn)
                                                   : propagate_diff2_err(Rout, Rside, dRout_dn, dRside_dn);
                }

                TGraphAsymmErrors* g_ikt = new TGraphAsymmErrors(NENERGIES, x_energy, y_overlay, xerr_low_s, xerr_high_s, yerrdn_overlay, yerrup_overlay);
                int col_ikt = ikt_colors[ik % (sizeof(ikt_colors)/sizeof(int))];
                int fillcol_ikt = TColor::GetColorTransparent(col_ikt, 0.35);
                g_ikt->SetFillColor(fillcol_ikt);
                g_ikt->SetFillStyle(1001);
                g_ikt->SetLineColor(col_ikt);
                g_ikt->SetLineWidth(3);
                g_ikt->Draw("3 same");
                g_ikt->SetMarkerStyle(ikt_markers[ik % (sizeof(ikt_markers)/sizeof(int))]);
                g_ikt->SetMarkerSize(0.0);
                g_ikt->SetMarkerColor(col_ikt);
                g_ikt->Draw("LPX same");

                g_ikt_overlay.push_back(g_ikt);
                double mt_edge_dn = sqrt(ktbins[ikt]*ktbins[ikt] + Mass2_pi); double mt_edge_up = sqrt(ktbins[ikt+1]*ktbins[ikt+1] + Mass2_pi);
                g_ikt_overlay_labels.push_back(Form("m_{T} = %.0f#minus%.0f MeV", mt_edge_dn*1000.0, mt_edge_up*1000.0));
                const char* der_label = (ider==0) ? "Rout_over_Rside" : "Rout2_minus_Rside2";
                write_graph_table(graph_txt, g_ikt, Form("sqrtS_derived_%s_ikt%d_mT=%.0f-%.0fMeV", der_label, ikt, mt_edge_dn*1000.0, mt_edge_up*1000.0), "sqrtSNN_GeV");
            }
        }

        TLegend* leg_overlay = new TLegend(frame_to_ndc_x(0.12), frame_to_ndc_y(0.60), frame_to_ndc_x(0.46), frame_to_ndc_y(0.90));
        leg_overlay->SetBorderSize(0);
        leg_overlay->SetFillStyle(0);
        leg_overlay->SetTextSize(kLegendTextSizeOverlay);
        if(plot_mtavg) leg_overlay->AddEntry(g_avg_overlay, mtchoice_syst ? "UrQMD <m_{T}>, incl. m_{T} choice sys.unc." : "UrQMD <m_{T}>", "f");
        for(size_t ig=0; ig<g_ikt_overlay.size(); ig++)
        {
            leg_overlay->AddEntry(g_ikt_overlay[ig], g_ikt_overlay_labels[ig].c_str(), "f");
        }
        leg_overlay->Draw("L same");
        // Add title
        {
            TLatex title;
            title.SetNDC();
            title.SetTextSize(kTitleTextSize);
            title.DrawLatex(frame_to_ndc_x(0.03), frame_to_ndc_y(0.95), Form("UrQMD Au+Au 0#minus10%%, #kern[-0.3]{#pi^{#pm}#pi^{#pm}}, L#acute{e}vy fit%s", (energy_to_plot>-1)? Form(", #sqrt{s_{NN}} = %s GeV", energies[energy_to_plot]) : ""));
        }

        if(ider==0) save_canvas_publication(can_overlay, "figs/syserr/sqrtS_overlay_Rout_over_Rside");
        else save_canvas_publication(can_overlay, "figs/syserr/sqrtS_overlay_Rout2_minus_Rside2");
        delete can_overlay;
    }

    if(graph_txt.is_open())
    {
        graph_txt.close();
    }


    // OUTPUT FILE to store results (keeps previous behavior)
    TFile* outfile = new TFile("syserr_results.root", "RECREATE");
    outfile->cd();
    for(int ie=0; ie<NENERGIES; ie++)
    {
        if(alpha_syserr[ie]) { alpha_syserr[ie]->SetName(Form("alpha_syserr_%s", energies[ie])); alpha_syserr[ie]->Write(); }
        if(R_syserr[ie]) { R_syserr[ie]->SetName(Form("R_syserr_%s", energies[ie])); R_syserr[ie]->Write(); }
        if(N_syserr[ie]) { N_syserr[ie]->SetName(Form("N_syserr_%s", energies[ie])); N_syserr[ie]->Write(); }
        if(R_out_syserr[ie]) { R_out_syserr[ie]->SetName(Form("R_out_syserr_%s", energies[ie])); R_out_syserr[ie]->Write(); }
        if(R_side_syserr[ie]) { R_side_syserr[ie]->SetName(Form("R_side_syserr_%s", energies[ie])); R_side_syserr[ie]->Write(); }
        if(R_long_syserr[ie]) { R_long_syserr[ie]->SetName(Form("R_long_syserr_%s", energies[ie])); R_long_syserr[ie]->Write(); }
    }
    outfile->Close();

    // Write CSV percent breakdown collected during the energy loop
    if(!csv_lines.empty()){
        std::ofstream csvout("syserr_percent_breakdown.csv");
        if(csvout.is_open()){
            for(const auto &ln : csv_lines) csvout << ln << "\n";
            csvout.close();
        } else {
            cerr << "Warning: could not open syserr_percent_breakdown.csv for writing" << endl;
        }
    }

    // Write LaTeX table for selected or all energies
    {
        const int NENERGIES_LATEX = 0; // set to 0 to include all energies
        const int latex_energy_indices[] = {0, 7, 10}; // 3.0, 11.5, 27.0 GeV
        std::vector<std::string> latex_lines;
        std::vector<int> energy_indices;
        if(NENERGIES_LATEX <= 0)
        {
            energy_indices.reserve(NENERGIES);
            for(int ie=0; ie<NENERGIES; ie++) energy_indices.push_back(ie);
        }
        else
        {
            energy_indices.assign(latex_energy_indices, latex_energy_indices + NENERGIES_LATEX);
        }

        auto latex_percent = [](double up, double dn) -> std::string
        {
            std::ostringstream os;
            os << std::fixed << std::setprecision(2) << " +" << up << " / -" << dn;
            return os.str();
        };

        latex_lines.push_back("\\begin{table*}[!htpb]");
        latex_lines.push_back("\\centering");
        latex_lines.push_back("\\tiny");
        latex_lines.push_back("\\caption{Average relative systematic uncertainties of the emission source parameters from different uncertainty sources at all collision energies, in percent. \\label{tab:sysunc_full}}");
        latex_lines.push_back("\\resizebox{\\linewidth}{!}{");
        latex_lines.push_back("\\begin{tabular}{ccccccc}");
        latex_lines.push_back("\\toprule");
        latex_lines.push_back("$\\sqrt{s_{NN}}$ & \\textbf{Source} & $\\alpha$ & $\\lambda^{*}$ & $R_{\\mathrm{out}}$ & $R_\\mathrm{side}$ & $R_\\mathrm{long}$\\\\");
        latex_lines.push_back("\\midrule");

        for(size_t idx=0; idx<energy_indices.size(); idx++)
        {
            int ie = energy_indices[idx];
            if(ie < 0 || ie >= NENERGIES) continue;

            std::ostringstream energylab;
            energylab << std::fixed << std::setprecision(1) << energydouble[ie];
            std::string e = energylab.str();

                        std::string multirow_align = (ie == 0) ? "c" : "m";

                        std::ostringstream row1;
                            row1 << "\\multirow[" << multirow_align << "]{4}{*}{" << e << " GeV} & $Q_{\\rm LCMS}^{\\rm max}$ & "
                  << latex_percent(alpha_qlcms_pct_up_arr[ie], alpha_qlcms_pct_dn_arr[ie]) << " & "
                  << latex_percent(N_qlcms_pct_up_arr[ie], N_qlcms_pct_dn_arr[ie]) << " & "
                  << latex_percent(R_out_qlcms_pct_up_arr[ie], R_out_qlcms_pct_dn_arr[ie]) << " & "
                  << latex_percent(R_side_qlcms_pct_up_arr[ie], R_side_qlcms_pct_dn_arr[ie]) << " & "
                 << latex_percent(R_long_qlcms_pct_up_arr[ie], R_long_qlcms_pct_dn_arr[ie]) << "\\\\";
            latex_lines.push_back(row1.str());

                        std::ostringstream row2;
                            row2 << " & $\\rho_{\\rm fit}^{\\rm max}$ & "
                  << latex_percent(alpha_rhofit_pct_up_arr[ie], alpha_rhofit_pct_dn_arr[ie]) << " & "
                  << latex_percent(N_rhofit_pct_up_arr[ie], N_rhofit_pct_dn_arr[ie]) << " & "
                  << latex_percent(R_out_rhofit_pct_up_arr[ie], R_out_rhofit_pct_dn_arr[ie]) << " & "
                  << latex_percent(R_side_rhofit_pct_up_arr[ie], R_side_rhofit_pct_dn_arr[ie]) << " & "
                 << latex_percent(R_long_rhofit_pct_up_arr[ie], R_long_rhofit_pct_dn_arr[ie]) << "\\\\";
            latex_lines.push_back(row2.str());

            std::ostringstream row3;
              row3 << " & $N_{\\rm events}$ & "//" & \\texttt{NEVT\\_AVG} & "
                  << latex_percent(alpha_nevt_pct_up_arr[ie], alpha_nevt_pct_dn_arr[ie]) << " & "
                  << latex_percent(N_nevt_pct_up_arr[ie], N_nevt_pct_dn_arr[ie]) << " & "
                  << latex_percent(R_out_nevt_pct_up_arr[ie], R_out_nevt_pct_dn_arr[ie]) << " & "
                  << latex_percent(R_side_nevt_pct_up_arr[ie], R_side_nevt_pct_dn_arr[ie]) << " & "
                 << latex_percent(R_long_nevt_pct_up_arr[ie], R_long_nevt_pct_dn_arr[ie]) << "\\\\";
            latex_lines.push_back(row3.str());

            std::ostringstream row4;
              row4 << " & Total & "
                  << latex_percent(alpha_total_pct_up_arr[ie], alpha_total_pct_dn_arr[ie]) << " & "
                  << latex_percent(N_total_pct_up_arr[ie], N_total_pct_dn_arr[ie]) << " & "
                  << latex_percent(R_out_total_pct_up_arr[ie], R_out_total_pct_dn_arr[ie]) << " & "
                  << latex_percent(R_side_total_pct_up_arr[ie], R_side_total_pct_dn_arr[ie]) << " & "
                 << latex_percent(R_long_total_pct_up_arr[ie], R_long_total_pct_dn_arr[ie]) << "\\\\";
            latex_lines.push_back(row4.str());

            if(idx + 1 < energy_indices.size())
            {
                latex_lines.push_back("\\midrule");
            }
        }

        latex_lines.push_back("\\bottomrule");
        latex_lines.push_back("\\end{tabular}}");
        latex_lines.push_back("\\end{table*}");

        std::ofstream texout("syserr_percent_table.tex");
        if(texout.is_open())
        {
            for(const auto& ln : latex_lines) texout << ln << "\n";
            texout.close();
        }
        else
        {
            cerr << "Warning: could not open syserr_percent_table.tex for writing" << endl;
        }
    }

    return 0;
}
