#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMarker.h>
#include <vector>
#include <string>
#include <cmath>
#include "header_for_all_emissionsource.h"

void circle_defaultnevtavg_point(TGraphErrors* param_vs_nevt_graph, int ienergy)
{
    if(!param_vs_nevt_graph || param_vs_nevt_graph->GetN() <= 0) {
        return;
    }

    double x = NEVT_AVGsyst[NEVT_AVG_DEFAULT[ienergy]];
    int ipoint = -1;
    
    // Find the default point in the graph
    for(int i=0; i<param_vs_nevt_graph->GetN(); i++)
    {
        if(std::abs(param_vs_nevt_graph->GetX()[i] - x) < 1e-9)
        {
            ipoint = i;
            break;
        }
    }

    if(ipoint < 0) {
        return;
    }

    double y = param_vs_nevt_graph->GetY()[ipoint];

    // Open-circle marker size is in screen space, so it stays circular on linear/log axes.
    TMarker* marker = new TMarker(x, y, 107); // 107 is an open circle marker, thicker
    marker->SetMarkerColor(kGreen-8);
    marker->SetMarkerSize(2.8);
    marker->Draw("same");
    
}

// ------------------ MAIN ------------------------
void plot_param_vs_nevt_avg(int ikt =-1) {
    // NEVT_AVG values defined in header_for_all_emissionsource.h as NEVT_AVGsyst array and NEVT_AVG_DEFAULT index
    std::vector<int> nevt_avg;
    for(int i=0; i<NEVTAVGS; i++)
    {
        nevt_avg.push_back(NEVT_AVGsyst[i]);
    }
        
    // Create vectors to store parameters for different energies
    std::vector<TGraphErrors*> alpha_graphs;
    std::vector<TGraphErrors*> R_graphs;
    std::vector<TGraphErrors*> N_graphs;
    std::vector<TGraphErrors*> conflev_graphs;
    std::vector<TGraphErrors*> R_out_graphs;
    std::vector<TGraphErrors*> R_side_graphs;
    std::vector<TGraphErrors*> R_long_graphs;
    std::vector<TGraphErrors*> chi2ndf_graphs;
    
    // Colors for different energies
    int colors[NENERGIES] = {kRed, kBlue, kTeal+3, kMagenta, kOrange+7, kGreen+2, kCyan+2, kViolet+2, kPink+7, kGray+2, kAzure+7};
    
    for(size_t e = 0; e < NENERGIES; e++) {
        std::vector<double> alphas, alphaErrs;
        std::vector<double> Rs, RErrs;
        std::vector<double> Ns, NErrs;
        std::vector<double> conflevs, conflevErrs;
        std::vector<double> nevt_points, nevt_errs;
        std::vector<double> R_outs, R_outErrs;
        std::vector<double> R_sides, R_sideErrs;
        std::vector<double> R_longs, R_longErrs;
        std::vector<double> chi2ndfs, chi2ndfErrs;
        std::vector<double> nevt_points_3d, nevt_errs_3d;
        
        for(int nevt : nevt_avg) {
            // Construct filename
            TString filename = Form("levyfit/results/UrQMD_onedfitresults_lcms_cent0-10_%s_AVG%d.root",
                                  energies[e], nevt);
            
            TFile* file = TFile::Open(filename);
            if(!file) {
                printf("Could not open file: %s\n", filename.Data());
                continue;
            }
            
            // Get the histogram (using alpha_vs_R_all for all kT bins combined)
            const char* histname = "alpha_vs_R_all";
            if(ikt>=0)
            {
                histname = Form("alpha_vs_R_ikt%i", ikt);
            }
            TH2F* h2 = (TH2F*)file->Get(histname);
            if(!h2) {
                printf("Could not find alpha_vs_R_all/alpha_vs_R_ikt# histogram in file: %s\n", filename.Data());
                file->Close();
                continue;
            }
            // Check if histogram has entries
            if(h2->GetEntries() == 0) {
                printf("Histogram alpha_vs_R_all/alpha_vs_R_ikt# has no entries in file: %s\n", filename.Data());
                file->Close();
                continue;
            }
            
            // Get mean values and errors
            double alpha = h2->GetMean(1);  // X axis = alpha
            double alphaErr = h2->GetRMS(1) / sqrt(h2->GetEntries());
            double R = h2->GetMean(2);      // Y axis = R
            double RErr = h2->GetRMS(2) / sqrt(h2->GetEntries());
            
            // Same for lambda (in other convention, displayed here as N) and confidence level
            // N is for diff. ikt only
            int ikt_N = ikt;
            if(ikt_N < 0) ikt_N = 0; // just to get the first histogram, later will be added together
            TH1F* hN = (TH1F*)file->Get(Form("Nhist_ikt%i", ikt_N));
            if(!hN) {
                printf("Could not find Nhist_ikt# histogram in file: %s\n", filename.Data());
                file->Close();
                continue;
            }
            // Check if combining all kT bins (original ikt was < 0)
            if(ikt < 0) {
                // Add all kT bins together
                for(int ikt_add = 1; ikt_add < NKT; ikt_add++) {
                    TH1F* hN_add = (TH1F*)file->Get(Form("Nhist_ikt%i", ikt_add));
                    if(hN_add) {
                        hN->Add(hN_add);
                    }
                }
                // Rename to avoid confusion
                hN->SetName("Nhist_allkt");
                hN->SetTitle("N, all kT bins combined");
            }
            
            // Check if hN has entries before computing mean/RMS
            if(hN->GetEntries() == 0) {
                printf("Nhist_ikt# is empty in file: %s\n", filename.Data());
                file->Close();
                continue;
            }
            
            double N = hN->GetMean();
            double NErr = hN->GetRMS() / sqrt(hN->GetEntries());
            // Confidence level is stored together and for all kT bins separately as well
            const char* conflev_histname = (ikt>=0) ? Form("confidencehist_ikt%i", ikt) : "confidencehist_all";
            TH1F* hConfLev = (TH1F*)file->Get(conflev_histname);
            if(!hConfLev) {
                printf("Could not find confidence histogram in file: %s\n", filename.Data());
                file->Close();
                continue;
            }
            if(hConfLev->GetEntries() == 0) {
                printf("Confidence histogram is empty in file: %s\n", filename.Data());
                file->Close();
                continue;
            }
            double conflev = hConfLev->GetMean();
            double conflevErr = hConfLev->GetRMS() / sqrt(hConfLev->GetEntries());

            // 3D parameters and chi2/NDF for a second 4-panel plot.
            bool has_all_3d_inputs = true;

            int ikt_3d = ikt;
            if(ikt_3d < 0) ikt_3d = 0;

            TH1F* hRout = (TH1F*)file->Get(Form("R_out_hist_ikt%i", ikt_3d));
            TH1F* hRside = (TH1F*)file->Get(Form("R_side_hist_ikt%i", ikt_3d));
            TH1F* hRlong = (TH1F*)file->Get(Form("R_long_hist_ikt%i", ikt_3d));

            if(hRout && hRside && hRlong && ikt < 0) {
                for(int ikt_add = 1; ikt_add < NKT; ikt_add++) {
                    TH1F* hRout_add = (TH1F*)file->Get(Form("R_out_hist_ikt%i", ikt_add));
                    TH1F* hRside_add = (TH1F*)file->Get(Form("R_side_hist_ikt%i", ikt_add));
                    TH1F* hRlong_add = (TH1F*)file->Get(Form("R_long_hist_ikt%i", ikt_add));
                    if(hRout_add && hRside_add && hRlong_add) {
                        hRout->Add(hRout_add);
                        hRside->Add(hRside_add);
                        hRlong->Add(hRlong_add);
                    }
                }
            }

            if(!hRout || !hRside || !hRlong) {
                has_all_3d_inputs = false;
            }

            if(has_all_3d_inputs) {
                if(hRout->GetEntries() == 0 || hRside->GetEntries() == 0 || hRlong->GetEntries() == 0) {
                    has_all_3d_inputs = false;
                }
            }

            const char* chi2ndf_histname = (ikt>=0) ? Form("chi2ndfhist_ikt%i", ikt) : "chi2ndfhist_all";
            TH1F* hChi2Ndf = (TH1F*)file->Get(chi2ndf_histname);
            if(!hChi2Ndf && ikt < 0) {
                hChi2Ndf = (TH1F*)file->Get("chi2ndfhist_ikt0");
                if(hChi2Ndf) {
                    for(int ikt_add = 1; ikt_add < NKT; ikt_add++) {
                        TH1F* hChi2Ndf_add = (TH1F*)file->Get(Form("chi2ndfhist_ikt%i", ikt_add));
                        if(hChi2Ndf_add) {
                            hChi2Ndf->Add(hChi2Ndf_add);
                        }
                    }
                }
            }

            if(!hChi2Ndf || hChi2Ndf->GetEntries() == 0) {
                has_all_3d_inputs = false;
            }

            // Push all quantities only after all required histograms were found.
            alphas.push_back(alpha);
            alphaErrs.push_back(alphaErr);
            Rs.push_back(R);
            RErrs.push_back(RErr);
            Ns.push_back(N);
            NErrs.push_back(NErr);
            conflevs.push_back(conflev);
            conflevErrs.push_back(conflevErr);
            nevt_points.push_back(nevt);
            nevt_errs.push_back(0);

            if(has_all_3d_inputs) {
                double R_out = hRout->GetMean();
                double R_outErr = hRout->GetRMS() / sqrt(hRout->GetEntries());
                double R_side = hRside->GetMean();
                double R_sideErr = hRside->GetRMS() / sqrt(hRside->GetEntries());
                double R_long = hRlong->GetMean();
                double R_longErr = hRlong->GetRMS() / sqrt(hRlong->GetEntries());
                double chi2ndf = hChi2Ndf->GetMean();
                double chi2ndfErr = hChi2Ndf->GetRMS() / sqrt(hChi2Ndf->GetEntries());

                R_outs.push_back(R_out);
                R_outErrs.push_back(R_outErr);
                R_sides.push_back(R_side);
                R_sideErrs.push_back(R_sideErr);
                R_longs.push_back(R_long);
                R_longErrs.push_back(R_longErr);
                // Skip dropped chi2/NDF points that appear at y=0 due to upstream overflow handling.
                // Keep other 3D points unaffected.
                if(chi2ndf > 1e-12) {
                    chi2ndfs.push_back(chi2ndf);
                    chi2ndfErrs.push_back(chi2ndfErr);
                    nevt_points_3d.push_back(nevt);
                    nevt_errs_3d.push_back(0);
                }
            }
            
            file->Close();
        }
        
        // Create graphs
        TGraphErrors* g_alpha = nullptr;
        TGraphErrors* g_R = nullptr;
        TGraphErrors* g_N = nullptr;
        TGraphErrors* g_conflev = nullptr;
        TGraphErrors* g_R_out = nullptr;
        TGraphErrors* g_R_side = nullptr;
        TGraphErrors* g_R_long = nullptr;
        TGraphErrors* g_chi2ndf = nullptr;
        if(!nevt_points.empty()) {
            g_alpha = new TGraphErrors(nevt_points.size(),
                nevt_points.data(), alphas.data(), nevt_errs.data(), alphaErrs.data());
            g_R = new TGraphErrors(nevt_points.size(),
                nevt_points.data(), Rs.data(), nevt_errs.data(), RErrs.data());
            g_N = new TGraphErrors(nevt_points.size(),
                nevt_points.data(), Ns.data(), nevt_errs.data(), NErrs.data());
            g_conflev = new TGraphErrors(nevt_points.size(),
                nevt_points.data(), conflevs.data(), nevt_errs.data(), conflevErrs.data());
        }
        if(!nevt_points_3d.empty()) {
            g_R_out = new TGraphErrors(nevt_points_3d.size(),
                nevt_points_3d.data(), R_outs.data(), nevt_errs_3d.data(), R_outErrs.data());
            g_R_side = new TGraphErrors(nevt_points_3d.size(),
                nevt_points_3d.data(), R_sides.data(), nevt_errs_3d.data(), R_sideErrs.data());
            g_R_long = new TGraphErrors(nevt_points_3d.size(),
                nevt_points_3d.data(), R_longs.data(), nevt_errs_3d.data(), R_longErrs.data());
            g_chi2ndf = new TGraphErrors(nevt_points_3d.size(),
                nevt_points_3d.data(), chi2ndfs.data(), nevt_errs_3d.data(), chi2ndfErrs.data());
        }
            
        // Use solid markers and slightly transparent, thicker lines to connect points
        int baseCol = colors[e];
        int transCol = TColor::GetColorTransparent(baseCol, 0.4); // 40% opacity

        if(g_alpha) {
            g_alpha->SetMarkerStyle(20);
            g_alpha->SetMarkerColor(baseCol);
            g_alpha->SetMarkerSize(1.0);
            g_alpha->SetLineColor(transCol);
            g_alpha->SetLineWidth(3);

            g_R->SetMarkerStyle(20);
            g_R->SetMarkerColor(baseCol);
            g_R->SetMarkerSize(1.0);
            g_R->SetLineColor(transCol);
            g_R->SetLineWidth(3);

            g_N->SetMarkerStyle(20);
            g_N->SetMarkerColor(baseCol);
            g_N->SetMarkerSize(1.0);
            g_N->SetLineColor(transCol);
            g_N->SetLineWidth(3);

            g_conflev->SetMarkerStyle(20);
            g_conflev->SetMarkerColor(baseCol);
            g_conflev->SetMarkerSize(1.0);
            g_conflev->SetLineColor(transCol);
            g_conflev->SetLineWidth(3);
        }
        if(g_R_out) {
            g_R_out->SetMarkerStyle(20);
            g_R_out->SetMarkerColor(baseCol);
            g_R_out->SetMarkerSize(1.0);
            g_R_out->SetLineColor(transCol);
            g_R_out->SetLineWidth(3);

            g_R_side->SetMarkerStyle(20);
            g_R_side->SetMarkerColor(baseCol);
            g_R_side->SetMarkerSize(1.0);
            g_R_side->SetLineColor(transCol);
            g_R_side->SetLineWidth(3);

            g_R_long->SetMarkerStyle(20);
            g_R_long->SetMarkerColor(baseCol);
            g_R_long->SetMarkerSize(1.0);
            g_R_long->SetLineColor(transCol);
            g_R_long->SetLineWidth(3);

            g_chi2ndf->SetMarkerStyle(20);
            g_chi2ndf->SetMarkerColor(baseCol);
            g_chi2ndf->SetMarkerSize(1.0);
            g_chi2ndf->SetLineColor(transCol);
            g_chi2ndf->SetLineWidth(3);
        }
        
        alpha_graphs.push_back(g_alpha);
        R_graphs.push_back(g_R);
        N_graphs.push_back(g_N);
        conflev_graphs.push_back(g_conflev);
        R_out_graphs.push_back(g_R_out);
        R_side_graphs.push_back(g_R_side);
        R_long_graphs.push_back(g_R_long);
        chi2ndf_graphs.push_back(g_chi2ndf);
    }
    
    // Create and divide canvas
    TCanvas* c1 = new TCanvas("c1", "Parameter vs NEVT_AVG", 2400, 2400);
    c1->Divide(2,2);
    
    // Plot alpha
    c1->cd(1);
    gPad->SetLogx();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.14);
    // Change legend position from (0.65, 0.75, 0.85, 0.85) to bottom left
    TLegend* leg1 = new TLegend(0.7, 0.15, 0.9, 0.45);
    bool first_alpha_drawn = false;
    
    for(size_t i = 0; i < alpha_graphs.size(); i++) {
        if(!alpha_graphs[i]) continue;
        if(!first_alpha_drawn) {
            alpha_graphs[i]->SetTitle("Levy #alpha vs NEVT_AVG");
            alpha_graphs[i]->GetXaxis()->SetTitle("NEVT_AVG");
            alpha_graphs[i]->GetYaxis()->SetTitle("#alpha");
            alpha_graphs[i]->GetXaxis()->SetLabelOffset(0.005);
            alpha_graphs[i]->GetYaxis()->SetTitleOffset(1.35);
            alpha_graphs[i]->GetYaxis()->SetRangeUser(1.3, 2.1);
            // assuming that lowest energies - plotted first - are highstat, this range setting should work
            //alpha_graphs[i]->GetXaxis()->SetRangeUser(nevt_avg[0], nevt_avg[NEVTAVGS-1]);
            alpha_graphs[i]->GetXaxis()->SetRangeUser(nevt_avg[0], nevt_avg.back());
            alpha_graphs[i]->Draw("APE");
            // overlay a connecting line so points remain visible with error bars
            alpha_graphs[i]->Draw("L same");
            first_alpha_drawn = true;
        } else {
            alpha_graphs[i]->Draw("PE same");
            alpha_graphs[i]->Draw("L same");
        }
        circle_defaultnevtavg_point(alpha_graphs[i], i); // circle the default NEVT_AVG point for this energy
        leg1->AddEntry(alpha_graphs[i], Form("%s GeV", energies[i]), "LP");
    }
    leg1->Draw();
    
    // Plot R
    c1->cd(2);
    gPad->SetLogx();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.14);
    // Change legend position from (0.65, 0.75, 0.85, 0.85) to bottom left
    //TLegend* leg2 = new TLegend(0.1, 0.1, 0.3, 0.4);
    // Place R legend in the upper-right with same size as alpha legend (width=0.2, height=0.3)
    TLegend* leg2 = new TLegend(0.7, 0.6, 0.9, 0.9);
    bool first_R_drawn = false;
    
    for(size_t i = 0; i < R_graphs.size(); i++) {
        if(!R_graphs[i]) continue;
        if(!first_R_drawn) {
            R_graphs[i]->SetTitle("Levy R vs NEVT_AVG");
            R_graphs[i]->GetXaxis()->SetTitle("NEVT_AVG");
            R_graphs[i]->GetYaxis()->SetTitle("R [fm]");
            R_graphs[i]->GetXaxis()->SetLabelOffset(0.005);
            R_graphs[i]->GetYaxis()->SetTitleOffset(1.35);
            R_graphs[i]->GetYaxis()->SetRangeUser(3, 7); // R between 0 and 10 fm
            R_graphs[i]->GetXaxis()->SetRangeUser(nevt_avg[0], nevt_avg.back());
            R_graphs[i]->Draw("APE");
            R_graphs[i]->Draw("L same");
            first_R_drawn = true;
        } else {
            R_graphs[i]->Draw("PE same");
            R_graphs[i]->Draw("L same");
        }
        circle_defaultnevtavg_point(R_graphs[i], i); // circle the default NEVT_AVG point for this energy
        leg2->AddEntry(R_graphs[i], Form("%s GeV", energies[i]), "LP");
    }
    leg2->Draw();

    // Plot N
    c1->cd(3);
    gPad->SetLogx();
    gPad->SetLeftMargin(0.17);
    gPad->SetBottomMargin(0.14);
    // Place N legend in the upper-right with same size as alpha/R legend (width=0.2, height=0.3)
    TLegend* leg3 = new TLegend(0.7, 0.6, 0.9, 0.9);
    bool first_N_drawn = false;

    for(size_t i = 0; i < N_graphs.size(); i++) {
        if(!N_graphs[i]) continue;
        if(!first_N_drawn) {
            N_graphs[i]->SetTitle("Levy #lambda vs NEVT_AVG");
            N_graphs[i]->GetXaxis()->SetTitle("NEVT_AVG");
            N_graphs[i]->GetYaxis()->SetTitle("#lambda");
            N_graphs[i]->GetXaxis()->SetLabelOffset(0.005);
            N_graphs[i]->GetYaxis()->SetTitleOffset(1.55);
            N_graphs[i]->GetYaxis()->SetRangeUser(0.92, 1.12); // N between 0.98 and 1.12
            N_graphs[i]->GetXaxis()->SetRangeUser(nevt_avg[0], nevt_avg.back());
            N_graphs[i]->Draw("APE");
            N_graphs[i]->Draw("L same");
            first_N_drawn = true;
        } else {
            N_graphs[i]->Draw("PE same");
            N_graphs[i]->Draw("L same");
        }
        circle_defaultnevtavg_point(N_graphs[i], i); // circle the default NEVT_AVG point for this energy
        leg3->AddEntry(N_graphs[i], Form("%s GeV", energies[i]), "LP");
    }
    leg3->Draw();

    // Plot confidence level
    c1->cd(4);
    gPad->SetLogx();
    gPad->SetLeftMargin(0.17);
    gPad->SetBottomMargin(0.14);
    // Place confidence legend in the lower-left with same size as alpha/R legend
    TLegend* leg4 = new TLegend(0.2, 0.18, 0.4, 0.48);
    bool first_conflev_drawn = false;

    for(size_t i = 0; i < conflev_graphs.size(); i++) {
        if(!conflev_graphs[i]) continue;
        if(!first_conflev_drawn) {
            conflev_graphs[i]->SetTitle("Fit Confidence Level vs NEVT_AVG");
            conflev_graphs[i]->GetXaxis()->SetTitle("NEVT_AVG");
            conflev_graphs[i]->GetYaxis()->SetTitle("Confidence level");
            conflev_graphs[i]->GetXaxis()->SetLabelOffset(0.005);
            conflev_graphs[i]->GetYaxis()->SetTitleOffset(1.55);
            conflev_graphs[i]->GetYaxis()->SetRangeUser(1e-15, 1.05); // Confidence level between 0 and 1
            conflev_graphs[i]->GetXaxis()->SetRangeUser(nevt_avg[0], nevt_avg.back());
            // Set log scale for y-axis
            gPad->SetLogy();
            conflev_graphs[i]->Draw("APE");
            conflev_graphs[i]->Draw("L same");
            first_conflev_drawn = true;
        } else {
            conflev_graphs[i]->Draw("PE same");
            conflev_graphs[i]->Draw("L same");
        }
        circle_defaultnevtavg_point(conflev_graphs[i], i); // circle the default NEVT_AVG point for this energy
        leg4->AddEntry(conflev_graphs[i], Form("%s GeV", energies[i]), "LP");
    }
    leg4->Draw();

    // Final, save
    const char* ikt_suffix = (ikt>=0) ? Form("_ikt%i", ikt) : "_allkt";
    c1->SaveAs(Form("figs/levy_params_vs_nevt_avg%s.png", ikt_suffix));

    // Auto-tuned y-ranges for the second canvas.
    double rmin = 1e99, rmax = -1e99;
    for(size_t i = 0; i < R_out_graphs.size(); i++) {
        TGraphErrors* gset[3] = {R_out_graphs[i], R_side_graphs[i], R_long_graphs[i]};
        for(int ig = 0; ig < 3; ig++) {
            if(!gset[ig]) continue;
            for(int ip = 0; ip < gset[ig]->GetN(); ip++) {
                double y = gset[ig]->GetY()[ip];
                double ey = gset[ig]->GetEY()[ip];
                if(y - ey < rmin) rmin = y - ey;
                if(y + ey > rmax) rmax = y + ey;
            }
        }
    }
    if(rmax <= rmin) {
        rmin = 2.5;
        rmax = 10.0;
    } else {
        double rpad = 0.08 * (rmax - rmin);
        rmin -= rpad;
        rmax += rpad;
        if(rmin < 0.0) rmin = 0.0;
    }

    double cmin = 1e99, cmax = -1e99;
    for(size_t i = 0; i < chi2ndf_graphs.size(); i++) {
        if(!chi2ndf_graphs[i]) continue;
        for(int ip = 0; ip < chi2ndf_graphs[i]->GetN(); ip++) {
            double y = chi2ndf_graphs[i]->GetY()[ip];
            double ey = chi2ndf_graphs[i]->GetEY()[ip];
            if(y - ey < cmin) cmin = y - ey;
            if(y + ey > cmax) cmax = y + ey;
        }
    }
    if(cmax <= cmin) {
        cmin = 0.0;
        cmax = 5.0;
    } else {
        double cpad = 0.10 * (cmax - cmin);
        cmin -= cpad;
        cmax += cpad;
        if(cmin < 0.0) cmin = 0.0;
    }

    // Second 4-panel plot: R_out, R_side, R_long, chi2/NDF.
    TCanvas* c2 = new TCanvas("c2", "3D Levy Parameters vs NEVT_AVG", 2400, 2400);
    c2->Divide(2,2);

    // Plot R_out
    c2->cd(1);
    gPad->SetLogx();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.14);
    TLegend* leg5 = new TLegend(0.7, 0.6, 0.9, 0.9);
    bool first_Rout_drawn = false;
    for(size_t i = 0; i < R_out_graphs.size(); i++) {
        if(!R_out_graphs[i]) continue;
        if(!first_Rout_drawn) {
            R_out_graphs[i]->SetTitle("Levy R_{out} vs NEVT_AVG");
            R_out_graphs[i]->GetXaxis()->SetTitle("NEVT_AVG");
            R_out_graphs[i]->GetYaxis()->SetTitle("R_{out} [fm]");
            R_out_graphs[i]->GetXaxis()->SetLabelOffset(0.005);
            R_out_graphs[i]->GetYaxis()->SetTitleOffset(1.35);
            R_out_graphs[i]->GetYaxis()->SetRangeUser(rmin, rmax);
            R_out_graphs[i]->GetXaxis()->SetRangeUser(nevt_avg[0], nevt_avg.back());
            R_out_graphs[i]->Draw("APE");
            R_out_graphs[i]->Draw("L same");
            first_Rout_drawn = true;
        } else {
            R_out_graphs[i]->Draw("PE same");
            R_out_graphs[i]->Draw("L same");
        }
        circle_defaultnevtavg_point(R_out_graphs[i], i);
        leg5->AddEntry(R_out_graphs[i], Form("%s GeV", energies[i]), "LP");
    }
    leg5->Draw();

    // Plot R_side
    c2->cd(2);
    gPad->SetLogx();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.14);
    TLegend* leg6 = new TLegend(0.7, 0.6, 0.9, 0.9);
    bool first_Rside_drawn = false;
    for(size_t i = 0; i < R_side_graphs.size(); i++) {
        if(!R_side_graphs[i]) continue;
        if(!first_Rside_drawn) {
            R_side_graphs[i]->SetTitle("Levy R_{side} vs NEVT_AVG");
            R_side_graphs[i]->GetXaxis()->SetTitle("NEVT_AVG");
            R_side_graphs[i]->GetYaxis()->SetTitle("R_{side} [fm]");
            R_side_graphs[i]->GetXaxis()->SetLabelOffset(0.005);
            R_side_graphs[i]->GetYaxis()->SetTitleOffset(1.35);
            R_side_graphs[i]->GetYaxis()->SetRangeUser(rmin, rmax);
            R_side_graphs[i]->GetXaxis()->SetRangeUser(nevt_avg[0], nevt_avg.back());
            R_side_graphs[i]->Draw("APE");
            R_side_graphs[i]->Draw("L same");
            first_Rside_drawn = true;
        } else {
            R_side_graphs[i]->Draw("PE same");
            R_side_graphs[i]->Draw("L same");
        }
        circle_defaultnevtavg_point(R_side_graphs[i], i);
        leg6->AddEntry(R_side_graphs[i], Form("%s GeV", energies[i]), "LP");
    }
    leg6->Draw();

    // Plot R_long
    c2->cd(3);
    gPad->SetLogx();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.14);
    TLegend* leg7 = new TLegend(0.7, 0.6, 0.9, 0.9);
    bool first_Rlong_drawn = false;
    for(size_t i = 0; i < R_long_graphs.size(); i++) {
        if(!R_long_graphs[i]) continue;
        if(!first_Rlong_drawn) {
            R_long_graphs[i]->SetTitle("Levy R_{long} vs NEVT_AVG");
            R_long_graphs[i]->GetXaxis()->SetTitle("NEVT_AVG");
            R_long_graphs[i]->GetYaxis()->SetTitle("R_{long} [fm]");
            R_long_graphs[i]->GetXaxis()->SetLabelOffset(0.005);
            R_long_graphs[i]->GetYaxis()->SetTitleOffset(1.35);
            R_long_graphs[i]->GetYaxis()->SetRangeUser(rmin, rmax);
            R_long_graphs[i]->GetXaxis()->SetRangeUser(nevt_avg[0], nevt_avg.back());
            R_long_graphs[i]->Draw("APE");
            R_long_graphs[i]->Draw("L same");
            first_Rlong_drawn = true;
        } else {
            R_long_graphs[i]->Draw("PE same");
            R_long_graphs[i]->Draw("L same");
        }
        circle_defaultnevtavg_point(R_long_graphs[i], i);
        leg7->AddEntry(R_long_graphs[i], Form("%s GeV", energies[i]), "LP");
    }
    leg7->Draw();

    // Plot chi2/NDF
    c2->cd(4);
    gPad->SetLogx();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.14);
    TLegend* leg8 = new TLegend(0.18, 0.6, 0.38, 0.9);
    bool first_chi2ndf_drawn = false;
    for(size_t i = 0; i < chi2ndf_graphs.size(); i++) {
        if(!chi2ndf_graphs[i]) continue;
        if(!first_chi2ndf_drawn) {
            chi2ndf_graphs[i]->SetTitle("Fit #chi^{2}/NDF vs NEVT_AVG");
            chi2ndf_graphs[i]->GetXaxis()->SetTitle("NEVT_AVG");
            chi2ndf_graphs[i]->GetYaxis()->SetTitle("#chi^{2}/NDF");
            chi2ndf_graphs[i]->GetXaxis()->SetLabelOffset(0.005);
            chi2ndf_graphs[i]->GetYaxis()->SetTitleOffset(1.35);
            chi2ndf_graphs[i]->GetYaxis()->SetRangeUser(cmin, cmax);
            chi2ndf_graphs[i]->GetXaxis()->SetRangeUser(nevt_avg[0], nevt_avg.back());
            chi2ndf_graphs[i]->Draw("APE");
            chi2ndf_graphs[i]->Draw("L same");
            first_chi2ndf_drawn = true;
        } else {
            chi2ndf_graphs[i]->Draw("PE same");
            chi2ndf_graphs[i]->Draw("L same");
        }
        circle_defaultnevtavg_point(chi2ndf_graphs[i], i);
        leg8->AddEntry(chi2ndf_graphs[i], Form("%s GeV", energies[i]), "LP");
    }
    leg8->Draw();

    c2->SaveAs(Form("figs/levy_3dparams_vs_nevt_avg%s.png", ikt_suffix));
}
