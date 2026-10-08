// After the urqmd .f19 output files have been converted with convert_f19.cc (thus, creating the ...._tree.root files):
// This code does not fit anything, it just makes histograms of the spatial separation of pion pairs in the LCMS frame, for different kT bins.
// TODO: add bim limits for centrality selection when minBias data will be available

#include <TH1.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>

#include <stdlib.h>

#include <iostream>
#include <cstring>
#include <vector>
#include <map>
#include <fstream>

#include "../header_for_all_emissionsource.h"
// const int NKT = 10;
// const double ktbins[NKT+1]= {0.175,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675};
int barWidth = 70;

#ifndef NCH
#define NCH 2
#endif
//const double Mass2_pi = 0.019479835;

//const double ktbins[NKT+1]= {0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50};
//const double qmax[NKT] =       {0.05, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26}; // deprecated I guess
//const int NCENT = 10; // number of centrality classes
//const char* centleg[NCENT+2] = {"0-5", "5-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-100","all","0-10"};

int kTBinNum(double KT)
{
  if (KT > ktbins[NKT])
    return NKT - 1;
  if (KT < ktbins[0])
    return 0;

  int ibin = 0;
  while (KT > ktbins[ibin + 1])
    ibin++;
  return ibin;
}

void progressbar(Int_t event, Int_t nEvents)
{
  float progress = event*1.0 / nEvents * 1.0;
  std::cout << "[";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
      if (i < pos) std::cout << "=";
      else if (i == pos) std::cout << ">";
      else std::cout << " ";
  }
  std::cout << "] " << int(progress * 100.0) << "% ( " << event << " / " << nEvents << " )" << "\r";
  std::cout.flush();
  
}

// Reimplemented in header_for_all_emissionsource.h:
// const char* qLCMS_cut[3] = {"default", "strict", "loose"};
// const double qLCMS_cut_values[3] = {0.15, 0.05, 0.25}; // in GeV/c
const char* qLCMS_outputtext[3] = {"", "_strictqLCMS", "_looseqLCMS"};


int main(int argc, char** argv)
{
  std::map<int, std::map<int, std::map<int, TH1F *> > > D_lcms;
  std::map<int, std::map<int, std::map<int, TH1F *> > > D_out_lcms;
  std::map<int, std::map<int, std::map<int, TH1F *> > > D_side_lcms;
  std::map<int, std::map<int, std::map<int, TH1F *> > > D_long_lcms;
  std::map<int, std::map<int, std::map<int, TH1F *> > > D_outlong_lcms;
  if(argc < 3)
  {
    std::cerr << "usage: pairsource_urqmd <energy_string> <qLCMS_cut>" << std::endl;
    return 0;
  }

  int qLCMS_cuttype = (int)atoi(argv[2]);
  if(qLCMS_cuttype < 0 || qLCMS_cuttype > 2)
  {
    std::cerr << "qLCMS_cuttype must be 0 (default), 1 (strict) or 2 (loose)" << std::endl;
    qLCMS_cuttype = 0;
    return 0;
  }
  std::cout << "Using qLCMS cut type: " << qLCMS_cut[qLCMS_cuttype] << std::endl;

  TFile *treefile = new TFile(Form("AuAu_%s_tree.root",argv[1]),"read");
  //TString inputFile = treefile->GetName(); 
  TFile* outFile = new TFile(Form("UrQMD_3d_source_%scent_all_%s%s.root",centleg[11],argv[1],
                                  qLCMS_outputtext[qLCMS_cuttype]),"RECREATE"); // FIXME make 11->ICENT for other centralities
  // File compression - no decrease in size, but was worth trying
  // outFile->SetCompressionLevel(9); // maximum compression - smaller output, but this might consume more memory and might become slower 

  TTree *event = (TTree*) treefile->Get("urqmd_tree");
  const int Nbins = 200;
  double xmin = 1.e-1;
  double xmax = 1.e3;
  double binfactor = TMath::Power(xmax/xmin,1.0/Nbins);
  double logbins[Nbins+1];
  logbins[0] = xmin;
  for(int ibin = 1; ibin <= Nbins; ibin++)
      logbins[ibin] = xmin * TMath::Power(binfactor,ibin);

  TH1D *KTdist = new TH1D("ktdist","K_{T} distribution",100,0,10);
  TH1D *rapdist = new TH1D("rapdist","rapidity distribution",100,-10,10);
  TH1D *post_rapdist = new TH1D("post_rapdist","post-cut rapidity distribution",100,-10,10);
  TH1D *etadist = new TH1D("etadist","pseudorapidity distribution",100,-10,10);
  TH1D *post_etadist = new TH1D("post_etadist","post-cut pseudorapidity distribution",100,-10,10);

  std::cout << "File " << " loaded" << std::endl;
  std::vector<double> *PX = 0;
  std::vector<double> *PY = 0;
  std::vector<double> *PZ = 0;
  std::vector<double> *energy = 0;
  std::vector<double> *rx = 0; 
  std::vector<double> *ry = 0; 
  std::vector<double> *rz = 0; 
  std::vector<double> *time = 0;
  std::vector<int> *ID = 0;
  event->SetBranchAddress("pX",&PX);
  event->SetBranchAddress("pY",&PY);
  event->SetBranchAddress("pZ",&PZ);
  event->SetBranchAddress("E",&energy);
  event->SetBranchAddress("rX",&rx);
  event->SetBranchAddress("rY",&ry);
  event->SetBranchAddress("rZ",&rz);
  event->SetBranchAddress("t",&time);
  event->SetBranchAddress("PID",&ID);

  Int_t nEvent = event->GetEntries();
  for (Int_t ievent = 0; ievent < nEvent; ievent++) {
    event->GetEntry(ievent);
    progressbar(ievent+1.0, nEvent);
    
    // Create histograms for this event only (incremental processing)
    for (int ich = 0; ich < NCH; ich++)
    {
      for (int iKT = 0; iKT < NKT; iKT++)
      {
        D_lcms[ievent][ich][iKT] = new TH1F(Form("D_lcms_ev%i_ch%i_KT%i", ievent, ich, iKT),"Spatial size of charged pion source in the LCMS",Nbins,logbins);
        D_out_lcms[ievent][ich][iKT] = new TH1F(Form("D_out_lcms_ev%i_ch%i_KT%i", ievent, ich, iKT),"Spatial size, out direction of charged pion source in the LCMS",Nbins,logbins);
        D_side_lcms[ievent][ich][iKT] = new TH1F(Form("D_side_lcms_ev%i_ch%i_KT%i", ievent, ich, iKT),"Spatial size, side direction of charged pion source in the LCMS",Nbins,logbins);
        D_long_lcms[ievent][ich][iKT] = new TH1F(Form("D_long_lcms_ev%i_ch%i_KT%i", ievent, ich, iKT),"Spatial size, long direction of charged pion source in the LCMS",Nbins,logbins);
        D_outlong_lcms[ievent][ich][iKT] = new TH1F(Form("D_outlong_lcms_ev%i_ch%i_KT%i", ievent, ich, iKT),"Spatial size, outlong direction of charged pion source in the LCMS",Nbins,logbins);
      }
    }
    
    int eventsize = PX->size();
    int pid = 0;
    for(int ich = 0; ich < NCH; ich++)
    {
      (ich == 0) ? pid = -211 : pid = 211;
      //std::cout << "Now pid: " << pid << std::endl;
      for (int j = 0; j < eventsize; j++)
      {
        if (ID->at(j) == pid /* || ID->at(j) == -211 */) {
          auto Ej = energy->at(j); //TMath::Sqrt(pj*pj + Mass2_pi);
          for (int k = j+1; k < eventsize; k++)
          {
            if (ID->at(k) == pid /* || ID->at(k) == -211 */) {
              //auto pTk = TMath::Sqrt(PX->at(k)*PX->at(k) + PY->at(k)*PY->at(k));
              //auto pk = TMath::Sqrt(pTk*pTk + PZ->at(k)*PZ->at(k));
              auto Ek = energy->at(k); //TMath::Sqrt(pk*pk + Mass2_pi);

              auto t = time->at(j) - time->at(k);
              auto Rx = rx->at(j) - rx->at(k);
              auto Ry = ry->at(j) - ry->at(k);
              auto Rz = rz->at(j) - rz->at(k);
                

              auto K0 = 0.5 * (Ej + Ek);
              auto Kx = 0.5 * (PX->at(j) + PX->at(k));
              auto Ky = 0.5 * (PY->at(j) + PY->at(k));
              auto Kz = 0.5 * (PZ->at(j) + PZ->at(k));
                
              auto Qx = PX->at(j) - PX->at(k);
              auto Qy = PY->at(j) - PY->at(k);

              auto qZLCMS2 = (PZ->at(j)*Ek-PZ->at(k)*Ej)*(PZ->at(j)*Ek-PZ->at(k)*Ej) / (K0*K0 - Kz*Kz);
              auto qLCMS = TMath::Sqrt(Qx*Qx + Qy*Qy + qZLCMS2);

              auto KT = TMath::Sqrt(Kx*Kx + Ky*Ky);
              auto cosphi = Kx / KT;
              auto sinphi = Ky / KT;

              auto rhoout = Rx * cosphi + Ry * sinphi - (KT / (K0*K0 - Kz*Kz)) * (K0 * t - Kz * Rz);
              auto rhoside = - Rx * sinphi + Ry * cosphi;
              auto rholong = (K0 * Rz - Kz * t) / TMath::Sqrt(K0*K0 - Kz*Kz);

              auto rhooutlong = std::abs(1.0/std::sqrt(2.0)*(rhoout + rholong));
              auto rho2 = rhoout*rhoout + rhoside*rhoside + rholong*rholong;

              KTdist->Fill(KT);
              int iKT = kTBinNum(KT);
              /*
              if((iKT == 0 || iKT == 15) && qLCMS < ktbins[iKT+1])
                std::cout << "qlcms: " << qLCMS << " kT: " << ktbins[iKT] << "-" << ktbins[iKT+1] << std::endl;
              */
              
              // CUTS
              double pTj = TMath::Sqrt(PX->at(j)*PX->at(j) + PY->at(j)*PY->at(j));
              double pTk = TMath::Sqrt(PX->at(k)*PX->at(k) + PY->at(k)*PY->at(k));
              if(pTj<0.15 || pTj>=1.0) continue;
              if(pTk<0.15 || pTk>=1.0) continue;
              double pzj = PZ->at(j);
              double pzk = PZ->at(k);
              double pj = TMath::Sqrt(pTj*pTj + pzj*pzj);
              double pk = TMath::Sqrt(pTk*pTk + pzk*pzk);
              double etaj = 0.5 * TMath::Log(TMath::Abs((pj + pzj) / (pj - pzj)));
              double etak = 0.5 * TMath::Log(TMath::Abs((pk + pzk) / (pk - pzk)));
              double yj = 0.5 * TMath::Log((Ej + pzj) / (Ej - pzj)) - 2.02172; // so far only manually, based on bporfy TODO boost automatically to CMS
              double yk = 0.5 * TMath::Log((Ek + pzk) / (Ek - pzk)) - 2.02172;
              if(k == j+1) // fill only once per pair
              {
                etadist->Fill(etaj);
                rapdist->Fill(yj);
              }
              if(fabs(etaj) > 1.) continue;
              if(fabs(etak) > 1.) continue;
              //if(fabs(yj) > 1.) continue;
              //if(fabs(yk) > 1.) continue;
              if(k == j+1) // fill only once per pair
              {
                post_etadist->Fill(etaj);
                post_rapdist->Fill(yj);
              }
              if(qLCMS < TMath::Sqrt(qLCMS_cut_values[qLCMS_cuttype] * TMath::Sqrt(KT*KT + Mass2_pi))) // 150 Mev * mT baseline; strict 50 MeV, loose 250 MeV
              //if(qLCMS < qmax[iKT]) 
              //if(qLCMS < ktbins[iKT+1]) 
              {
                D_lcms[ievent][ich][iKT]->Fill(TMath::Sqrt(rho2));
                // absolute value is needed, otherwise all negative values go into underflow bin !!!
                D_out_lcms[ievent][ich][iKT]->Fill(TMath::Abs(rhoout));
                D_side_lcms[ievent][ich][iKT]->Fill(TMath::Abs(rhoside));
                D_long_lcms[ievent][ich][iKT]->Fill(TMath::Abs(rholong));
                D_outlong_lcms[ievent][ich][iKT]->Fill(TMath::Abs(rhooutlong));
              }
            } // only particles with matching PID
          } // Loop over pairs - loop for particle k
        } // only particles with matching PID
      } // end of track loop - loop for particle j
    } // end of charge loop
    
    // Write and delete histograms for this event immediately (incremental processing)
    // Only write non-empty histograms to save disk space
    outFile->cd();
    for(int ich = 0; ich < NCH; ich++)
    {
      for(int iKT = 0; iKT < NKT; iKT++)
      {
        //if(D_lcms[ievent][ich][iKT]->GetEntries() == 0) std::cerr << "Event " << ievent << ", charge " << ich << ", kT bin " << iKT << ": empty histogram." << std::endl;
        D_lcms[ievent][ich][iKT]->Write();
        //if(D_out_lcms[ievent][ich][iKT]->GetEntries() == 0) std::cerr << "Event " << ievent << ", out, charge " << ich << ", kT bin " << iKT << ": empty histogram." << std::endl;
        D_out_lcms[ievent][ich][iKT]->Write();
        //if(D_side_lcms[ievent][ich][iKT]->GetEntries() == 0) std::cerr << "Event " << ievent << ", side, charge " << ich << ", kT bin " << iKT << ": empty histogram." << std::endl;
        D_side_lcms[ievent][ich][iKT]->Write();
        //if(D_long_lcms[ievent][ich][iKT]->GetEntries() == 0) std::cerr << "Event " << ievent << ", long, charge " << ich << ", kT bin " << iKT << ": empty histogram." << std::endl;
        D_long_lcms[ievent][ich][iKT]->Write();
        //if(D_outlong_lcms[ievent][ich][iKT]->GetEntries() == 0) std::cerr << "Event " << ievent << ", charge " << ich << ", kT bin " << iKT << ": empty histogram." << std::endl;
        //D_outlong_lcms[ievent][ich][iKT]->Write(); // not really used

        // Free memory immediately
        delete D_lcms[ievent][ich][iKT];
        delete D_out_lcms[ievent][ich][iKT];
        delete D_side_lcms[ievent][ich][iKT];
        delete D_long_lcms[ievent][ich][iKT];
        delete D_outlong_lcms[ievent][ich][iKT];
      }
    }
    // Clear map entries for this event to release map memory
    D_lcms.erase(ievent);
    D_out_lcms.erase(ievent);
    D_side_lcms.erase(ievent);
    D_long_lcms.erase(ievent);
    D_outlong_lcms.erase(ievent);
    
  } // end of event loop
  std::cout << std::endl;
  outFile->cd();
  KTdist->Write();
  rapdist->Write();
  post_rapdist->Write();
  etadist->Write();
  post_etadist->Write();

  outFile->Close();

  return 0;
}
