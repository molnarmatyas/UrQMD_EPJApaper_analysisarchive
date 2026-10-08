#include<algorithm>
#include<cmath>
#include<complex>
#include<iomanip>
#include<iostream>
#include<string>
#include<vector>
#include <sstream> // std::stringstream
#include <fstream>
#include <TFile.h>
#include <TTree.h>

int main(int argc, char** argv) {
    // massive speedup reading stdin
    std::cin.sync_with_stdio(false);
    std::cout << std::scientific << std::setprecision(10);
    if(argc < 2)
    {
      std::cerr << "usage: converter_f19 <energy_string> " << std::endl;
      return 0;
    }
    //std::string inputfilename = Form("../auau-%s_test_1000evt_b0-4p73_tmax1000.f19",argv[1]);
    std::string inputfilename = Form("../urqmd-3.4/auau-%s_test_10000evt_b0-4p73_tmax1000.f19",argv[1]);
    std::string line;
    
    // ROOT setup
    TFile* outputFile = new TFile(Form("AuAu_%s_tree.root",argv[1]), "RECREATE");
    TTree* tree = new TTree("urqmd_tree", "UrQMD Tree");
    
    Int_t event = 0;
    Int_t nParts = 0;
    Double_t imp_par = 0.0;
    std::vector<int> id;
    std::vector<double> px;
    std::vector<double> py;
    std::vector<double> pz;
    std::vector<double> E;
    std::vector<double> rx;
    std::vector<double> ry;
    std::vector<double> rz;
    std::vector<double> t;

    // Branch creation
    tree->Branch("eventId", &event);
    tree->Branch("nParts", &nParts);
    tree->Branch("imp", &imp_par);
    tree->Branch("PID", &id);
    tree->Branch("pX", &px);
    tree->Branch("pY", &py);
    tree->Branch("pZ", &pz);
    tree->Branch("E", &E);
    tree->Branch("rX", &rx);
    tree->Branch("rY", &ry);
    tree->Branch("rZ", &rz);
    tree->Branch("t", &t);

    int ipart = 0;
    bool newEvent = false;

    // read stdin
  	std::ifstream file(inputfilename);
    std::cout << inputfilename << std::endl;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        // if header
        if (line.size() == 42) {
            std::cout << "eventId: " << event << std::endl;
            // if event was processed
            if (newEvent && ipart > 0) {
                tree->Fill();
            }
            
            // clear vectors for new event
            id.clear();
            px.clear();
            py.clear();
            pz.clear();
            E.clear();
            rx.clear();
            ry.clear();
            rz.clear();
            t.clear();
            
            event += 1;
            nParts = std::stoi(line.substr(19, 3));
            imp_par = std::stod(line.substr(27, 5));
            ipart = 0;
            newEvent = true;
        }
        else if (line.size() == 148) {
            // particle data
            id.push_back(std::stoi(line.substr(18, 4)));

            px.push_back(std::stod(line.substr(24, 12)));
            py.push_back(std::stod(line.substr(38, 12)));
            pz.push_back(std::stod(line.substr(52, 12)));
            E.push_back(std::stod(line.substr(66, 12)));

            rx.push_back(std::stod(line.substr(94, 12)));
            ry.push_back(std::stod(line.substr(108, 12)));
            rz.push_back(std::stod(line.substr(122, 12)));
            t.push_back(std::stod(line.substr(136, 12)));
            ipart++;
        }
    }

    // Fill the last event if it exists
    if (newEvent && ipart > 0) {
        tree->Fill();
    }

    tree->Write();
    outputFile->Close();
    return 0;
}
