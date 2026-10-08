#ifndef HEADER_FOR_ALL_EMISSIONSOURCE_H
#define HEADER_FOR_ALL_EMISSIONSOURCE_H
#include <cmath> // for taking square root here

const int NFRAME = 3; // lcms, pcms, lab
const int NOSL = 3; // out, side, long, plus rho
const int NCH = 2; // number of charge combinations
const int NCENT = 10; // number of centrality classes

//#include <array>    // std::array
//#include <cstddef>  // std::size_t

const bool is3Dfit = true; // global flag to indicate if 3D fit is being performed

// Mass^2 values
const double Mass2_pi = 0.019479835;
const double Mass2_ka = 0.24371698032;

// kT bins
const int NKT = 10;
const double ktbins[NKT+1] = {0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675};
const double kT_center[NKT] = {0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65};

// Energies
const char* energies[] = {"3p0","3p2","3p5","3p9","4p2","4p5","5p2","6p2","7p2","7p7","9p2","11p5","14p5","19p6","27"};
const int NENERGIES = static_cast<int>(sizeof(energies) / sizeof(energies[0]));
const double energydouble[NENERGIES] = {3.0, 3.2, 3.5, 3.9, 4.2, 4.5, 5.2, 6.2, 7.2, 7.7, 9.2, 11.5, 14.5, 19.6, 27.0};
const int NENERGIES_highstat = 4; // number of the lowest energies with higher statistics (100k instead of 10k evts)
                                  // change it if you want to use higher statistics for more or fewer energies

// from pairsource_urqmd.cc (avoid leading underscore in global scope)
const char* qLCMS_cut[3] = {"default", "strict", "loose"};
const double qLCMS_cut_values[3] = {0.15, 0.05, 0.25}; // GeV/c

// from onedim_EbE_or_Eavg_fit.cc
// TODO rfitmax_systlimits centrality dependence maybe in the future?
//const double B[3] = {2500.0, 1600.0, 3600.0}; // rho_fitmax limits: default, strict, loose - NOT USED
const double rfitmax_def = 80.0;//72.0;
double rfitmax_systlimits[NENERGIES][NKT][3]; // to be calculated in fitting code

// Default number of events to be averaged over for each fit, for each energy; indices correspond to NEVT_AVGsyst array below
const int NEVT_AVGsyst[] = {10, 25, 50, 100, 200, 500, 1000, 5000, 10000, 20000, 25000, 50000, 100000}; // removed 1, not meaningful for low energies, runs too long. Above 10k: highstat
const int NEVTAVGS = sizeof(NEVT_AVGsyst) / sizeof(int);
const int NEVT_AVG_DEFAULT[NENERGIES] = {9,8,7,7,6,6,5,5,5,5,4,4,4,3,3}; // indices in NEVT_AVGsyst corresponding to default

//{60.0, 25.0, 95.0}; // simpler limits for 3D
//const double rfitmax_systlimits[3] = {100.0, 50.0, 150.0}; // for 1D this unified stuff was sufficient;

// bporfy kT bins and centers
// const int NKT = 8;
// const double kT_center[NKT] = {0.0452, 0.0904, 0.1299, 0.1649, 0.2044, 0.2633, 0.3820, 0.6221};
// const double ktbins[NKT + 1] = {0.0, 0.07, 0.11, 0.15, 0.18, 0.23, 0.3, 0.5, 1.0}; // bporfy
// const char* energies[] = {"30A"};
// const int NENERGIES = sizeof(energies) / sizeof(energies[0]);
// const double energydouble[NENERGIES] = {30.0};

const char* centleg[NCENT+2] = {"0-5", "5-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-100","all","0-10"};

#endif // HEADER_FOR_ALL_EMISSIONSOURCE_H
