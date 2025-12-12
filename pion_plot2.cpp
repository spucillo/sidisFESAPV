//--- Authors: Lorenzo Polizzi (lorenzo.polizzi@unife.it), Sara Pucillo (sara.pucillo@cern.ch), Nicolò Valle (nicolo.valle@cern.ch)
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "TText.h"
#include <TLatex.h>
#include <TPaveStats.h>
#include <TPaveStatsEditor.h>
#include "TPaletteAxis.h"
#include "TPolyLine.h"
#include "TStyle.h"
#include "TColor.h"
#include "ePIC_style.C"

using namespace std;
namespace fs = std::filesystem;
gROOT->SetBatch(kTRUE);
// to download the data
// rsync -avz -e "ssh -J lpolizzi@login.jlab.org" lpolizzi@ifarm:/lustre24/expphy/volatile/clas12/lpolizzi/sidis/eic/NC_25.10/ /Users/lorenzopolizzi/Desktop/PhD/epic/25.10_10x100_kaon
// 
// 25.10 bad file:
// hiDiv_2.0003, hiDiv_2.0033, hiDiv_3.0006


// 
vector<double> CreateLogBinning(int nbins, double xmin, double xmax) {
    vector<double> bin_edges(nbins + 1);
    double logxmin = log10(xmin);
    double logxmax = log10(xmax);
    double bin_width = (logxmax - logxmin) / nbins;
    for (int i = 0; i <= nbins; ++i) {
        bin_edges[i] = pow(10, logxmin + i * bin_width);
    }
    return bin_edges;
}

void SetStatsBox(TH2* hist) {
    TPaveStats* stats0 = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
    if (stats0) {
        stats0->SetX1NDC(0.75);  // Posizione pannello (sinistra)
        stats0->SetX2NDC(0.89);  // Posizione pannello (destra)
        stats0->SetY1NDC(0.68);  // Posizione pannello (basso)
        stats0->SetY2NDC(0.88);  // Posizione pannello (alto)
    }
}
void SetStatsBox2(TH2* hist) {
    TPaveStats* stats0 = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
    if (stats0) {
        stats0->SetX1NDC(0.12);  
        stats0->SetX2NDC(0.26);  
        stats0->SetY1NDC(0.68);  
        stats0->SetY2NDC(0.88);  
    }
}

int getBinIndex_xQ2(double xB, double Q2){
    // bin_xB 
    double bin_xB[][2] = {{1e-4, 1e-3}, {1e-3, 3e-3}, {3e-3, 1e-2}, {1e-2, 4e-2}, {4e-2, 1}};
    // bin Q2
    vector<vector<array<double,2>>> binning_Q2_for_xB = {
        // in sequenza i bin di Q2 per i rispettivi bin di xB
        {{1,2}, {2,100}},
        {{1,2}, {2,5}, {5,100}},
        {{1,2}, {2,5}, {5,12}, {12,100}},
        {{1,2}, {2,5}, {5,12}, {12,100}},
        {{1,5}, {5,12}, {12,100}}
    };
    int binIndex = 1;
    for (int ix = 0; ix < 5; ix++){
        if (xB >= bin_xB[ix][0] && xB <= bin_xB[ix][1]){
            for (size_t iq = 0; iq < binning_Q2_for_xB[ix].size(); iq++){
                double Q2min = binning_Q2_for_xB[ix][iq][0];
                double Q2max = binning_Q2_for_xB[ix][iq][1];
                if (Q2 >= Q2min && Q2 < Q2max){
                    return binIndex;
                }
                binIndex++;
            }
            // nel caso non si trovi un valore di Q2
            return -1;
        } else {
            binIndex += binning_Q2_for_xB[ix].size();
        }
    }

    return -1;
}

int getBinIndex_zPt(double z, double Pt){
    // bin_z
    double bin_z[][2] = {{0, 0.05}, {0.05, 0.1}, {0.1, 0.2}, {0.2, 0.3}, {0.3, 0.4}, {0.4, 1}};
    // bin_Pt
    vector<vector<array<double,2>>> binning_Pt_for_z = {
        {{0, 0.2}, {0.2, 0.4}, {0.4, 0.6}, {0.6, 0.8}, {0.8, 2}},
        {{0, 0.2}, {0.2, 0.4}, {0.4, 0.6}, {0.6, 0.8}, {0.8, 2}},
        {{0, 0.2}, {0.2, 0.4}, {0.4, 0.6}, {0.6, 0.8}, {0.8, 2}},
        {{0, 0.2}, {0.2, 0.4}, {0.4, 0.6}, {0.6, 0.8}, {0.8, 2}},
        {{0, 0.2}, {0.2, 0.4}, {0.4, 0.6}, {0.6, 0.8}, {0.8, 2}},
        {{0, 0.2}, {0.2, 0.4}, {0.4, 0.6}, {0.6, 0.8}, {0.8, 2}}
    };
    int binIndex = 1;
    for (int iz = 0; iz < 6; iz++){
        if (z >= bin_z[iz][0] && z <= bin_z[iz][1]){
            for (size_t ip = 0; ip < binning_Pt_for_z[iz].size(); ip++){
                double Ptmin = binning_Pt_for_z[iz][ip][0];
                double Ptmax = binning_Pt_for_z[iz][ip][1];
                if (Pt >= Ptmin && Pt < Ptmax){
                    return binIndex;
                }
                binIndex++;
            }

            return -1;
        } else {
            binIndex += binning_Pt_for_z[iz].size();
        }
    }

    return -1;
}

int getBinIndex_z(double z){
    // bin_z
    double bin_z[][2] = {{0, 0.05}, {0.05, 0.1}, {0.1, 0.2}, {0.2, 0.3}, {0.3, 0.4}, {0.4, 1}};
    int binIndex = 1;
    for (int iz = 0; iz < 6; iz++){
        if (z >= bin_z[iz][0] && z <= bin_z[iz][1]){
            return iz+1;
        }
    }
    return -1;
}
int getBinIndex_Pt(double Pt){
    // bin_Pt
    double bin_Pt[][2] = {{0, 0.2}, {0.2, 0.4}, {0.4, 0.6}, {0.6, 0.8}, {0.8, 2}};
    int binIndex = 1;
    for (int ip = 0; ip < 5; ip++){
        if (Pt >= bin_Pt[ip][0] && Pt <= bin_Pt[ip][1]){
            return ip+1;
        }
    }
    return -1;
}



void pion_plot2() {
    // applico lo stile di ePIC
    //set_ePIC_style();
    gROOT->ProcessLine("set_ePIC_style()");
    // (opzionale) se vuoi assicurarti che sia forzato su tutti i canvas
    //gROOT->ForceStyle();
    // Variables
    // elettrone
    double electron_px, electron_py, electron_pz, electron_mom, electron_Theta, electron_Phi, electron_ThetaDeg, electron_E, electron_W, electron_Q2, electron_ass;
    double electron_eta, electron_y;
    // pione +
    double pionp_mom, pionp_Q2, pionp_xB, pionp_xF, pionp_z, pionp_PhT, pionp_Phi_h, pionp_Phi_s, pionp_Phi_lab, pionp_Theta, pionp_eta, pionp_y, pionp_W, pionp_Mx;
    double helicity, eps, pionp_px, pionp_py, pionp_pz, el_px, el_py, el_pz, el_theta, el_phi, el_eta, el_mom, pr_mom, pr_px, pr_py, pr_pz, pr_phi, pr_theta, pr_eta;
    double rec_pdg, good_PID, el_rec_pdg;
    double el_ass_rec_pdg, el_ass_px, el_ass_py, el_ass_pz, el_ass_theta, el_ass_phi, el_ass_eta, el_ass_mom;
    double pionp_E;
    //
    double pionp_mom_mc, pionp_Q2_mc, pionp_xB_mc, pionp_xF_mc, pionp_z_mc, pionp_PhT_mc;
    double pionp_Phi_h_mc, pionp_Phi_s_mc, pionp_Phi_lab_mc, pionp_Theta_mc, pionp_eta_mc, pionp_y_mc, pionp_W_mc, pionp_Mx_mc;
    double hel_mc, eps_mc, pionp_px_mc, pionp_py_mc, pionp_pz_mc;
    double pionp_mc_index, index_mc;
    // 
    double pionp_mom_all, pionp_Q2_all, pionp_xB_all, pionp_xF_all, pionp_z_all, pionp_PhT_all;
    double pionp_Phi_h_all, pionp_Phi_s_all, pionp_Phi_lab_all, pionp_Theta_all, pionp_eta_all, pionp_y_all, pionp_W_all, pionp_Mx_all;
    double hel_all, eps_all, pionp_px_all, pionp_py_all, pionp_pz_all;
    double good_PID_all, pdg_all;
    double pionp_all_index, index_all;

    // input file
    string inputDir = "25.10_10x100_kaon";
    
    TTree treePionP("Pion+", "");
    TTree treePionP_MC("MC Pion+", "");
    TChain chainElectron("Electron");
    TChain chainHadron_Reco("Hadron Reco");
    TChain chainHadron_MC("Hadron MC");
    int fileCount = 0;
    for (const auto &entry : fs::directory_iterator(inputDir)) {
        if (entry.path().extension() == ".root") {
            string filePath = entry.path().string();
            chainElectron.Add(Form("%s/ElectronTree_MC", filePath.c_str()));
            chainHadron_MC.Add(Form("%s/HadronTree_MC", filePath.c_str()));
            chainHadron_Reco.Add(Form("%s/HadronTree_RECO", filePath.c_str()));
            fileCount++;
        }
    }
    if (fileCount == 0) {
        cerr << "Nessun file .root trovato in " << inputDir << endl;
        return;
    }
    // creo un output root 
    const char* outputFile  = "new_plot_2510_epic_pion.root"; 
    TFile outFile(outputFile, "RECREATE");  // File di output ROOT

    // To save all the variables
    // here we collect all the variables from the ttree
    // Electron
    chainElectron.SetBranchAddress("el_px_mc", &electron_px);
    chainElectron.SetBranchAddress("el_py_mc", &electron_py);
    chainElectron.SetBranchAddress("el_pz_mc", &electron_pz);
    chainElectron.SetBranchAddress("el_mom_mc", &electron_mom);
    chainElectron.SetBranchAddress("el_theta_mc", &electron_Theta);
    chainElectron.SetBranchAddress("el_phi_mc", &electron_Phi);
    chainElectron.SetBranchAddress("el_eta_mc", &electron_eta);
    chainElectron.SetBranchAddress("el_y_mc", &electron_y);
    // Pion +
    chainHadron_Reco.SetBranchAddress("hadron_index", &pionp_mc_index);
    chainHadron_Reco.SetBranchAddress("hadron_pdg", &rec_pdg);
    chainHadron_Reco.SetBranchAddress("hadron_good_PID", &good_PID);
    chainHadron_Reco.SetBranchAddress("hadron_px", &pionp_px);
    chainHadron_Reco.SetBranchAddress("hadron_py", &pionp_py);
    chainHadron_Reco.SetBranchAddress("hadron_pz", &pionp_pz);
    chainHadron_Reco.SetBranchAddress("hadron_mom", &pionp_mom);
    //chainHadron_Reco.SetBranchAddress("W", &pionp_W);
    chainHadron_Reco.SetBranchAddress("hadron_Q2", &pionp_Q2);
    //chainHadron_Reco.SetBranchAddress("hadron_xF", &pionp_xF);
    chainHadron_Reco.SetBranchAddress("hadron_xB", &pionp_xB);
    chainHadron_Reco.SetBranchAddress("hadron_y", &pionp_y);
    chainHadron_Reco.SetBranchAddress("hadron_z", &pionp_z);
    chainHadron_Reco.SetBranchAddress("hadron_PhT", &pionp_PhT);
    chainHadron_Reco.SetBranchAddress("hadron_Phi_lab", &pionp_Phi_lab);
    chainHadron_Reco.SetBranchAddress("hadron_Theta", &pionp_Theta);
    chainHadron_Reco.SetBranchAddress("hadron_eta", &pionp_eta);
    chainHadron_Reco.SetBranchAddress("hadron_Phi_h", &pionp_Phi_h);
    //chainHadron_Reco.SetBranchAddress("Phi_s", &pionp_Phi_s);
    //chainHadron_Reco.SetBranchAddress("helicity", &helicity);
    //chainHadron_Reco.SetBranchAddress("hadron_Mx", &pionp_Mx);
    //
    // PionMC +
    chainHadron_MC.SetBranchAddress("hadron_index", &index_mc);
    chainHadron_MC.SetBranchAddress("hadron_mom_mc", &pionp_mom_mc);
    chainHadron_MC.SetBranchAddress("hadron_Q2_mc", &pionp_Q2_mc);
    chainHadron_MC.SetBranchAddress("hadron_xB_mc", &pionp_xB_mc);
    //chainHadron_MC.SetBranchAddress("hadron_xF_mc", &pionp_xF_mc);
    chainHadron_MC.SetBranchAddress("hadron_z_mc", &pionp_z_mc);
    chainHadron_MC.SetBranchAddress("hadron_PhT_mc", &pionp_PhT_mc);
    chainHadron_MC.SetBranchAddress("hadron_Phi_lab_mc", &pionp_Phi_lab_mc);
    chainHadron_MC.SetBranchAddress("hadron_Phi_h_mc", &pionp_Phi_h_mc);
    //chainHadron_MC.SetBranchAddress("Phi_s_mc", &pionp_Phi_s_mc);
    chainHadron_MC.SetBranchAddress("hadron_Theta_mc", &pionp_Theta_mc);
    chainHadron_MC.SetBranchAddress("hadron_eta_mc", &pionp_eta_mc);
    chainHadron_MC.SetBranchAddress("hadron_y_mc", &pionp_y_mc);
    //chainHadron_MC.SetBranchAddress("hadron_W_mc", &pionp_W_mc);
    //chainHadron_MC.SetBranchAddress("hadron_Mx_mc", &pionp_Mx_mc);
    //chainHadron_MC.SetBranchAddress("helicity_mc", &hel_mc);
    //chainHadron_MC.SetBranchAddress("hadron_epsilon_mc", &eps_mc);
    chainHadron_MC.SetBranchAddress("hadron_px_mc", &pionp_px_mc);
    chainHadron_MC.SetBranchAddress("hadron_py_mc", &pionp_py_mc);
    chainHadron_MC.SetBranchAddress("hadron_pz_mc", &pionp_pz_mc);
    // tree
    // Pion +
    treePionP.Branch("mc_index", &pionp_mc_index, "mc_index/I");
    treePionP.Branch("rec_pdg", &rec_pdg, "rec_pdg/D");
    treePionP.Branch("good_PID", &good_PID, "good_PID/D");
    treePionP.Branch("px", &pionp_px, "pionp_px/D");
    treePionP.Branch("py", &pionp_py, "pionp_py/D");
    treePionP.Branch("pz", &pionp_pz, "pionp_pz/D");
    treePionP.Branch("E", &pionp_E, "E/D");
    treePionP.Branch("Mom", &pionp_mom, "Mom/D");
    treePionP.Branch("Q2", &pionp_Q2, "Q2/D");
    treePionP.Branch("xB", &pionp_xB, "xB/D");
    treePionP.Branch("xF", &pionp_xF, "xF/D");
    treePionP.Branch("z", &pionp_z, "z/D");
    treePionP.Branch("PhT", &pionp_PhT, "PhT/D");
    treePionP.Branch("Phi_lab", &pionp_Phi_lab, "Phi_h/D");
    treePionP.Branch("Phi_h", &pionp_Phi_h, "Phi_h/D");
    treePionP.Branch("Phi_s", &pionp_Phi_s, "Phi_s/D");
    treePionP.Branch("theta", &pionp_Theta, "theta/D");
    treePionP.Branch("eta", &pionp_eta, "eta/D");
    treePionP.Branch("y", &pionp_y, "y/D");
    treePionP.Branch("W", &pionp_W, "W/D");
    treePionP.Branch("Mx", &pionp_Mx, "Mx/D");
    treePionP.Branch("helicity", &helicity, "hel/D");
    treePionP.Branch("epsilon", &eps, "eps/D");
    //
    // PionMC +
    treePionP_MC.Branch("index", &index_mc, "index/I");
    treePionP_MC.Branch("Mom_mc", &pionp_mom_mc, "Mom_mc/D");
    treePionP_MC.Branch("Q2_mc", &pionp_Q2_mc, "Q2_mc/D");
    treePionP_MC.Branch("xB_mc", &pionp_xB_mc, "xB_mc/D");
    treePionP_MC.Branch("xF_mc", &pionp_xF_mc, "xF_mc/D");
    treePionP_MC.Branch("z_mc", &pionp_z_mc, "z_mc/D");
    treePionP_MC.Branch("PhT_mc", &pionp_PhT_mc, "PhT_mc/D");
    treePionP_MC.Branch("Phi_lab_mc", &pionp_Phi_lab_mc, "Phi_lab_mc/D");
    treePionP_MC.Branch("Phi_h_mc", &pionp_Phi_h_mc, "Phi_h_mc/D");
    treePionP_MC.Branch("Phi_s_mc", &pionp_Phi_s_mc, "Phi_s_mc/D");
    treePionP_MC.Branch("theta_mc", &pionp_Theta_mc, "theta_mc/D");
    treePionP_MC.Branch("eta_mc", &pionp_eta_mc, "eta_mc/D");
    treePionP_MC.Branch("y_mc", &pionp_y_mc, "y_mc/D");
    treePionP_MC.Branch("W_mc", &pionp_W_mc, "W_mc/D");
    treePionP_MC.Branch("Mx_mc", &pionp_Mx_mc, "Mx_mc/D");
    treePionP_MC.Branch("helicity_mc", &hel_mc, "helicity_mc/D");
    treePionP_MC.Branch("epsilon_mc", &eps_mc, "epsilon_mc/D");
    treePionP_MC.Branch("pion_px_mc", &pionp_px_mc, "pion_px_mc/D");
    treePionP_MC.Branch("pion_py_mc", &pionp_py_mc, "pion_py_mc/D");
    treePionP_MC.Branch("pion_pz_mc", &pionp_pz_mc, "pion_pz_mc/D");
    //  
    /*
    treePionP_all.Branch("index", &index_all, "index/I");
    treePionP_all.Branch("Mom_all", &pionp_mom_all, "Mom_all/D");
    treePionP_all.Branch("Q2_all", &pionp_Q2_all, "Q2_all/D");
    treePionP_all.Branch("xB_all", &pionp_xB_all, "xB_all/D");
    treePionP_all.Branch("xF_all", &pionp_xF_all, "xF_all/D");
    treePionP_all.Branch("z_all", &pionp_z_all, "z_all/D");
    treePionP_all.Branch("PhT_all", &pionp_PhT_all, "PhT_all/D");
    treePionP_all.Branch("Phi_lab_all", &pionp_Phi_lab_all, "Phi_lab_all/D");
    treePionP_all.Branch("Phi_h_all", &pionp_Phi_h_all, "Phi_h_all/D");
    treePionP_all.Branch("Phi_s_all", &pionp_Phi_s_all, "Phi_s_all/D");
    treePionP_all.Branch("theta_all", &pionp_Theta_all, "theta_all/D");
    treePionP_all.Branch("eta_all", &pionp_eta_all, "eta_all/D");
    treePionP_all.Branch("y_all", &pionp_y_all, "y_all/D");
    treePionP_all.Branch("W_all", &pionp_W_all, "W_all/D");
    treePionP_all.Branch("Mx_all", &pionp_Mx_all, "Mx_all/D");
    treePionP_all.Branch("helicity_all", &hel_all, "helicity_all/D");
    treePionP_all.Branch("epsilon_all", &eps_all, "epsilon_all/D");
    treePionP_all.Branch("pion_px_all", &pionp_px_all, "pion_px_all/D");
    treePionP_all.Branch("pion_py_all", &pionp_py_all, "pion_py_all/D");
    treePionP_all.Branch("pion_pz_all", &pionp_pz_all, "pion_pz_all/D");
    */

    //
    // PLOT
    //dirPionp->cd();
    double bin = 150;
    double chi_min = 0.1;
    double chi_max = 2000;
    const double xmin_xB = 1e-4, xmax_xB = 1;
    const double xmin_Q2 = 1, xmax_Q2 = 100.;
    auto make_bins = [](int bins, double min, double max) {
        return CreateLogBinning(bins, min, max);
      };
    const auto log_chi2 = make_bins(bin, chi_min, chi_max);
    const auto log_bins_Q2 = make_bins(bin, xmin_Q2, xmax_Q2);
    const auto log_bins_xB = make_bins(bin, xmin_xB, xmax_xB);
    TH2D pip_Q2VsXb_MC ("_Q2VsXb_MC", "Correlation Q^{2} vs x_{B}  |  MC #pi+ ; x_{B}; Q^{2} [GeV^{2}]", bin, log_bins_xB.data(), bin, log_bins_Q2.data());
    TH2D pip_PhTvsZ_MC ("_PhTvsZ_MC", "Correlation P_{hT} vs Z  |  MC #pi+ ; z; P_{hT} [GeV]", bin, 0, 1, bin, 0, 2);
    // Mom
    TH1D pip_evnt_chi2 ("_evnt_chi2", "#chi^{2} EventBuilder PID | only EventBuilder | 1.2 < Mom < 8 GeV ; #chi^{2}; count", bin, -8, 8);
    TH1D pip_m ("_best_mass", "m extracted from #beta | 1.2 < Mom < 8 GeV ; m [GeV]; count", bin, 0, 1);
    TH1D pip_deltaB ("_delta_beta", "#beta_{meas} - #beta_{th} | 1.2 < Mom < 8 GeV ; #Delta_{#beta}; count", bin, -0.05, 0.05);
    TH1D pip_Mx ("_missing_mass", "missing mass | 1.2 < Mom < 8 GeV ; M_{x} [GeV]; count", bin, 0, 60);
    //TH1D pip_rich_chi2 ("_rich_chi2", "#chi^{2} RICH PID ; #chi^{2}; count", 300, -3, 3);
    TH2D pip_MomVsPhT ("_MomVsPhT", "Correlation Mom vs P_{hT}  |  #pi+ ; P_{hT} [GeV]; Mom [GeV]", bin, 0, 2, bin, 0, 20);
    TH2D pip_MomVsXb ("_MomVsXb", "Correlation Mom vs x_{B}  |  #pi+ ; x_{B}; Mom [GeV]", bin, log_bins_xB.data(), bin, 0, 20);
    TH2D pip_MomVsXf ("_MomVsXf", "Correlation Mom vs x_{F}  |  #pi+ ; x_{F}; Mom [GeV]", bin, -0.5, 0.5, bin, 0, 20);
    TH2D pip_MomVsZ ("_MomVsZ", "Correlation Mom vs Z  |  #pi+ ; z; Mom [GeV]", bin, 0, 1, bin, 0, 20);
    TH2D pip_MomVsY ("_MomVsY", "Correlation Mom vs Y  |  #pi+ ; y; Mom [GeV]", bin, 0.0, 1.0, bin, 0, 20);
    TH2D pip_MomVsEta ("_MomVsEta", "Correlation Mom vs Eta  |  #pi+ ; Eta; Mom [GeV]", bin, -3, 3.5, bin, 0, 20);
    TH2D pip_MomVsTheta ("_MomVsTheta", "Correlation Mom vs Theta  |  #pi+ ; Mom [GeV]; #theta [Rad]", bin, 0, 20, bin, 0, 180);
    TH2D pip_MomVsPhi_h ("_MomVsPhi_h", "Correlation Mom vs #Phi_{h}  |  #pi+ ; #Phi_{h} [Rad]; Mom [GeV]", bin, -TMath::Pi(), TMath::Pi(), bin, 0, 20);
    TH2D pip_MomVsMx ("_MomVsMx", "correlation Mom vs M^{2}_{x} | #pi+ |; Mom [GeV]; M_{x} [GeV]", bin, 0, 20, bin, 0, 62);
    // Q2
    TH2D pip_Q2VsXb ("_Q2VsXb", "Correlation Q^{2} vs x_{B}  |  #pi+ ; x_{B}; Q^{2} [GeV^{2}]", bin, log_bins_xB.data(), bin, log_bins_Q2.data());
    TH2D pip_Q2VsXf ("_Q2VsXf", "Correlation Q^{2} vs x_{F}  |  #pi+ ; x_{F}; Q^{2} [GeV^{2}]", bin, -0.5, 0.5, bin, log_bins_Q2.data());
    TH2D pip_Q2VsMom ("_Q2VsMom", "Correlation Q^{2} vs Mom  |  #pi+ ; Mom [GeV]; Q^{2} [GeV^{2}]", bin, 0, 20, bin, log_bins_Q2.data());
    TH2D pip_Q2VsPhT ("_Q2VsPhT", "Correlation Q^{2} vs P_{hT}  |  #pi+ ; P_{hT} [GeV]; Q^{2} [GeV^{2}]", bin, 0, 2, bin, log_bins_Q2.data());
    TH2D pip_Q2VsZ ("_Q2VsZ", "Correlation Q^{2} vs Z  |  #pi+ ; z; Q^{2} [GeV^{2}]", bin, 0, 1, bin, log_bins_Q2.data());
    TH2D pip_Q2VsY ("_Q2VsY", "Correlation Q^{2} vs Y  |  #pi+ ; y; Q^{2} [GeV^{2}]", bin, 0.0, 1.0, bin, log_bins_Q2.data());
    TH2D pip_Q2VsEta ("_Q2VsEta", "Correlation Q^{2} vs Eta  |  #pi+ ; Eta; Q^{2} [GeV^{2}]", bin, -3, 3.5, bin, log_bins_Q2.data());
    TH2D pip_Q2VsPhi_h ("_Q2VsPhi_h", "Correlation Q^{2} vs #Phi_{h}  |  #pi+ ; #Phi_{h} [Rad]; Q^{2} [GeV^{2}]", bin, -TMath::Pi(), TMath::Pi(), bin, log_bins_Q2.data());
    // PhT
    TH2D pip_PhTvsZ ("_PhTvsZ", "Correlation P_{hT} vs Z  |  #pi+ ; z; P_{hT} [GeV]", bin, 0, 1, bin, 0, 2);
    TH2D pip_PhTvsXb ("_PhTvsXb", "Correlation P_{hT} vs x_{B}  |  #pi+ ; x_{B}; P_{hT} [GeV]", bin, log_bins_xB.data(), bin, 0, 2);
    TH2D pip_PhTvsEta ("_PhTvsEta", "Correlation P_{hT} vs Eta  |  #pi+ ; Eta; P_{hT} [GeV]", bin, -3, 3.5, bin, 0, 2);
    TH2D pip_PhTvsPhi_h ("_PhTvsPhi_h", "Correlation P_{hT} vs #Phi_{h}  |  #pi+ ; #Phi_{h} [Rad]; P_{hT} [GeV]", bin, -TMath::Pi(), TMath::Pi(), bin, 0, 2);
    // Z
    TH2D pip_zVsXb ("_zVsXb", "Correlation Z vs x_{B}  |  #pi+ ; x_{B}; z", bin, log_bins_xB.data(), bin, 0, 1);
    TH2D pip_zVsXf ("_zVsXf", "Correlation Z vs x_{F}  |  #pi+ ; x_{F}; z", bin, -0.5, 0.5, bin, 0, 1);
    TH2D pip_zVsEta ("_zVsEta", "Correlation Z vs Eta  |  #pi+ ; Eta; z", bin, -3, 3.5, bin, 0, 1);
    TH2D pip_zVsPhi_h ("_zVsPhi_h", "Correlation Z vs #Phi_{h}  |  #pi+ ; #Phi_{h} [Rad]; z", bin, -TMath::Pi(), TMath::Pi(), bin, 0, 1);
    //
    TH2D pip_xBvsY ("_xBvsY", "Correlation y vs x_{B}  |  #pi+ ; x_{B}; y", bin, log_bins_xB.data(), bin, 0.0, 1.0);
    // Angles
    TH2D pip_ThetaVsPhi_h ("_ThetaVsPhi_h", "Correlation Theta vs #Phi_{h}  |  #pi+ ; #Phi_{h} [Rad]; #theta [Rad]", bin, -TMath::Pi(), TMath::Pi(), bin, 0, 180);
    TH2D pip_ThetaVsPhi_Lab ("_ThetaVsPhi_Lab", "Correlation #theta vs #Phi_{Lab} | #pi+; #Phi_{Lab} [Rad]; #theta [Rad]", bin, -TMath::Pi(), TMath::Pi(), bin, 0, 180);
    //



    // ARRAY 2D AND 4D FOR THE EFFICIENCY AND PURITY (I'LL JUST INSERT THE STATISTICAL UNCERTAINTY INSIDE SINCE THE EFFICIENCY IS NOT CALCULATED YET)
    int nBin_xQ2 = 16;
    int nBin_zPt = 30;
    double bin_xB[][2] = {{1e-4, 1e-3}, {1e-3, 3e-3}, {3e-3, 1e-2}, {1e-2, 4e-2}, {4e-2, 1}};
    // bin_Q2 generali
    double bin_Q2_full[][2] = {{1, 2}, {2, 5}, {5, 12}, {12, 100}};
    double bin_xB_plot[] = {1e-4, 1e-3, 3e-3, 1e-2, 4e-2, 1};
    double bin_Q2_plot[] = {1, 2, 5, 12, 100};
    //
    double bin_z[][2] = {{0, 0.05}, {0.05, 0.1}, {0.1, 0.2}, {0.2, 0.3}, {0.3, 0.4}, {0.4, 1}};
    double bin_z_plot[] = {0, 0.05, 0.1, 0.2, 0.3, 0.4, 1};
    double bin_Pt_plot[] = {0, 0.2, 0.4, 0.6, 0.8, 2};
    // bin_Q2 
    vector<vector<array<double,2>>> binning_Pt_for_z = {
        {{0, 0.2}, {0.2, 0.4}, {0.4, 0.6}, {0.6, 0.8}, {0.8, 2}},
        {{0, 0.2}, {0.2, 0.4}, {0.4, 0.6}, {0.6, 0.8}, {0.8, 2}},
        {{0, 0.2}, {0.2, 0.4}, {0.4, 0.6}, {0.6, 0.8}, {0.8, 2}},
        {{0, 0.2}, {0.2, 0.4}, {0.4, 0.6}, {0.6, 0.8}, {0.8, 2}},
        {{0, 0.2}, {0.2, 0.4}, {0.4, 0.6}, {0.6, 0.8}, {0.8, 2}},
        {{0, 0.2}, {0.2, 0.4}, {0.4, 0.6}, {0.6, 0.8}, {0.8, 2}}
    };
    // struttura asimmetrica dei bin
    vector<vector<array<double,2>>> binning_Q2_for_xB = {
        // in sequenza i bin di Q2 per i rispettivi bin di xB
        {{1,2}, {2,100}},
        {{1,2}, {2,5}, {5,100}},
        {{1,2}, {2,5}, {5,12}, {12,100}},
        {{1,2}, {2,5}, {5,12}, {12,100}},
        {{1,5}, {5,12}, {12,100}}
    };
    vector<vector<double>> pip_efficiency_xQ2(nBin_xQ2);
    vector<vector<double>> pip_efficiency_zPt(nBin_zPt);
    vector<vector<vector<double>>> pip_efficiency_xQ2_zPt(nBin_xQ2, vector<vector<double>> (nBin_zPt));
    vector<vector<vector<double>>> pip_efficiency_xQ2_zPt_2(nBin_xQ2, vector<vector<double>> (nBin_zPt));
    vector<vector<vector<double>>> pip_efficiency_xQ2_zPt_mc(nBin_xQ2, vector<vector<double>> (nBin_zPt));
    vector<vector<vector<double>>> pip_purity_xQ2_zPt_num(nBin_xQ2, vector<vector<double>> (nBin_zPt));
    vector<vector<vector<double>>> pip_purity_xQ2_zPt_den(nBin_xQ2, vector<vector<double>> (nBin_zPt));
    
    // relative plot
    //vector<TH1D*> hist_efficiency_xQ2(nBin_xQ2);
    //vector<TH1D*> hist_efficiency_zPt(nBin_zPt);
    //TH2D hist_efficiency_xQ2 = new TH2D("hist_efficiency_xQ2", "Efficiency pionp in x_{B},Q^{2} bin", 5, bin_xB_plot, 4, bin_Q2_plot);
    vector<TH2D*> hist_efficiency_xQ2_zPt(nBin_xQ2);
    vector<TH2D*> hist_purity_xQ2_zPt(nBin_xQ2);
    for (int ix = 0; ix < nBin_xQ2; ix++){
        hist_efficiency_xQ2_zPt[ix] = new TH2D(Form("hist_efficiency4D_xQ2_%d", ix+1),
            Form("#pi+ efficiency as P_{hT} vs z for bin (%d, x_{B}-Q^{2}); z; P_{hT} [GeV]", ix+1), 6, bin_z_plot, 5, bin_Pt_plot);
        hist_purity_xQ2_zPt[ix] = new TH2D(Form("hist_purity4D_xQ2_%d", ix+1),
            Form("#pi+ purity as P_{hT} vs z for bin (%d, x_{B}-Q^{2}); z; P_{hT} [GeV]", ix+1), 6, bin_z_plot, 5, bin_Pt_plot);
    }

    //


    // ORA RIEMPI I GRAFICI
    // Pion+
    Long64_t nEntries_pip = chainHadron_Reco.GetEntries();
    Long64_t nEntries_pipMC = chainHadron_MC.GetEntries();
    //
    for (Long64_t i = 0; i < nEntries_pipMC; i++) {
        chainHadron_MC.GetEntry(i);
        if (i % 100000 == 0) cout << "MC entry: " << i << "/" << nEntries_pipMC << endl;
        if(pionp_y_mc < 0.99 && pionp_y_mc >= 0.01){
            double bin_xQ2 = getBinIndex_xQ2(pionp_xB_mc, pionp_Q2_mc);
            double bin_zPt = getBinIndex_zPt(pionp_z_mc, pionp_PhT_mc);
            //double bin_z = getBinIndex_z(pionp_z);
            //double bin_Pt = getBinIndex_Pt(pionp_PhT);
            //
            if(bin_xQ2 >= 0){
                if(bin_zPt >= 0){
                    pip_efficiency_xQ2_zPt_mc[bin_xQ2-1][bin_zPt-1].push_back(pionp_y_mc); // compenso il fatto che i bin partano da 1
                }
            }
            pip_Q2VsXb_MC.Fill(pionp_xB_mc, pionp_Q2_mc);
            pip_PhTvsZ_MC.Fill(pionp_z_mc, pionp_PhT_mc);
            treePionP_MC.Fill();
        }
    }
    //
    /*
    for (Long64_t i = 0; i < nEntries_pipAll; i++) {
        chainPionP_all.GetEntry(i);
        if(good_PID_all != 0) continue;
        if(pionp_y_all < 0.99 && pionp_y_all >= 0.01){
            double bin_xQ2 = getBinIndex_xQ2(pionp_xB_all, pionp_Q2_all);
            double bin_zPt = getBinIndex_zPt(pionp_z_all, pionp_PhT_all);
            //
            if(bin_xQ2 >= 0){
                if(bin_zPt >= 0){
                    pip_purity_xQ2_zPt_den[bin_xQ2-1][bin_zPt-1].push_back(pionp_y_all);
                    if(pdg_all == 211) pip_purity_xQ2_zPt_num[bin_xQ2-1][bin_zPt-1].push_back(pionp_y_all); // compenso il fatto che i bin partano da 1
                }
            }
            //treePionP_all.Fill();
        }
    }
        */
    //
    for (Long64_t i = 0; i < nEntries_pip; i++) {
        chainHadron_Reco.GetEntry(i);
        if (i % 100000 == 0) cout << "RECO entry: " << i << "/" << nEntries_pip << endl;
        //if(rec_pdg == 211 && pionp_z < 1 && pionp_Q2 >= 1){
            if(pionp_y < 0.99 && pionp_y >= 0.01){
                double bin_xQ2 = getBinIndex_xQ2(pionp_xB, pionp_Q2);
                double bin_zPt = getBinIndex_zPt(pionp_z, pionp_PhT);
                //if(bin_xQ2 != 14) continue;
                //
                if(bin_xQ2 >= 0){
                    if(bin_zPt >= 0){
                        pip_efficiency_xQ2_zPt_2[bin_xQ2-1][bin_zPt-1].push_back(pionp_y); // compenso il fatto che i bin partano da 1
                    }
                }
            }
            if(pionp_y < 0.99 && pionp_y >= 0.01 && good_PID == 0 && rec_pdg == 211 && pionp_z < 1 && pionp_Q2 >= 1){
                double bin_xQ2 = getBinIndex_xQ2(pionp_xB, pionp_Q2);
                double bin_zPt = getBinIndex_zPt(pionp_z, pionp_PhT);
                double bin_z = getBinIndex_z(pionp_z);
                double bin_Pt = getBinIndex_Pt(pionp_PhT);
                //if(bin_xQ2 != 14) continue;
                //
                if(bin_xQ2 >= 0){
                    if(bin_zPt >= 0){
                        pip_efficiency_xQ2_zPt[bin_xQ2-1][bin_zPt-1].push_back(pionp_y); 
                    }
                }
                //
                // Mom
                pip_MomVsPhT.Fill(pionp_PhT, pionp_mom);
                pip_MomVsEta.Fill(pionp_eta, pionp_mom);
                pip_MomVsMx.Fill(pionp_mom, pionp_Mx);
                pip_MomVsPhi_h.Fill(pionp_Phi_h, pionp_mom);
                pip_MomVsTheta.Fill(pionp_mom, pionp_Theta);
                pip_MomVsXb.Fill(pionp_xB, pionp_mom);
                pip_MomVsXf.Fill(pionp_xF, pionp_mom);
                pip_MomVsY.Fill(pionp_y, pionp_mom);
                pip_MomVsZ.Fill(pionp_z, pionp_mom);
                // Q2
                pip_Q2VsEta.Fill(pionp_eta, pionp_Q2);
                pip_Q2VsMom.Fill(pionp_mom, pionp_Q2);
                pip_Q2VsPhi_h.Fill(pionp_Phi_h, pionp_Q2);
                pip_Q2VsPhT.Fill(pionp_PhT, pionp_Q2);
                pip_Q2VsXb.Fill(pionp_xB, pionp_Q2);
                pip_Q2VsXf.Fill(pionp_xF, pionp_Q2);
                pip_Q2VsY.Fill(pionp_y, pionp_Q2);
                pip_Q2VsZ.Fill(pionp_z, pionp_Q2);
                // PhT
                pip_PhTvsEta.Fill(pionp_eta, pionp_PhT);
                pip_PhTvsXb.Fill(pionp_xB, pionp_PhT);
                pip_PhTvsZ.Fill(pionp_z, pionp_PhT);
                pip_PhTvsPhi_h.Fill(pionp_Phi_h, pionp_PhT);
                // z
                pip_zVsEta.Fill(pionp_eta, pionp_z);
                pip_zVsPhi_h.Fill(pionp_Phi_h, pionp_z);
                pip_zVsXb.Fill(pionp_xB, pionp_z);
                pip_zVsXf.Fill(pionp_xF, pionp_z);
                //
                pip_xBvsY.Fill(pionp_xB, pionp_y);
                pip_ThetaVsPhi_h.Fill(pionp_Phi_h, pionp_Theta);
                pip_ThetaVsPhi_Lab.Fill(pionp_Phi_lab, pionp_Theta);
                treePionP.Fill();
            }

        }

        treePionP.Write();
        treePionP_MC.Write();

        for (int ix = 0; ix < nBin_xQ2; ix++) {
            for (int iz = 0; iz < nBin_zPt; iz++) {
                double n = pip_efficiency_xQ2_zPt[ix][iz].size();
                double n_mc = pip_efficiency_xQ2_zPt_mc[ix][iz].size();
                double n2 = pip_efficiency_xQ2_zPt_2[ix][iz].size();
                double n_all = pip_purity_xQ2_zPt_num[ix][iz].size();
                double n_all_den = pip_purity_xQ2_zPt_den[ix][iz].size();
                //double eff = (n > 0) ? 1.0 / sqrt(n) : 0.0;
                double eff = (n > 0 || n2 > 0) ? n/n2 : 0.0;
                double err = (n_mc > 0) ? sqrt(eff * (1.0 - eff) / n_mc) : 0.0;
                double purity = (n > 0 || n_all > 0) ? n/n_all_den : 0.0;
                //if(n < 300) eff = 0;
                //hist_efficiency_xQ2_zPt[ix]->SetBinError(..., err);
                if (eff >= 1){
                    //cout << "eff > 1 in bin xQ2: " << ix+1 << " and bin zPt: " << iz+1 << ". Eff = " << eff << endl;
                    eff = 1.0;
                }
                if (purity >= 1){
                    //cout << "purity > 1 in bin xQ2: " << ix+1 << " and bin zPt: " << iz+1 << ". Purity = " << purity << endl;
                    purity = 1.0;
                }
                if(iz < 5) {
                    hist_efficiency_xQ2_zPt[ix]->SetBinContent(1, iz+1, eff);
                    hist_purity_xQ2_zPt[ix]->SetBinContent(1, iz+1, purity);
                } else if (iz < 10) {
                    hist_efficiency_xQ2_zPt[ix]->SetBinContent(2, iz+1-5, eff);
                    hist_purity_xQ2_zPt[ix]->SetBinContent(2, iz+1-5, purity);
                } else if (iz < 15) {
                    hist_efficiency_xQ2_zPt[ix]->SetBinContent(3, iz+1-10, eff);
                    hist_purity_xQ2_zPt[ix]->SetBinContent(3, iz+1-10, purity);
                } else if (iz < 20) {
                    hist_efficiency_xQ2_zPt[ix]->SetBinContent(4, iz+1-15, eff);
                    hist_purity_xQ2_zPt[ix]->SetBinContent(4, iz+1-15, purity);
                } else if (iz < 25) {
                    hist_efficiency_xQ2_zPt[ix]->SetBinContent(5, iz+1-20, eff);
                    hist_purity_xQ2_zPt[ix]->SetBinContent(5, iz+1-20, purity);
                } else if (iz < 30) {
                    hist_efficiency_xQ2_zPt[ix]->SetBinContent(6, iz+1-25, eff);
                    hist_purity_xQ2_zPt[ix]->SetBinContent(6, iz+1-25, purity);
                }
            }
        //}
    }


    vector<TH2D*> hists_pip = {
        &pip_MomVsPhT, &pip_MomVsXb, &pip_MomVsXf, &pip_MomVsZ, &pip_MomVsY, &pip_MomVsEta,
        &pip_MomVsTheta, &pip_MomVsPhi_h, &pip_MomVsMx,  
        &pip_Q2VsXb, &pip_Q2VsXb_MC, &pip_Q2VsXf, &pip_Q2VsMom, &pip_Q2VsPhT, &pip_Q2VsZ, &pip_Q2VsY, &pip_Q2VsEta, &pip_Q2VsPhi_h,
        &pip_PhTvsZ, &pip_PhTvsZ_MC, &pip_PhTvsXb, &pip_PhTvsEta, &pip_PhTvsPhi_h,
        &pip_zVsXb, &pip_zVsXf, &pip_zVsEta, &pip_zVsPhi_h,
        &pip_xBvsY, &pip_ThetaVsPhi_h, &pip_ThetaVsPhi_Lab,
    };

    // list to set the statbox2 on the canvas
    set<string> id_box2 = { "_MomVsXf", "_MomVsZ", "_MomVsY", "_MomVsEta", "_Q2VsXb", "_Q2VsY", "_zVsEta", "_zVsPhi_h", "_MomVsMass_RICH"};
    set<string> id_box3 = { "_MomVsXb", "_Q2VsXb", "_Q2VsXb_MC", "_PhTvsXb", "_zVsXb", "_xBvsY"};
    set<string> id_box4 = { "_Q2VsXb", "_Q2VsXb_MC", "_Q2VsXf", "_Q2VsMom", "_Q2VsPhT", "_Q2VsZ", "_Q2VsY", "_Q2VsEta", "_Q2VsPhi_h"};
    for (size_t i = 0; i < hists_pip.size(); ++i) {
        const char* histName = hists_pip[i]->GetName();
        //string cname_pip = string("c_pip") + histName;
        TCanvas *c = new TCanvas(histName, histName, 800, 600);
        c->SetLogz();   
        hists_pip[i]->SetStats(0);
        hists_pip[i]->SetTitle("");
        hists_pip[i]->Draw("");
        //gPad->Update();
        if (id_box3.count(histName)) c->SetLogx();
        if (id_box4.count(histName)) c->SetLogy();
        //if (id_box2.count(histName)) SetStatsBox2(hists_pip[i]);
        //else SetStatsBox(hists_pip[i]);
        TLatex Text_ePIC;
        Text_ePIC.SetTextSize(0.05);
        Text_ePIC.SetTextFont(62);
        //Text_ePIC.DrawLatexNDC(.15,.88,"ePIC Performance");  // performance plot
        //Text_ePIC.DrawLatexNDC(.15,.88,"ePIC Internal");  // for internal use only
        Text_ePIC.DrawLatexNDC(.15,.88,"ePIC Preliminary"); // preliminary released version 
        //Text_ePIC.DrawLatexNDC(.15,.88,"ePIC Work in Progress"); // work in progress to be shown outside
        //Text_ePIC.DrawLatexNDC(.15,.88,"ePIC"); // final published version
        TLatex Text_com;
        Text_com.SetTextAlign(13);  //align at top
        Text_com.DrawLatexNDC(.15,.85,"e+p, #sqrt{s} = 63.2 GeV");
        //Text_com.DrawLatexNDC(.15,.8,"L_{proj} = 10 fb^{-1}");
        TLatex Text_date;
        Text_date.SetTextSize(0.035);
        Text_date.SetTextFont(52);
        Text_date.DrawLatexNDC(.65,.96,"Simu campaign: 10/2025");

        TImage *logo = TImage::Open("EPIC-logo_black_small.png");
        TPad *pad2 = new TPad("pad2", "Pad 2", 0.8, 0.8, 0.93, 0.93); // Create a new pad and then draw the image in it
        pad2->Draw();
        pad2->cd(); // Enter the new pad
        logo->Draw();

        c->Write();
    }
    
    // plot della griglia
    TCanvas *c_bin_xQ2 = new TCanvas("Q2_vs_xB_Bin", "Q^{2} vs x_{B} bin", 800, 700);
    c_bin_xQ2->SetLogz(); c_bin_xQ2->SetLogx(); c_bin_xQ2->SetLogy();
    pip_Q2VsXb.Draw("COLZ");

    vector<TPolyLine*> rectangles;
    vector<TText*> labels;
    int bin_index = 1;
    // Loop sui bin 
    for (int ix = 0; ix < 5; ++ix) {
        for (size_t iq = 0; iq < binning_Q2_for_xB[ix].size(); ++iq) {
            double xB[5] = {bin_xB[ix][0], bin_xB[ix][1], bin_xB[ix][1], bin_xB[ix][0], bin_xB[ix][0]};
            double Q2[5] = {binning_Q2_for_xB[ix][iq][0], binning_Q2_for_xB[ix][iq][0], binning_Q2_for_xB[ix][iq][1], binning_Q2_for_xB[ix][iq][1], binning_Q2_for_xB[ix][iq][0]};

            TPolyLine *rect = new TPolyLine(5, xB, Q2);
            rect->SetLineColor(kBlack);
            rect->SetLineWidth(2);
            rect->Draw("same");
            rectangles.push_back(rect);

            // posizione indici
            double x_center = sqrt(bin_xB[ix][0] * bin_xB[ix][1]);
            double Q2_center = sqrt(binning_Q2_for_xB[ix][iq][0] * binning_Q2_for_xB[ix][iq][1]);

            TText *label = new TText(x_center, Q2_center, Form("%d", bin_index++));
            label->SetTextAlign(22);
            label->SetTextSize(0.03);
            label->SetTextColor(kRed+1); //forse rosso si vede di più
            label->Draw("same");
            labels.push_back(label);
        }
    }

    TLatex Text_ePIC;
    Text_ePIC.SetTextSize(0.05);
    Text_ePIC.SetTextFont(62);
    //Text_ePIC.DrawLatexNDC(.15,.88,"ePIC Performance");  // performance plot
    //Text_ePIC.DrawLatexNDC(.15,.88,"ePIC Internal");  // for internal use only
    Text_ePIC.DrawLatexNDC(.15,.88,"ePIC Preliminary"); // preliminary released version 
    //Text_ePIC.DrawLatexNDC(.15,.88,"ePIC Work in Progress"); // work in progress to be shown outside
    //Text_ePIC.DrawLatexNDC(.15,.88,"ePIC"); // final published version
    TLatex Text_com;
    Text_com.SetTextAlign(13);  //align at top
    Text_com.DrawLatexNDC(.15,.85,"e+p, #sqrt{s} = 63.2 GeV");
    //Text_com.DrawLatexNDC(.15,.8,"L_{proj} = 10 fb^{-1}");
    TLatex Text_date;
    Text_date.SetTextSize(0.035);
    Text_date.SetTextFont(52);
    Text_date.DrawLatexNDC(.65,.96,"Simu campaign: 08/2025");

    TImage *logo = TImage::Open("EPIC-logo_black_small.png");
    TPad *pad2 = new TPad("pad2", "Pad 2", 0.8, 0.8, 0.93, 0.93); // Create a new pad and then draw the image in it
    pad2->Draw();
    pad2->cd(); // Enter the new pad
    logo->Draw();
    c_bin_xQ2->Update();
    c_bin_xQ2->Write();

    TCanvas *c_bin_zPt = new TCanvas("z_vs_Pt_Bin", "z vs P_{hT} bin", 800, 700);
    c_bin_zPt->SetLogz();
    pip_PhTvsZ.Draw("COLZ");

    vector<TPolyLine*> gridLines_zp;
    vector<TText*> labels_zp;

    int bin_index_zp = 1;

    // Loop generale
    for (int iz = 0; iz < 6; ++iz) {
        for (size_t ip = 0; ip < binning_Pt_for_z[iz].size(); ++ip) {

            double z[5] = {bin_z[iz][0], bin_z[iz][1], bin_z[iz][1], bin_z[iz][0], bin_z[iz][0]};
            double Pt[5] = {
                binning_Pt_for_z[iz][ip][0],
                binning_Pt_for_z[iz][ip][0],
                binning_Pt_for_z[iz][ip][1],
                binning_Pt_for_z[iz][ip][1],
                binning_Pt_for_z[iz][ip][0]
            };

            // Rettangolo
            TPolyLine *rect_zp = new TPolyLine(5, z, Pt);
            rect_zp->SetLineWidth(2);
            rect_zp->SetLineColor(kBlack);
            rect_zp->Draw("same");
            gridLines_zp.push_back(rect_zp);

            // Centro geometrico (utile se log-scale in Pt)
            double z_center = 0.5 * (bin_z[iz][0] + bin_z[iz][1]);
            double Pt_center = 0.5 * (binning_Pt_for_z[iz][ip][0] + binning_Pt_for_z[iz][ip][1]);

            TText *label = new TText(z_center, Pt_center, Form("%d", bin_index_zp++));
            label->SetTextAlign(22);
            label->SetTextSize(0.03);
            label->SetTextColor(kRed+1);
            label->Draw("same");
            labels_zp.push_back(label);
        }
    }

    c_bin_zPt->Update();
    c_bin_zPt->Write();


    double global_min = 1e9;
    double global_max = -1e9;
    for (int ix = 0; ix < nBin_xQ2; ++ix) {
        double min_tmp = hist_efficiency_xQ2_zPt[ix]->GetMinimum();
        double max_tmp = hist_efficiency_xQ2_zPt[ix]->GetMaximum();
        if (min_tmp < global_min) global_min = min_tmp;
        if (max_tmp > global_max) global_max = max_tmp;
        //if(global_max > 0.8) global_max = 0.8;
    }

    for (int ix = 0; ix < nBin_xQ2; ++ix) {
        hist_efficiency_xQ2_zPt[ix]->SetMinimum(global_min);
        hist_efficiency_xQ2_zPt[ix]->SetMaximum(global_max);
    }

    // Evita zmin = zmax
    if (global_min == global_max) global_max = global_min + 1e-6;


    for (int ixQ2 = 0; ixQ2 < nBin_xQ2; ++ixQ2) {
        TCanvas *c_bin_zPt = new TCanvas(Form("c_efficiency4D_xQ2_%d", ixQ2+1), Form("Canvas %d", ixQ2+1), 800, 800);
        //c_bin_zPt->SetLogz();
        hist_efficiency_xQ2_zPt[ixQ2]->SetStats(0);
        hist_efficiency_xQ2_zPt[ixQ2]->Draw("colz");

        vector<TPolyLine*> gridLines_zp;
        vector<TText*> labels_zp;

        int bin_index_zp = 1;

        // Loop generale
        for (int iz = 0; iz < 6; ++iz) {
            for (size_t ip = 0; ip < binning_Pt_for_z[iz].size(); ++ip) {

                double z[5] = {bin_z[iz][0], bin_z[iz][1], bin_z[iz][1], bin_z[iz][0], bin_z[iz][0]};
                double Pt[5] = {
                    binning_Pt_for_z[iz][ip][0],
                    binning_Pt_for_z[iz][ip][0],
                    binning_Pt_for_z[iz][ip][1],
                    binning_Pt_for_z[iz][ip][1],
                    binning_Pt_for_z[iz][ip][0]
                };

                // Rettangolo
                TPolyLine *rect_zp = new TPolyLine(5, z, Pt);
                rect_zp->SetLineWidth(2);
                rect_zp->SetLineColor(kBlack);
                rect_zp->Draw("same");
                gridLines_zp.push_back(rect_zp);

                // Centro geometrico (utile se log-scale in Pt)
                double z_center = 0.5 * (bin_z[iz][0] + bin_z[iz][1]);
                double Pt_center = 0.5 * (binning_Pt_for_z[iz][ip][0] + binning_Pt_for_z[iz][ip][1]);

                TText *label = new TText(z_center, Pt_center, Form("%d", bin_index_zp++));
                label->SetTextAlign(22);
                label->SetTextSize(0.03);
                label->SetTextColor(kRed+1);
                label->Draw("same");
                labels_zp.push_back(label);
            }
        }

        c_bin_zPt->Update();
        c_bin_zPt->Write();
    }

    const int nRows = 4;
    const int nCols = 5;    
    int layout[nRows][nCols] = {
        { 0,  0,  9, 13, 16},
        { 0,  5,  8, 12, 15},
        { 2,  4,  7, 11, 14},
        { 1,  3,  6, 10,  0}
    };

    // plot multiplo

    TCanvas *c_layout = new TCanvas("c_4D_efficiency", "All efficiency bins", 1800, 1200);
    c_layout->cd();

    // === OVERLAY PAD per assi log ===
    TPad *pad_axes = new TPad("pad_axes", "Global Axes", 0, 0, 1, 1);
    //pad_axes->SetFillStyle(4000);
    //pad_axes->SetFrameFillStyle(0);
    pad_axes->SetFrameLineWidth(0);
    pad_axes->Draw();
    pad_axes->cd();

    pad_axes->SetTickx(0);
    pad_axes->SetTicky(0);

    // coordinate NDC dove vogliamo gli assi
    double x_ndc_min = 0.05, x_ndc_max = 0.96;
    double y_ndc_min = 0.05, y_ndc_max = 0.96;

    // range fisico (log)
    double x_min = 1e-4, x_max = 1;
    double y_min = 1, y_max = 100;

    // disegna le due linee principali
    TLine *xAxisLine = new TLine(x_ndc_min, y_ndc_min, x_ndc_max, y_ndc_min);
    xAxisLine->SetLineWidth(2);
    xAxisLine->SetNDC(true);
    xAxisLine->Draw();

    TLine *yAxisLine = new TLine(x_ndc_min, y_ndc_min, x_ndc_min, y_ndc_max);
    yAxisLine->SetLineWidth(2);
    yAxisLine->SetNDC(true);
    yAxisLine->Draw();

    // asse x
    for (int i = -4; i <= 0; ++i) {
        double val = pow(10, i);
        double pos = x_ndc_min + (log10(val) - log10(x_min)) / (log10(x_max) - log10(x_min)) * (x_ndc_max - x_ndc_min);

        // tick
        TLine *tick = new TLine(pos, y_ndc_min, pos, y_ndc_min - 0.01);
        tick->SetNDC(true);
        tick->Draw();

        // label
        TLatex *lab = new TLatex(pos, y_ndc_min - 0.03, Form("10^{%d}", i));
        lab->SetTextFont(42);
        lab->SetTextSize(0.02);
        lab->SetTextAlign(22);
        lab->SetNDC(true);
        lab->Draw();
    }

    // titolo asse X
    TLatex *xlabel = new TLatex(0.9, y_ndc_min - 0.025, "x_{B}");
    xlabel->SetTextSize(0.025);
    xlabel->SetTextAlign(22);
    xlabel->SetNDC(true);
    xlabel->Draw();

    // asse y
    for (int i = 0; i <= 2; ++i) {
        double val = pow(10, i);
        double pos = y_ndc_min + (log10(val) - log10(y_min)) / (log10(y_max) - log10(y_min)) * (y_ndc_max - y_ndc_min);

        // tick
        TLine *tick = new TLine(x_ndc_min, pos, x_ndc_min - 0.01, pos);
        tick->SetNDC(true);
        tick->Draw();

        // label
        TLatex *lab = new TLatex(x_ndc_min - 0.01, pos, Form("10^{%d}", i));
        lab->SetTextFont(42);
        lab->SetTextSize(0.02);
        lab->SetTextAlign(32);
        lab->SetNDC(true);
        lab->Draw();
    }

    // titolo asse Y
    TLatex *ylabel = new TLatex(x_ndc_min - 0.025, 0.85, "Q^{2} (GeV^{2})");
    ylabel->SetTextSize(0.025);
    ylabel->SetTextAngle(90);
    ylabel->SetTextAlign(22);
    ylabel->SetNDC(true);
    ylabel->Draw();

    pad_axes->Modified();
    pad_axes->Update();


    c_layout->cd();

    // margini e dimensioni dei pad 
    double xMargin = 0.06, yMargin = 0.06;
    double padW = (1.0 - 2*xMargin) / nCols;
    double padH = (1.0 - 2*yMargin) / nRows;

    std::vector<TPolyLine*> grid_zp;

    for (int iz = 0; iz < 6; ++iz) {
        for (size_t ip = 0; ip < binning_Pt_for_z[iz].size(); ++ip) {
            double z[5] = {bin_z[iz][0], bin_z[iz][1], bin_z[iz][1], bin_z[iz][0], bin_z[iz][0]};
            double Pt[5] = {
                binning_Pt_for_z[iz][ip][0],
                binning_Pt_for_z[iz][ip][0],
                binning_Pt_for_z[iz][ip][1],
                binning_Pt_for_z[iz][ip][1],
                binning_Pt_for_z[iz][ip][0]
            };
            TPolyLine *rect = new TPolyLine(5, z, Pt);
            rect->SetLineWidth(1);
            rect->SetLineColor(kBlack);
            grid_zp.push_back(rect);
        }
    }
    // Loop e creazione pad 
    for (int iRow = 0; iRow < nRows; ++iRow) {
        for (int iCol = 0; iCol < nCols; ++iCol) {
            int idx = layout[iRow][iCol]; // idx è il numero del bin, quindi 0 = vuoto
            if (idx == 0) continue;

            // coordinate normalized: (x1,y1) bottom-left, (x2,y2) top-right
            double x1 = xMargin + iCol * padW;
            double x2 = x1 + padW;
            // ATTENZIONE: iRow = 0 significa riga in alto nella tua matrice, per cui
            // dobbiamo invertire l'ordine per le y in coordinate normalized
            double y2 = 1.0 - yMargin - iRow * padH;
            double y1 = y2 - padH;

            // nome unico del pad 
            TString padName = Form("pad_r%d_c%d", iRow, iCol);
            // ritorno alla canvas
            c_layout->cd();

            TPad *pad = new TPad(padName, padName, x1, y1, x2, y2);
            pad->SetRightMargin(0.0);
            pad->SetLeftMargin(0.0);
            pad->SetBottomMargin(0.0);
            pad->SetTopMargin(0.0);
            pad->Draw();
            pad->cd();

            // Disegna l'istogramma (senza colorbar multipla)
            hist_efficiency_xQ2_zPt[idx - 1]->SetTitle("");
            hist_efficiency_xQ2_zPt[idx - 1]->Draw("col"); // idx-1 perché vettore 0-based
            for (auto &rect : grid_zp) rect->DrawClone("same");
            // Torna al canvas principale per il prossimo pad
            c_layout->cd();
        }
    }
    c_layout->Modified();

    TPad *pad_palette = new TPad("pad_palette", "", 0.94, 0.1, 1, 0.9);
    pad_palette->SetRightMargin(0.5);  // spazio per la palette
    pad_palette->Draw();
    pad_palette->cd();

    // copia istogramma
    TH2F *h_ref = (TH2F*) hist_efficiency_xQ2_zPt[2]->Clone("h_ref_for_palette");
    h_ref->SetStats(0);
    h_ref->Draw("COLZ");
    gPad->Update();

    // recupero la palette
    TPaletteAxis *pal = (TPaletteAxis*) h_ref->GetListOfFunctions()->FindObject("palette");
    if (pal) {
        pal->SetLabelSize(0.22);
        pal->SetTitleOffset(0.0);
        //pal->SetTitle("Purity");
        pal->SetX1NDC(0.05);
        pal->SetX2NDC(0.5);
        pal->SetY1NDC(0.05);
        pal->SetY2NDC(0.95);
    }
    pad_palette->Modified();
    c_layout->Update();

    c_layout->cd();
    // Scritta in alto a sinistra
    TLatex *globalTitle = new TLatex(0.10, 0.95, "4D Efficiency Distribution");
    globalTitle->SetNDC(true);         // coordinate normalizzate (relative alla canvas)
    //globalTitle->SetTextFont(62);    // font bold
    globalTitle->SetTextSize(0.04);    // dimensione
    globalTitle->SetTextAlign(13);     // 13 = left-top aligned
    globalTitle->Draw();
    TLatex *subTitle = new TLatex(0.15, 0.9, "RECO #pi+ / MC #pi+");
    subTitle->SetNDC(true);
    subTitle->SetTextSize(0.03);
    subTitle->SetTextAlign(13);
    subTitle->Draw();
    TLatex *subsub = new TLatex(0.17, 0.86, "(x_{B}, Q^{2}, z, P_{hT})");
    subsub->SetNDC(true);
    subsub->SetTextSize(0.03);
    subsub->SetTextAlign(13);
    subsub->Draw();
    // Fine
    c_layout->Update();
    c_layout->Write();






    // purity 



    double global_min_all = 1e9;
    double global_max_all = -1e9;
    for (int ix = 0; ix < nBin_xQ2; ++ix) {
        double min_tmp = hist_purity_xQ2_zPt[ix]->GetMinimum();
        double max_tmp = hist_purity_xQ2_zPt[ix]->GetMaximum();
        if (min_tmp < global_min_all) global_min_all = min_tmp;
        if (max_tmp > global_max_all) global_max_all = max_tmp;
        if(global_max_all > 1) global_max_all = 1;
        //if(global_min_all < 0.9) global_min_all = 0.9;
    }

    for (int ix = 0; ix < nBin_xQ2; ++ix) {
        hist_purity_xQ2_zPt[ix]->SetMinimum(0.8);
        hist_purity_xQ2_zPt[ix]->SetMaximum(global_max_all);
    }

    // Evita zmin = zmax
    if (global_min_all == global_max_all) global_max_all = global_min_all + 1e-6;


    for (int ixQ2 = 0; ixQ2 < nBin_xQ2; ++ixQ2) {
        TCanvas *c_bin_zPt = new TCanvas(Form("c_purity4D_xQ2_%d", ixQ2+1), Form("Canvas %d", ixQ2+1), 800, 800);
        //c_bin_zPt->SetLogz();
        hist_purity_xQ2_zPt[ixQ2]->SetStats(0);
        hist_purity_xQ2_zPt[ixQ2]->Draw("colz");

        vector<TPolyLine*> gridLines_zp;
        vector<TText*> labels_zp;

        int bin_index_zp = 1;

        // Loop generale
        for (int iz = 0; iz < 6; ++iz) {
            for (size_t ip = 0; ip < binning_Pt_for_z[iz].size(); ++ip) {

                double z[5] = {bin_z[iz][0], bin_z[iz][1], bin_z[iz][1], bin_z[iz][0], bin_z[iz][0]};
                double Pt[5] = {
                    binning_Pt_for_z[iz][ip][0],
                    binning_Pt_for_z[iz][ip][0],
                    binning_Pt_for_z[iz][ip][1],
                    binning_Pt_for_z[iz][ip][1],
                    binning_Pt_for_z[iz][ip][0]
                };

                // Rettangolo
                TPolyLine *rect_zp = new TPolyLine(5, z, Pt);
                rect_zp->SetLineWidth(2);
                rect_zp->SetLineColor(kBlack);
                rect_zp->Draw("same");
                gridLines_zp.push_back(rect_zp);

                // Centro geometrico (utile se log-scale in Pt)
                double z_center = 0.5 * (bin_z[iz][0] + bin_z[iz][1]);
                double Pt_center = 0.5 * (binning_Pt_for_z[iz][ip][0] + binning_Pt_for_z[iz][ip][1]);

                TText *label = new TText(z_center, Pt_center, Form("%d", bin_index_zp++));
                label->SetTextAlign(22);
                label->SetTextSize(0.03);
                label->SetTextColor(kRed+1);
                label->Draw("same");
                labels_zp.push_back(label);
            }
        }

        c_bin_zPt->Update();
        c_bin_zPt->Write();
    }


    TCanvas *c_layout_all = new TCanvas("c_4D_purity", "All purity bins", 1800, 1200);
    c_layout_all->cd();

    pad_axes->Draw();
    pad_axes->cd();

    pad_axes->SetTickx(0);
    pad_axes->SetTicky(0);

    xAxisLine->Draw();
    yAxisLine->Draw();

    // asse x
    for (int i = -4; i <= 0; ++i) {
        double val = pow(10, i);
        double pos = x_ndc_min + (log10(val) - log10(x_min)) / (log10(x_max) - log10(x_min)) * (x_ndc_max - x_ndc_min);

        // tick
        TLine *tick = new TLine(pos, y_ndc_min, pos, y_ndc_min - 0.01);
        tick->SetNDC(true);
        tick->Draw();

        // label
        TLatex *lab = new TLatex(pos, y_ndc_min - 0.03, Form("10^{%d}", i));
        lab->SetTextFont(42);
        lab->SetTextSize(0.02);
        lab->SetTextAlign(22);
        lab->SetNDC(true);
        lab->Draw();
    }

    xlabel->Draw();

    // asse y
    for (int i = 0; i <= 2; ++i) {
        double val = pow(10, i);
        double pos = y_ndc_min + (log10(val) - log10(y_min)) / (log10(y_max) - log10(y_min)) * (y_ndc_max - y_ndc_min);

        // tick
        TLine *tick = new TLine(x_ndc_min, pos, x_ndc_min - 0.01, pos);
        tick->SetNDC(true);
        tick->Draw();

        // label
        TLatex *lab = new TLatex(x_ndc_min - 0.01, pos, Form("10^{%d}", i));
        lab->SetTextFont(42);
        lab->SetTextSize(0.02);
        lab->SetTextAlign(32);
        lab->SetNDC(true);
        lab->Draw();
    }

    ylabel->Draw();

    pad_axes->Modified();
    pad_axes->Update();
    


    c_layout_all->cd();


    for (int iz = 0; iz < 6; ++iz) {
        for (size_t ip = 0; ip < binning_Pt_for_z[iz].size(); ++ip) {
            double z[5] = {bin_z[iz][0], bin_z[iz][1], bin_z[iz][1], bin_z[iz][0], bin_z[iz][0]};
            double Pt[5] = {
                binning_Pt_for_z[iz][ip][0],
                binning_Pt_for_z[iz][ip][0],
                binning_Pt_for_z[iz][ip][1],
                binning_Pt_for_z[iz][ip][1],
                binning_Pt_for_z[iz][ip][0]
            };
            TPolyLine *rect = new TPolyLine(5, z, Pt);
            rect->SetLineWidth(1);
            rect->SetLineColor(kBlack);
            grid_zp.push_back(rect);
        }
    }
    // Loop e creazione pad 
    for (int iRow = 0; iRow < nRows; ++iRow) {
        for (int iCol = 0; iCol < nCols; ++iCol) {
            int idx = layout[iRow][iCol]; // idx è il numero del bin, quindi 0 = vuoto
            if (idx == 0) continue;

            // coordinate normalized: (x1,y1) bottom-left, (x2,y2) top-right
            double x1 = xMargin + iCol * padW;
            double x2 = x1 + padW;
            // ATTENZIONE: iRow = 0 significa riga in alto nella tua matrice, per cui
            // dobbiamo invertire l'ordine per le y in coordinate normalized
            double y2 = 1.0 - yMargin - iRow * padH;
            double y1 = y2 - padH;

            // nome unico del pad 
            TString padName = Form("pad_r%d_c%d", iRow, iCol);
            // ritorno alla canvas
            c_layout_all->cd();

            TPad *pad = new TPad(padName, padName, x1, y1, x2, y2);
            pad->SetRightMargin(0.0);
            pad->SetLeftMargin(0.0);
            pad->SetBottomMargin(0.0);
            pad->SetTopMargin(0.0);
            pad->Draw();
            pad->cd();

            // Disegna l'istogramma (senza colorbar multipla)
            hist_purity_xQ2_zPt[idx - 1]->SetTitle("");
            hist_purity_xQ2_zPt[idx - 1]->Draw("col"); // idx-1 perché vettore 0-based
            for (auto &rect : grid_zp) rect->DrawClone("same");
            // Torna al canvas principale per il prossimo pad
            c_layout_all->cd();
        }
    }
    c_layout_all->Modified();

    TPad *pad_palette_all = new TPad("pad_palette_all", "", 0.94, 0.1, 1, 0.9);
    pad_palette_all->SetRightMargin(0.5);  // spazio per la palette
    pad_palette_all->Draw();
    pad_palette_all->cd();

    // copia histo
    TH2F *h_ref_all = (TH2F*) hist_purity_xQ2_zPt[2]->Clone("h_ref_all_for_palette");
    h_ref_all->SetStats(0);
    h_ref_all->Draw("COLZ");
    gPad->Update();

    // palette
    TPaletteAxis *pal_all = (TPaletteAxis*) h_ref_all->GetListOfFunctions()->FindObject("palette");
    if (pal_all) {
        pal_all->SetLabelSize(0.22);
        pal_all->SetTitleOffset(0.0);
        //pal->SetTitle("Purity");
        pal_all->SetX1NDC(0.05);
        pal_all->SetX2NDC(0.5);
        pal_all->SetY1NDC(0.05);
        pal_all->SetY2NDC(0.95);
    }


    pad_palette_all->Modified();
    c_layout_all->Update();

    c_layout_all->cd();
    // Scritta in alto a sinistra
    TLatex *globalTitle_all = new TLatex(0.12, 0.95, "4D Purity Distribution");
    globalTitle_all->SetNDC(true);         // coordinate normalizzate (relative alla canvas)
    //globalTitle->SetTextFont(62);      // font bold
    globalTitle_all->SetTextSize(0.04);   // dimensione
    globalTitle_all->SetTextAlign(13);     // 13 = left-top aligned
    globalTitle_all->Draw();
    TLatex *subTitle_all = new TLatex(0.08, 0.9, "(MC #pi+ -> RECO #pi+) / (MC all -> RECO #pi+)");
    subTitle_all->SetNDC(true);
    subTitle_all->SetTextSize(0.025);
    subTitle_all->SetTextAlign(13);
    subTitle_all->Draw();
    TLatex *subsub_all = new TLatex(0.17, 0.86, "(x_{B}, Q^{2}, z, P_{hT})");
    subsub_all->SetNDC(true);
    subsub_all->SetTextSize(0.03);
    subsub_all->SetTextAlign(13);
    subsub_all->Draw();

    // Fine
    c_layout_all->Update();
    c_layout_all->Write();










    //outFile.Write();
    outFile.Close();
    //chain.Close();

    cout << "ROOT output file: " << outputFile << endl;
}
