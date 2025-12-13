//-------------- Plotting relevant graphs
//--- Authors: Lorenzo Polizzi (lorenzo.polizzi@unife.i), Sara Pucillo (sara.pucillo@cern.ch), Nicolò Valle (nicolo.valle@cern.ch)
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

//------------- Generation of logarithmically spaced bin edges between xmin and xmax
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

//------------- Adjusts the position of the statistics box on a 2D histogram (top-right)
void SetStatsBox(TH2* hist) {
    TPaveStats* stats0 = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
    if (stats0) {
        stats0->SetX1NDC(0.75);  //left
        stats0->SetX2NDC(0.89);  //right
        stats0->SetY1NDC(0.68);  //bottom
        stats0->SetY2NDC(0.88);  //top
    }
}
void SetStatsBox2(TH2* hist) {
    TPaveStats* stats0 = (TPaveStats*)hist->GetListOfFunctions()->FindObject("stats");
    if (stats0) {
        stats0->SetX1NDC(0.12);  //left
        stats0->SetX2NDC(0.26);  //right
        stats0->SetY1NDC(0.68);  //bottom
        stats0->SetY2NDC(0.88);  //top
    }
}

//------------- Helper: binning functions in xB–Q2 and z–P_hT
int getBinIndex_xQ2(double xB, double Q2){
    // bin_xB
    double bin_xB[][2] = {{1e-4, 1e-3}, {1e-3, 3e-3}, {3e-3, 1e-2}, {1e-2, 4e-2}, {4e-2, 1}}; // if 10x100
    //double bin_xB[][2] = {{5e-5, 3e-4}, {3e-4, 1e-3}, {1e-3, 3e-3}, {3e-3, 1e-2}, {1e-2, 1}}; // if 18x275
    // bin Q2
    vector<vector<array<double,2>>> binning_Q2_for_xB = {
        {{1,2}, {2,100}},
        {{1,2}, {2,5}, {5,100}},
        {{1,2}, {2,5}, {5,20}, {20,100}},
        {{1,2}, {2,5}, {5,20}, {20,1000}},
        {{1,5}, {5,20}, {20,1000}}
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
            return -1;
        } else {
            binIndex += binning_Q2_for_xB[ix].size();
        }
    }

    return -1;
}

int getBinIndex_zPt(double z, double Pt){
    // bin_z
    double bin_z[][2] = {{0, 0.1}, {0.1, 0.2}, {0.2, 0.3}, {0.3, 0.5}, {0.5, 0.7}, {0.7, 1}};
    // bin_Pt
    vector<vector<array<double,2>>> binning_Pt_for_z = {
        {{0, 0.25}, {0.25, 0.5}, {0.5, 0.75}, {0.75, 1.25}, {1.25, 2.5}},
        {{0, 0.25}, {0.25, 0.5}, {0.5, 0.75}, {0.75, 1.25}, {1.25, 2.5}},
        {{0, 0.25}, {0.25, 0.5}, {0.5, 0.75}, {0.75, 1.25}, {1.25, 2.5}},
        {{0, 0.25}, {0.25, 0.5}, {0.5, 0.75}, {0.75, 1.25}, {1.25, 2.5}},
        {{0, 0.25}, {0.25, 0.5}, {0.5, 0.75}, {0.75, 1.25}, {1.25, 2.5}},
        {{0, 0.25}, {0.25, 0.5}, {0.5, 0.75}, {0.75, 1.25}, {1.25, 2.5}}
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
    double bin_z[][2] = {{0, 0.1}, {0.1, 0.2}, {0.2, 0.3}, {0.3, 0.5}, {0.5, 0.7}, {0.7, 1}};
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
    double bin_Pt[][2] = {{0, 0.25}, {0.25, 0.5}, {0.5, 0.75}, {0.75, 1.25}, {1.25, 2.5}};
    int binIndex = 1;
    for (int ip = 0; ip < 5; ip++){
        if (Pt >= bin_Pt[ip][0] && Pt <= bin_Pt[ip][1]){
            return ip+1;
        }
    }
    return -1;
}


//------------- Main function
void relevant_plots(int target_pdg = 211, const char* inputDir = "25.10_10x100") {

    //---set_ePIC_style();
    gROOT->ProcessLine("set_ePIC_style()");
    //gROOT->ForceStyle();

    // --- Hadron selection (PDG code)
    TString tag, label;
    switch (target_pdg) {
        case  211: tag = "pipos"; label = "#pi^{+}"; break;
        case -211: tag = "pineg"; label = "#pi^{-}"; break;
        case  321: tag = "kaonpos";  label = "K^{+}";    break;
        case -321: tag = "kaonneg";  label = "K^{-}";    break;
        default:
            tag   = Form("pdg%d", target_pdg);
            label = Form("PDG %d", target_pdg);
            break;
    }

    // Variables
    // --- electron
    double electron_px, electron_py, electron_pz, electron_mom, electron_Theta, electron_Phi, electron_ThetaDeg, electron_E, electron_W, electron_Q2, electron_ass;
    double electron_eta, electron_y;
    // -- hadron
    double hadron_mom, hadron_Q2, hadron_xB, hadron_xF, hadron_z, hadron_PhT, hadron_Phi_h, hadron_Phi_s, hadron_Phi_lab, hadron_Theta, hadron_eta, hadron_y, hadron_W, hadron_Mx;
    double helicity, eps, hadron_px, hadron_py, hadron_pz, el_px, el_py, el_pz, el_theta, el_phi, el_eta, el_mom, pr_mom, pr_px, pr_py, pr_pz, pr_phi, pr_theta, pr_eta;
    double rec_pdg, good_PID, el_rec_pdg, rec_pdg_mc;
    double el_ass_rec_pdg, el_ass_px, el_ass_py, el_ass_pz, el_ass_theta, el_ass_phi, el_ass_eta, el_ass_mom;
    double hadron_E;
    //
    double hadron_mom_mc, hadron_Q2_mc, hadron_xB_mc, hadron_xF_mc, hadron_z_mc, hadron_PhT_mc;
    double hadron_Phi_h_mc, hadron_Phi_s_mc, hadron_Phi_lab_mc, hadron_Theta_mc, hadron_eta_mc, hadron_y_mc, hadron_W_mc, hadron_Mx_mc;
    double hel_mc, eps_mc, hadron_px_mc, hadron_py_mc, hadron_pz_mc;
    double hadron_mc_index, index_mc;
    //
    double hadron_mom_all, hadron_Q2_all, hadron_xB_all, hadron_xF_all, hadron_z_all, hadron_PhT_all;
    double hadron_Phi_h_all, hadron_Phi_s_all, hadron_Phi_lab_all, hadron_Theta_all, hadron_eta_all, hadron_y_all, hadron_W_all, hadron_Mx_all;
    double hel_all, eps_all, hadron_px_all, hadron_py_all, hadron_pz_all;
    double good_PID_all, pdg_all;
    double hadron_index, index_all;

    //--- input file
    string inputDirStr = inputDir;
    TTree treeHadron(Form("%s_RECO", tag.Data()), Form("RECO %s", label.Data()));
    TTree treeHadron_MC(Form("%s_MC", tag.Data()), Form("MC %s", label.Data()));
    TChain chainElectron("Electron");
    TChain chainHadron_Reco("Hadron Reco");
    TChain chainHadron_MC("Hadron MC");

    int fileCount = 0;

    fs::path inPath(inputDirStr);
    if (fs::exists(inPath) && fs::is_regular_file(inPath) && inPath.extension() == ".root") {
        string filePath = inPath.string();
        chainElectron.Add(Form("%s/ElectronTree_MC", filePath.c_str()));
        chainHadron_MC.Add(Form("%s/HadronTree_MC", filePath.c_str()));
        chainHadron_Reco.Add(Form("%s/HadronTree_RECO", filePath.c_str()));
        fileCount = 1;
    } else {
        for (const auto &entry : fs::directory_iterator(inputDirStr)) {
            if (entry.path().extension() == ".root") {
                string filePath = entry.path().string();
                chainElectron.Add(Form("%s/ElectronTree_MC", filePath.c_str()));
                chainHadron_MC.Add(Form("%s/HadronTree_MC", filePath.c_str()));
                chainHadron_Reco.Add(Form("%s/HadronTree_RECO", filePath.c_str()));
                fileCount++;
            }
        }
    }

    if (fileCount == 0) {
        cerr << "No .root file found in " << inputDirStr << endl;
        return;
    }

    //--- output file
    TString outputFile = Form("%s/Relevant_plots_2510_epic_%s.root", inputDirStr.c_str(), tag.Data());
    TFile outFile(outputFile, "RECREATE");

    //--- Here we collect all the variables from the ttree
    if (chainElectron.GetNtrees() > 0) {
        chainElectron.SetBranchAddress("el_px_mc", &electron_px);
        chainElectron.SetBranchAddress("el_py_mc", &electron_py);
        chainElectron.SetBranchAddress("el_pz_mc", &electron_pz);
        chainElectron.SetBranchAddress("el_mom_mc", &electron_mom);
        chainElectron.SetBranchAddress("el_theta_mc", &electron_Theta);
        chainElectron.SetBranchAddress("el_phi_mc", &electron_Phi);
        chainElectron.SetBranchAddress("el_eta_mc", &electron_eta);
        chainElectron.SetBranchAddress("el_y_mc", &electron_y);
    }

    //Hadron reco
    chainHadron_Reco.SetBranchAddress("hadron_index", &hadron_index);
    chainHadron_Reco.SetBranchAddress("hadron_pdg", &rec_pdg);
    chainHadron_Reco.SetBranchAddress("hadron_pdg_mc", &rec_pdg_mc);
    chainHadron_Reco.SetBranchAddress("hadron_good_PID", &good_PID);
    chainHadron_Reco.SetBranchAddress("hadron_px", &hadron_px);
    chainHadron_Reco.SetBranchAddress("hadron_py", &hadron_py);
    chainHadron_Reco.SetBranchAddress("hadron_pz", &hadron_pz);
    chainHadron_Reco.SetBranchAddress("hadron_mom", &hadron_mom);
    //chainHadron_Reco.SetBranchAddress("W", &hadron_W);
    chainHadron_Reco.SetBranchAddress("hadron_Q2", &hadron_Q2);
    //chainHadron_Reco.SetBranchAddress("hadron_xF", &hadron_xF);
    chainHadron_Reco.SetBranchAddress("hadron_xB", &hadron_xB);
    chainHadron_Reco.SetBranchAddress("hadron_y", &hadron_y);
    chainHadron_Reco.SetBranchAddress("hadron_z", &hadron_z);
    chainHadron_Reco.SetBranchAddress("hadron_PhT", &hadron_PhT);
    chainHadron_Reco.SetBranchAddress("hadron_Phi_lab", &hadron_Phi_lab);
    chainHadron_Reco.SetBranchAddress("hadron_Theta", &hadron_Theta);
    chainHadron_Reco.SetBranchAddress("hadron_eta", &hadron_eta);
    chainHadron_Reco.SetBranchAddress("hadron_Phi_h", &hadron_Phi_h);
    //chainHadron_Reco.SetBranchAddress("Phi_s", &hadron_Phi_s);
    //chainHadron_Reco.SetBranchAddress("helicity", &helicity);
    //chainHadron_Reco.SetBranchAddress("hadron_Mx", &hadron_Mx);

    //Hadron MC
    chainHadron_MC.SetBranchAddress("hadron_index_mc", &index_mc);
    // Optional: if present, use MC PDG to select the hadron species
    double mc_pdg = 0;
    //bool has_mc_pdg = (chainHadron_MC.SetBranchAddress("hadron_pdg_mc", &mc_pdg) == 0);
    //if (!has_mc_pdg) has_mc_pdg = (chainHadron_MC.SetBranchAddress("hadron_pdg", &mc_pdg) == 0);
    chainHadron_MC.SetBranchAddress("hadron_pdg_mc", &mc_pdg);
    chainHadron_MC.SetBranchAddress("hadron_mom_mc", &hadron_mom_mc);
    chainHadron_MC.SetBranchAddress("hadron_Q2_mc", &hadron_Q2_mc);
    chainHadron_MC.SetBranchAddress("hadron_xB_mc", &hadron_xB_mc);
    //chainHadron_MC.SetBranchAddress("hadron_xF_mc", &hadron_xF_mc);
    chainHadron_MC.SetBranchAddress("hadron_z_mc", &hadron_z_mc);
    chainHadron_MC.SetBranchAddress("hadron_PhT_mc", &hadron_PhT_mc);
    chainHadron_MC.SetBranchAddress("hadron_Phi_lab_mc", &hadron_Phi_lab_mc);
    chainHadron_MC.SetBranchAddress("hadron_Phi_h_mc", &hadron_Phi_h_mc);
    //chainHadron_MC.SetBranchAddress("Phi_s_mc", &hadron_Phi_s_mc);
    chainHadron_MC.SetBranchAddress("hadron_Theta_mc", &hadron_Theta_mc);
    chainHadron_MC.SetBranchAddress("hadron_eta_mc", &hadron_eta_mc);
    chainHadron_MC.SetBranchAddress("hadron_y_mc", &hadron_y_mc);
    //chainHadron_MC.SetBranchAddress("hadron_W_mc", &hadron_W_mc);
    //chainHadron_MC.SetBranchAddress("hadron_Mx_mc", &hadron_Mx_mc);
    //chainHadron_MC.SetBranchAddress("helicity_mc", &hel_mc);
    //chainHadron_MC.SetBranchAddress("hadron_epsilon_mc", &eps_mc);
    chainHadron_MC.SetBranchAddress("hadron_px_mc", &hadron_px_mc);
    chainHadron_MC.SetBranchAddress("hadron_py_mc", &hadron_py_mc);
    chainHadron_MC.SetBranchAddress("hadron_pz_mc", &hadron_pz_mc);

    // Interesting hadron --- reco info
    treeHadron.Branch("mc_index", &hadron_index, "mc_index/I");
    treeHadron.Branch("rec_pdg", &rec_pdg, "rec_pdg/D");
    treeHadron.Branch("good_PID", &good_PID, "good_PID/D");
    treeHadron.Branch("px", &hadron_px, "hadron_px/D");
    treeHadron.Branch("py", &hadron_py, "hadron_py/D");
    treeHadron.Branch("pz", &hadron_pz, "hadron_pz/D");
    treeHadron.Branch("E", &hadron_E, "E/D");
    treeHadron.Branch("Mom", &hadron_mom, "Mom/D");
    treeHadron.Branch("Q2", &hadron_Q2, "Q2/D");
    treeHadron.Branch("xB", &hadron_xB, "xB/D");
    treeHadron.Branch("xF", &hadron_xF, "xF/D");
    treeHadron.Branch("z", &hadron_z, "z/D");
    treeHadron.Branch("PhT", &hadron_PhT, "PhT/D");
    treeHadron.Branch("Phi_lab", &hadron_Phi_lab, "Phi_h/D");
    treeHadron.Branch("Phi_h", &hadron_Phi_h, "Phi_h/D");
    treeHadron.Branch("Phi_s", &hadron_Phi_s, "Phi_s/D");
    treeHadron.Branch("theta", &hadron_Theta, "theta/D");
    treeHadron.Branch("eta", &hadron_eta, "eta/D");
    treeHadron.Branch("y", &hadron_y, "y/D");
    treeHadron.Branch("W", &hadron_W, "W/D");
    treeHadron.Branch("Mx", &hadron_Mx, "Mx/D");
    treeHadron.Branch("helicity", &helicity, "hel/D");
    treeHadron.Branch("epsilon", &eps, "eps/D");

    // Interesting hadron --- MC info
    treeHadron_MC.Branch("index", &index_mc, "index/I");
    treeHadron_MC.Branch("Mom_mc", &hadron_mom_mc, "Mom_mc/D");
    treeHadron_MC.Branch("Q2_mc", &hadron_Q2_mc, "Q2_mc/D");
    treeHadron_MC.Branch("xB_mc", &hadron_xB_mc, "xB_mc/D");
    treeHadron_MC.Branch("xF_mc", &hadron_xF_mc, "xF_mc/D");
    treeHadron_MC.Branch("z_mc", &hadron_z_mc, "z_mc/D");
    treeHadron_MC.Branch("PhT_mc", &hadron_PhT_mc, "PhT_mc/D");
    treeHadron_MC.Branch("Phi_lab_mc", &hadron_Phi_lab_mc, "Phi_lab_mc/D");
    treeHadron_MC.Branch("Phi_h_mc", &hadron_Phi_h_mc, "Phi_h_mc/D");
    treeHadron_MC.Branch("Phi_s_mc", &hadron_Phi_s_mc, "Phi_s_mc/D");
    treeHadron_MC.Branch("theta_mc", &hadron_Theta_mc, "theta_mc/D");
    treeHadron_MC.Branch("eta_mc", &hadron_eta_mc, "eta_mc/D");
    treeHadron_MC.Branch("y_mc", &hadron_y_mc, "y_mc/D");
    treeHadron_MC.Branch("W_mc", &hadron_W_mc, "W_mc/D");
    treeHadron_MC.Branch("Mx_mc", &hadron_Mx_mc, "Mx_mc/D");
    treeHadron_MC.Branch("helicity_mc", &hel_mc, "helicity_mc/D");
    treeHadron_MC.Branch("epsilon_mc", &eps_mc, "epsilon_mc/D");
    treeHadron_MC.Branch("pion_px_mc", &hadron_px_mc, "pion_px_mc/D");
    treeHadron_MC.Branch("pion_py_mc", &hadron_py_mc, "pion_py_mc/D");
    treeHadron_MC.Branch("pion_pz_mc", &hadron_pz_mc, "pion_pz_mc/D");

    /*
    treeHadron_all.Branch("index", &index_all, "index/I");
    treeHadron_all.Branch("Mom_all", &hadron_mom_all, "Mom_all/D");
    treeHadron_all.Branch("Q2_all", &hadron_Q2_all, "Q2_all/D");
    treeHadron_all.Branch("xB_all", &hadron_xB_all, "xB_all/D");
    treeHadron_all.Branch("xF_all", &hadron_xF_all, "xF_all/D");
    treeHadron_all.Branch("z_all", &hadron_z_all, "z_all/D");
    treeHadron_all.Branch("PhT_all", &hadron_PhT_all, "PhT_all/D");
    treeHadron_all.Branch("Phi_lab_all", &hadron_Phi_lab_all, "Phi_lab_all/D");
    treeHadron_all.Branch("Phi_h_all", &hadron_Phi_h_all, "Phi_h_all/D");
    treeHadron_all.Branch("Phi_s_all", &hadron_Phi_s_all, "Phi_s_all/D");
    treeHadron_all.Branch("theta_all", &hadron_Theta_all, "theta_all/D");
    treeHadron_all.Branch("eta_all", &hadron_eta_all, "eta_all/D");
    treeHadron_all.Branch("y_all", &hadron_y_all, "y_all/D");
    treeHadron_all.Branch("W_all", &hadron_W_all, "W_all/D");
    treeHadron_all.Branch("Mx_all", &hadron_Mx_all, "Mx_all/D");
    treeHadron_all.Branch("helicity_all", &hel_all, "helicity_all/D");
    treeHadron_all.Branch("epsilon_all", &eps_all, "epsilon_all/D");
    treeHadron_all.Branch("pion_px_all", &hadron_px_all, "pion_px_all/D");
    treeHadron_all.Branch("pion_py_all", &hadron_py_all, "pion_py_all/D");
    treeHadron_all.Branch("pion_pz_all", &hadron_pz_all, "pion_pz_all/D");
    */

    //
    // --- Distributions
    double bin = 200;
    double chi_min = 0.1;
    double chi_max = 2000;
    const double xmin_xB = 5e-5, xmax_xB = 1;
    const double xmin_Q2 = 1, xmax_Q2 = 2000.;
    auto make_bins = [](int bins, double min, double max){
        return CreateLogBinning(bins, min, max);
    };
    const auto log_chi2 = make_bins(bin, chi_min, chi_max);
    const auto log_bins_Q2 = make_bins(bin, xmin_Q2, xmax_Q2);
    const auto log_bins_xB = make_bins(bin, xmin_xB, xmax_xB);

    TH2D had_Q2VsXb_MC (Form("%s_Q2VsXb_MC", tag.Data()), Form("Correlation Q^{2} vs x_{B}  |  MC %s ; x_{B}; Q^{2} [GeV^{2}]",label.Data()), bin, log_bins_xB.data(), bin, log_bins_Q2.data());
    TH2D had_PhTvsZ_MC (Form("%s_PhTvsZ_MC", tag.Data()), Form("Correlation P_{hT} vs Z  |  MC %s ; z; P_{hT} [GeV]",label.Data()), bin, 0, 1, bin, 0, 3);

    TH1D had_evnt_chi2 (Form("%s_evnt_chi2", tag.Data()), Form("#chi^{2} EventBuilder PID | %s ; only EventBuilder | 1.2 < Mom < 8 GeV ; #chi^{2}; count",label.Data()), bin, -8, 8);
    TH1D had_m (Form("%s_best_mass", tag.Data()), Form("m extracted from #beta | %s ; 1.2 < Mom < 8 GeV ; m [GeV]; count",label.Data()), bin, 0, 1);
    TH1D had_deltaB (Form("%s_delta_beta", tag.Data()), Form("#beta_{meas} - #beta_{th} | %s ; 1.2 < Mom < 8 GeV ; #Delta_{#beta}; count",label.Data()), bin, -0.05, 0.05);
    TH1D had_Mx (Form("%s_missing_mass", tag.Data()), Form("missing mass | %s ; 1.2 < Mom < 8 GeV ; M_{x} [GeV]; count",label.Data()), bin, 0, 60);

    //--- Mom
    TH2D had_MomVsPhT (Form("%s_MomVsPhT", tag.Data()), Form("Correlation Mom vs P_{hT}  |  %s ; P_{hT} [GeV]; Mom [GeV]",label.Data()), bin, 0, 3, bin, 0, 20);
    TH2D had_MomVsXb (Form("%s_MomVsXb", tag.Data()), Form("Correlation Mom vs x_{B}  |  %s ; x_{B}; Mom [GeV]",label.Data()), bin, log_bins_xB.data(), bin, 0, 20);
    TH2D had_MomVsXf (Form("%s_MomVsXf", tag.Data()), Form("Correlation Mom vs x_{F}  |  %s ; x_{F}; Mom [GeV]",label.Data()), bin, -0.5, 0.5, bin, 0, 20);
    TH2D had_MomVsZ (Form("%s_MomVsZ", tag.Data()), Form("Correlation Mom vs Z  |  %s ; z; Mom [GeV]",label.Data()), bin, 0, 1, bin, 0, 20);
    TH2D had_MomVsY (Form("%s_MomVsY", tag.Data()), Form("Correlation Mom vs Y  |  %s ; y; Mom [GeV]",label.Data()), bin, 0.0, 1.0, bin, 0, 20);
    TH2D had_MomVsEta (Form("%s_MomVsEta", tag.Data()), Form("Correlation Mom vs Eta  |  %s ; Eta; Mom [GeV]",label.Data()), bin, -3, 3.5, bin, 0, 20);
    TH2D had_MomVsTheta (Form("%s_MomVsTheta", tag.Data()), Form("Correlation Mom vs Theta  |  %s ; Mom [GeV]; #theta [Rad]",label.Data()), bin, 0, 20, bin, 0, 180);
    TH2D had_MomVsPhi_h (Form("%s_MomVsPhi_h", tag.Data()), Form("Correlation Mom vs #Phi_{h}  |  %s ; #Phi_{h} [Rad]; Mom [GeV]",label.Data()), bin, -TMath::Pi(), TMath::Pi(), bin, 0, 20);
    TH2D had_MomVsMx (Form("%s_MomVsMx", tag.Data()), Form("Correlation Mom vs M^{2}_{x} | %s |; Mom [GeV]; M_{x} [GeV]",label.Data()), bin, 0, 20, bin, 0, 62);
    //--- Q2
    TH2D had_Q2VsXb (Form("%s_Q2VsXb", tag.Data()), Form("Correlation Q^{2} vs x_{B}  |  %s ; x_{B}; Q^{2} [GeV^{2}]",label.Data()), bin, log_bins_xB.data(), bin, log_bins_Q2.data());
    TH2D had_Q2VsXf (Form("%s_Q2VsXf", tag.Data()), Form("Correlation Q^{2} vs x_{F}  |  %s ; x_{F}; Q^{2} [GeV^{2}]",label.Data()), bin, -0.5, 0.5, bin, log_bins_Q2.data());
    TH2D had_Q2VsMom (Form("%s_Q2VsMom", tag.Data()), Form("Correlation Q^{2} vs Mom  |  %s ; Mom [GeV]; Q^{2} [GeV^{2}]",label.Data()), bin, 0, 20, bin, log_bins_Q2.data());
    TH2D had_Q2VsPhT (Form("%s_Q2VsPhT", tag.Data()), Form("Correlation Q^{2} vs P_{hT}  |  %s ; P_{hT} [GeV]; Q^{2} [GeV^{2}]",label.Data()), bin, 0, 3, bin, log_bins_Q2.data());
    TH2D had_Q2VsZ (Form("%s_Q2VsZ", tag.Data()), Form("Correlation Q^{2} vs Z  |  %s ; z; Q^{2} [GeV^{2}]",label.Data()), bin, 0, 1, bin, log_bins_Q2.data());
    TH2D had_Q2VsY (Form("%s_Q2VsY", tag.Data()), Form("Correlation Q^{2} vs Y  |  %s ; y; Q^{2} [GeV^{2}]",label.Data()), bin, 0.0, 1.0, bin, log_bins_Q2.data());
    TH2D had_Q2VsEta (Form("%s_Q2VsEta", tag.Data()), Form("Correlation Q^{2} vs Eta  |  %s ; Eta; Q^{2} [GeV^{2}]",label.Data()), bin, -3, 3.5, bin, log_bins_Q2.data());
    TH2D had_Q2VsPhi_h (Form("%s_Q2VsPhi_h", tag.Data()), Form("Correlation Q^{2} vs #Phi_{h}  |  %s ; #Phi_{h} [Rad]; Q^{2} [GeV^{2}]",label.Data()), bin, -TMath::Pi(), TMath::Pi(), bin, log_bins_Q2.data());
    //--- PhT
    TH2D had_PhTvsZ (Form("%s_PhTvsZ", tag.Data()), Form("Correlation P_{hT} vs Z  |  %s ; z; P_{hT} [GeV]",label.Data()), bin, 0, 1, bin, 0, 3);
    TH2D had_PhTvsXb (Form("%s_PhTvsXb", tag.Data()), Form("Correlation P_{hT} vs x_{B}  |  %s ; x_{B}; P_{hT} [GeV]",label.Data()), bin, log_bins_xB.data(), bin, 0, 3);
    TH2D had_PhTvsEta (Form("%s_PhTvsEta", tag.Data()), Form("Correlation P_{hT} vs Eta  |  %s ; Eta; P_{hT} [GeV]",label.Data()), bin, -3, 3.5, bin, 0, 3);
    TH2D had_PhTvsPhi_h (Form("%s_PhTvsPhi_h", tag.Data()), Form("Correlation P_{hT} vs #Phi_{h}  |  %s ; #Phi_{h} [Rad]; P_{hT} [GeV]",label.Data()), bin, -TMath::Pi(), TMath::Pi(), bin, 0, 3);
    //--- z
    TH2D had_zVsXb (Form("%s_zVsXb", tag.Data()), Form("Correlation Z vs x_{B}  |  %s ; x_{B}; z",label.Data()), bin, log_bins_xB.data(), bin, 0, 1);
    TH2D had_zVsXf (Form("%s_zVsXf", tag.Data()), Form("Correlation Z vs x_{F}  |  %s ; x_{F}; z",label.Data()), bin, -0.5, 0.5, bin, 0, 1);
    TH2D had_zVsEta (Form("%s_zVsEta", tag.Data()), Form("Correlation Z vs Eta  |  %s ; Eta; z",label.Data()), bin, -3, 3.5, bin, 0, 1);
    TH2D had_zVsPhi_h (Form("%s_zVsPhi_h", tag.Data()), Form("Correlation Z vs #Phi_{h}  |  %s ; #Phi_{h} [Rad]; z",label.Data()), bin, -TMath::Pi(), TMath::Pi(), bin, 0, 1);

    TH2D had_xBvsY (Form("%s_xBvsY", tag.Data()), Form("Correlation y vs x_{B}  |  %s ; x_{B}; y",label.Data()), bin, log_bins_xB.data(), bin, 0.0, 1.0);
    //--- Angles
    TH2D had_ThetaVsPhi_h (Form("%s_ThetaVsPhi_h", tag.Data()), Form("Correlation Theta vs #Phi_{h}  |  %s ; #Phi_{h} [Rad]; #theta [Rad]",label.Data()), bin, -TMath::Pi(), TMath::Pi(), bin, 0, 180);
    TH2D had_ThetaVsPhi_Lab (Form("%s_ThetaVsPhi_Lab", tag.Data()), Form("Correlation #theta vs #Phi_{Lab} | %s; #Phi_{Lab} [Rad]; #theta [Rad]",label.Data()), bin, -TMath::Pi(), TMath::Pi(), bin, 0, 180);


    //--- Array 2D and 4D for efficiency and purity
    int nBin_xQ2 = 16;
    int nBin_zPt = 30;
    double bin_xB[][2] = {{1e-4, 1e-3}, {1e-3, 3e-3}, {3e-3, 1e-2}, {1e-2, 4e-2}, {4e-2, 1}}; // 10x100
    //double bin_xB[][2] = {{5e-5, 3e-4}, {3e-4, 1e-3}, {1e-3, 3e-3}, {3e-3, 1e-2}, {1e-2, 1}}; // 18x275
    // bin_Q2
    double bin_Q2_full[][2] = {{1, 2}, {2, 5}, {5, 20}, {20, 1000}};
    double bin_xB_plot[] = {1e-4, 1e-3, 3e-3, 1e-2, 4e-2, 1}; // 10x100
    //double bin_xB_plot[] = {5e-5, 3e-4, 1e-3, 3e-3, 1e-2, 1}; // 18x275
    double bin_Q2_plot[] = {1, 2, 5, 20, 1000};
    // bin_z
    double bin_z[][2] = {{0, 0.1}, {0.1, 0.2}, {0.2, 0.3}, {0.3, 0.5}, {0.5, 0.7}, {0.7, 1}};
    double bin_z_plot[] = {0, 0.1, 0.2, 0.3, 0.5, 0.7, 1};
    double bin_Pt_plot[] = {0, 0.25, 0.5, 0.75, 1.25, 2.5};
    // bin_Q2
    vector<vector<array<double,2>>> binning_Pt_for_z = {
        {{0, 0.25}, {0.25, 0.5}, {0.5, 0.75}, {0.75, 1.25}, {1.25, 2.5}},
        {{0, 0.25}, {0.25, 0.5}, {0.5, 0.75}, {0.75, 1.25}, {1.25, 2.5}},
        {{0, 0.25}, {0.25, 0.5}, {0.5, 0.75}, {0.75, 1.25}, {1.25, 2.5}},
        {{0, 0.25}, {0.25, 0.5}, {0.5, 0.75}, {0.75, 1.25}, {1.25, 2.5}},
        {{0, 0.25}, {0.25, 0.5}, {0.5, 0.75}, {0.75, 1.25}, {1.25, 2.5}},
        {{0, 0.25}, {0.25, 0.5}, {0.5, 0.75}, {0.75, 1.25}, {1.25, 2.5}}
    };
    vector<vector<array<double,2>>> binning_Q2_for_xB = {
        {{1,2}, {2,5}},
        {{1,2}, {2,5}, {5,20}},
        {{1,2}, {2,5}, {5,20}, {20,100}},
        {{1,2}, {2,5}, {5,20}, {20,1000}},
        {{1,5}, {5,20}, {20,1000}}
    };
    vector<vector<double>> had_efficiency_xQ2(nBin_xQ2);
    vector<vector<double>> had_efficiency_zPt(nBin_zPt);
    vector<vector<vector<double>>> had_efficiency_xQ2_zPt(nBin_xQ2, vector<vector<double>> (nBin_zPt));
    vector<vector<vector<double>>> had_efficiency_xQ2_zPt_2(nBin_xQ2, vector<vector<double>> (nBin_zPt));
    vector<vector<vector<double>>> had_efficiency_xQ2_zPt_mc(nBin_xQ2, vector<vector<double>> (nBin_zPt));
    vector<vector<vector<double>>> had_purity_xQ2_zPt_num(nBin_xQ2, vector<vector<double>> (nBin_zPt));
    vector<vector<vector<double>>> had_purity_xQ2_zPt_den(nBin_xQ2, vector<vector<double>> (nBin_zPt));

    // --- relative plots
    vector<TH2D*> hist_efficiency_xQ2_zPt(nBin_xQ2);
    vector<TH2D*> hist_purity_xQ2_zPt(nBin_xQ2);
    for (int ix = 0; ix < nBin_xQ2; ix++){
      hist_efficiency_xQ2_zPt[ix] = new TH2D(Form("hist_efficiency4D_xQ2_%d", ix+1), Form("%s efficiency as P_{hT} vs z for bin (%d, x_{B}-Q^{2}); z; P_{hT} [GeV]", label.Data(), ix+1), 6, bin_z_plot, 5, bin_Pt_plot);
      hist_purity_xQ2_zPt[ix] = new TH2D(Form("hist_purity4D_xQ2_%d", ix+1),Form("%s purity as P_{hT} vs z for bin (%d, x_{B}-Q^{2}); z; P_{hT} [GeV]", label.Data(), ix+1), 6, bin_z_plot, 5, bin_Pt_plot);
    }
    // empty plots for inlet
    TH2D* hist_efficiency_xQ2_zPt_inlet = new TH2D("hist_efficiency4D_inlet", "hist_efficiency4D_inlet; z; P_{hT} [GeV]", 6, bin_z_plot, 5, bin_Pt_plot);
    TH2D* hist_purity_xQ2_zPt_inlet = new TH2D("hist_purity4D_inlet", "hist_purity4D_inlet; z; P_{hT} [GeV]", 6, bin_z_plot, 5, bin_Pt_plot);
    hist_efficiency_xQ2_zPt_inlet->SetStats(0);
    hist_purity_xQ2_zPt_inlet->SetStats(0);

    // --- Filling of the graphs
    // Hadron of interest
    Long64_t nEntries_had = chainHadron_Reco.GetEntries();
    Long64_t nEntries_hadMC = chainHadron_MC.GetEntries();

    for (Long64_t i = 0; i < nEntries_hadMC; i++) {
        chainHadron_MC.GetEntry(i);
        if (i % 100000 == 0) cout << "MC entry: " << i << "/" << nEntries_hadMC << endl;
        if(hadron_y_mc <= 0.99 && hadron_y_mc >= 0.01 && hadron_Q2_mc >= 1 && hadron_z_mc < 1){
            if (mc_pdg != target_pdg) continue;
            double bin_xQ2 = getBinIndex_xQ2(hadron_xB_mc, hadron_Q2_mc);
            double bin_zPt = getBinIndex_zPt(hadron_z_mc, hadron_PhT_mc);
            if(bin_xQ2 >= 0){
                if(bin_zPt >= 0){
                    had_efficiency_xQ2_zPt_mc[bin_xQ2-1][bin_zPt-1].push_back(hadron_y_mc); // compenso il fatto che i bin partano da 1
                }
            }
            had_Q2VsXb_MC.Fill(hadron_xB_mc, hadron_Q2_mc);
            had_PhTvsZ_MC.Fill(hadron_z_mc, hadron_PhT_mc);
            treeHadron_MC.Fill();
        }
    }

    /*
    for (Long64_t i = 0; i < nEntries_hadAll; i++) {
        chainHadron_all.GetEntry(i);
        if(good_PID_all != 0) continue;
        if(hadron_y_all <= 0.99 && hadron_y_all >= 0.01){ //same request as in epic_studies
            double bin_xQ2 = getBinIndex_xQ2(hadron_xB_all, hadron_Q2_all);
            double bin_zPt = getBinIndex_zPt(hadron_z_all, hadron_PhT_all);
            //
            if(bin_xQ2 >= 0){
                if(bin_zPt >= 0){
                    had_purity_xQ2_zPt_den[bin_xQ2-1][bin_zPt-1].push_back(hadron_y_all);
                    if(pdg_all == 211) had_purity_xQ2_zPt_num[bin_xQ2-1][bin_zPt-1].push_back(hadron_y_all); // compenso il fatto che i bin partano da 1
                }
            }
            //treeHadron_all.Fill();
        }
    }
        */

    for (Long64_t i = 0; i < nEntries_had; i++) {
        chainHadron_Reco.GetEntry(i);
        if (i % 100000 == 0) cout << "RECO entry: " << i << "/" << nEntries_had << endl;
            if(rec_pdg == target_pdg && hadron_z < 1 && hadron_Q2 >= 1){
                if(hadron_y <= 0.99 && hadron_y >= 0.01){ //same request as in epic_studies
                    double bin_xQ2 = getBinIndex_xQ2(hadron_xB, hadron_Q2);
                    double bin_zPt = getBinIndex_zPt(hadron_z, hadron_PhT);
                    if(bin_xQ2 >= 0){
                        if(bin_zPt >= 0){ // here we have reconstructed pion, but we do not ask if they are also MC pion -> contamination
                            had_efficiency_xQ2_zPt_2[bin_xQ2-1][bin_zPt-1].push_back(hadron_y); // compenso il fatto che i bin partano da 1
                            had_purity_xQ2_zPt_den[bin_xQ2-1][bin_zPt-1].push_back(hadron_y_all);
                        }
                    }
                }
            }
            if(hadron_y <= 0.99 && hadron_y >= 0.01 && good_PID == 0 && rec_pdg == target_pdg && rec_pdg_mc == target_pdg && hadron_z < 1 && hadron_Q2 >= 1){ //same request as in epic_studies + goodPID + specific pdg
                double bin_xQ2 = getBinIndex_xQ2(hadron_xB, hadron_Q2);
                double bin_zPt = getBinIndex_zPt(hadron_z, hadron_PhT);
                double bin_z = getBinIndex_z(hadron_z);
                double bin_Pt = getBinIndex_Pt(hadron_PhT);
                if(bin_xQ2 >= 0){
                    if(bin_zPt >= 0){
                        had_efficiency_xQ2_zPt[bin_xQ2-1][bin_zPt-1].push_back(hadron_y);
                        had_purity_xQ2_zPt_num[bin_xQ2-1][bin_zPt-1].push_back(hadron_y_all);
                    }
                }
                // Mom
                had_MomVsPhT.Fill(hadron_PhT, hadron_mom);
                had_MomVsEta.Fill(hadron_eta, hadron_mom);
                had_MomVsMx.Fill(hadron_mom, hadron_Mx);
                had_MomVsPhi_h.Fill(hadron_Phi_h, hadron_mom);
                had_MomVsTheta.Fill(hadron_mom, hadron_Theta);
                had_MomVsXb.Fill(hadron_xB, hadron_mom);
                had_MomVsXf.Fill(hadron_xF, hadron_mom);
                had_MomVsY.Fill(hadron_y, hadron_mom);
                had_MomVsZ.Fill(hadron_z, hadron_mom);
                // Q2
                had_Q2VsEta.Fill(hadron_eta, hadron_Q2);
                had_Q2VsMom.Fill(hadron_mom, hadron_Q2);
                had_Q2VsPhi_h.Fill(hadron_Phi_h, hadron_Q2);
                had_Q2VsPhT.Fill(hadron_PhT, hadron_Q2);
                had_Q2VsXb.Fill(hadron_xB, hadron_Q2);
                had_Q2VsXf.Fill(hadron_xF, hadron_Q2);
                had_Q2VsY.Fill(hadron_y, hadron_Q2);
                had_Q2VsZ.Fill(hadron_z, hadron_Q2);
                // PhT
                had_PhTvsEta.Fill(hadron_eta, hadron_PhT);
                had_PhTvsXb.Fill(hadron_xB, hadron_PhT);
                had_PhTvsZ.Fill(hadron_z, hadron_PhT);
                had_PhTvsPhi_h.Fill(hadron_Phi_h, hadron_PhT);
                // z
                had_zVsEta.Fill(hadron_eta, hadron_z);
                had_zVsPhi_h.Fill(hadron_Phi_h, hadron_z);
                had_zVsXb.Fill(hadron_xB, hadron_z);
                had_zVsXf.Fill(hadron_xF, hadron_z);
                //
                had_xBvsY.Fill(hadron_xB, hadron_y);
                had_ThetaVsPhi_h.Fill(hadron_Phi_h, hadron_Theta);
                had_ThetaVsPhi_Lab.Fill(hadron_Phi_lab, hadron_Theta);
                treeHadron.Fill();
            }
        }

        treeHadron.Write();
        treeHadron_MC.Write();

        for (int ix = 0; ix < nBin_xQ2; ix++) {
            for (int iz = 0; iz < nBin_zPt; iz++) {
                double n = had_efficiency_xQ2_zPt[ix][iz].size();
                double n_mc = had_efficiency_xQ2_zPt_mc[ix][iz].size();
                double n2 = had_efficiency_xQ2_zPt_2[ix][iz].size();
                double n_all = had_purity_xQ2_zPt_num[ix][iz].size();
                double n_all_den = had_purity_xQ2_zPt_den[ix][iz].size();
                double eff = (n > 0 || n_mc > 0) ? n/n_mc : 0.0;
                double err = (n_mc > 0) ? sqrt(eff * (1.0 - eff) / n_mc) : 0.0;
                double purity = (n > 0 || n_all > 0) ? n/n_all_den : 0.0;
                //hist_efficiency_xQ2_zPt[ix]->SetBinError(..., err);
                if (eff >= 1) eff = 1.0;
                if (purity >= 1) purity = 1.0;
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
    }


    vector<TH2D*> hists_had = {
        &had_MomVsPhT, &had_MomVsXb, &had_MomVsXf, &had_MomVsZ, &had_MomVsY, &had_MomVsEta,
        &had_MomVsTheta, &had_MomVsPhi_h, &had_MomVsMx,
        &had_Q2VsXb, &had_Q2VsXb_MC, &had_Q2VsXf, &had_Q2VsMom, &had_Q2VsPhT, &had_Q2VsZ, &had_Q2VsY, &had_Q2VsEta, &had_Q2VsPhi_h,
        &had_PhTvsZ, &had_PhTvsZ_MC, &had_PhTvsXb, &had_PhTvsEta, &had_PhTvsPhi_h,
        &had_zVsXb, &had_zVsXf, &had_zVsEta, &had_zVsPhi_h,
        &had_xBvsY, &had_ThetaVsPhi_h, &had_ThetaVsPhi_Lab,
    };


    // list to set the statbox2 on the canvas (suffix-based IDs)
    set<string> id_box2 = { "_MomVsXf", "_MomVsZ", "_MomVsY", "_MomVsEta", "_Q2VsXb", "_Q2VsY", "_zVsEta", "_zVsPhi_h", "_MomVsMass_RICH" };
    set<string> id_box3 = { "_MomVsXb", "_Q2VsXb", "_Q2VsXb_MC", "_PhTvsXb", "_zVsXb", "_xBvsY" };
    set<string> id_box4 = { "_Q2VsXb", "_Q2VsXb_MC", "_Q2VsXf", "_Q2VsMom", "_Q2VsPhT", "_Q2VsZ", "_Q2VsY", "_Q2VsEta", "_Q2VsPhi_h" };

    for (size_t i = 0; i < hists_had.size(); ++i) {

        const std::string histName = hists_had[i]->GetName();

        // Convert "pip_Q2VsXb_MC" -> "_Q2VsXb_MC" (so it matches your suffix sets)
        std::string key = histName;
        const std::string prefix = std::string(tag.Data()) + "_";
        if (key.rfind(prefix, 0) == 0) {               // starts_with(prefix)
            key = key.substr(prefix.size() - 1);       // keep leading "_"
        }

        TCanvas *c = new TCanvas(histName.c_str(), histName.c_str(), 800, 600);
        c->SetLogz();

        hists_had[i]->SetStats(0);
        hists_had[i]->SetTitle("");
        hists_had[i]->Draw("");

        if (id_box3.count(key)) c->SetLogx();
        if (id_box4.count(key)) c->SetLogy();
        // if (id_box2.count(key)) SetStatsBox2(hists_had[i]);
        // else SetStatsBox(hists_had[i]);

        TLatex Text_ePIC;
        Text_ePIC.SetTextSize(0.05);
        Text_ePIC.SetTextFont(62);
        //Text_ePIC.DrawLatexNDC(.15,.88,"ePIC Performance");  // performance plot
        //Text_ePIC.DrawLatexNDC(.15,.88,"ePIC Internal");  // for internal use only
        Text_ePIC.DrawLatexNDC(.15,.88,"ePIC Preliminary");
        //Text_ePIC.DrawLatexNDC(.15,.88,"ePIC Work in Progress"); // work in progress to be shown outside
        //Text_ePIC.DrawLatexNDC(.15,.88,"ePIC"); // final published version

        TLatex Text_com;
        Text_com.SetTextAlign(13);
        if (inputDirStr == "25.10_10x100") Text_com.DrawLatexNDC(.15,.85,"e+p, #sqrt{s} = 63.2 GeV");
        if (inputDirStr == "25.10_18x275") Text_com.DrawLatexNDC(.15,.85,"e+p, #sqrt{s} = 141 GeV");

        TLatex Text_date;
        Text_date.SetTextSize(0.035);
        Text_date.SetTextFont(52);
        Text_date.DrawLatexNDC(.65,.96,"Simu campaign: 10/2025");
        //Text_com.DrawLatexNDC(.15,.8,"L_{proj} = 10 fb^{-1}");

        TImage *logo = TImage::Open("EPIC-logo_black_small.png");
        TPad *pad2 = new TPad(Form("pad2_%s", histName.c_str()),Form("pad2_%s", histName.c_str()), 0.8, 0.8, 0.93, 0.93);
        pad2->Draw();
        pad2->cd();
        logo->Draw();

        c->Write();
    }


    // -- bin grid (Q2 vs xB)
    TCanvas *c_bin_xQ2 = new TCanvas(Form("%s_Q2_vs_xB_Bin",tag.Data()), Form("%s Q^{2} vs x_{B} bin",tag.Data()), 800, 700);
    c_bin_xQ2->SetLogz(); c_bin_xQ2->SetLogx(); c_bin_xQ2->SetLogy();
    had_Q2VsXb.Draw("COLZ");

    vector<TPolyLine*> rectangles;
    vector<TText*> labels;
    int bin_index = 1;
    // Loop over the bins
    for (int ix = 0; ix < 5; ++ix) {
        for (size_t iq = 0; iq < binning_Q2_for_xB[ix].size(); ++iq) {
            double xB[5] = {bin_xB[ix][0], bin_xB[ix][1], bin_xB[ix][1], bin_xB[ix][0], bin_xB[ix][0]};
            double Q2[5] = {binning_Q2_for_xB[ix][iq][0], binning_Q2_for_xB[ix][iq][0], binning_Q2_for_xB[ix][iq][1], binning_Q2_for_xB[ix][iq][1], binning_Q2_for_xB[ix][iq][0]};

            TPolyLine *rect = new TPolyLine(5, xB, Q2);
            rect->SetLineColor(kBlack);
            rect->SetLineWidth(2);
            rect->Draw("same");
            rectangles.push_back(rect);

            double x_center = sqrt(bin_xB[ix][0] * bin_xB[ix][1]);
            double Q2_center = sqrt(binning_Q2_for_xB[ix][iq][0] * binning_Q2_for_xB[ix][iq][1]);

            TText *label = new TText(x_center, Q2_center, Form("%d", bin_index++));
            label->SetTextAlign(22);
            label->SetTextSize(0.03);
            label->SetTextColor(kRed+1);
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
    Text_com.SetTextAlign(13);
    if (inputDirStr == "25.10_10x100") Text_com.DrawLatexNDC(.15,.85,"e+p, #sqrt{s} = 63.2 GeV");
    if (inputDirStr == "25.10_18x275") Text_com.DrawLatexNDC(.15,.85,"e+p, #sqrt{s} = 141 GeV");
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
    c_bin_xQ2->Update();
    c_bin_xQ2->Write();


    // -- bin grid (PhT vs z)
    TCanvas *c_bin_zPt = new TCanvas(Form("%s Pt_vs_z_Bin",tag.Data()), Form("%s P_{hT} vs z bin",tag.Data()), 800, 700);
    c_bin_zPt->SetLogz();
    had_PhTvsZ.Draw("COLZ");

    vector<TPolyLine*> gridLines_zp;
    vector<TText*> labels_zp;

    int bin_index_zp = 1;

    // General loop
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

            // rectangle
            TPolyLine *rect_zp = new TPolyLine(5, z, Pt);
            rect_zp->SetLineWidth(2);
            rect_zp->SetLineColor(kBlack);
            rect_zp->Draw("same");
            gridLines_zp.push_back(rect_zp);

            // centre
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



    //--- Efficiency
    double global_min = 1e9;
    double global_max = -1e9;
    for (int ix = 0; ix < nBin_xQ2; ++ix) {
        double min_tmp = hist_efficiency_xQ2_zPt[ix]->GetMinimum();
        double max_tmp = hist_efficiency_xQ2_zPt[ix]->GetMaximum();
        if (min_tmp < global_min) global_min = min_tmp;
        if (max_tmp > global_max) global_max = max_tmp;
    }

    for (int ix = 0; ix < nBin_xQ2; ++ix) {
        hist_efficiency_xQ2_zPt[ix]->SetMinimum(0); // now we use 0 and 1
        hist_efficiency_xQ2_zPt[ix]->SetMaximum(1);
    }

    // exclude zmin = zmax
    if (global_min == global_max) global_max = global_min + 1e-6;


    //--- Efficiency plots
    for (int ixQ2 = 0; ixQ2 < nBin_xQ2; ++ixQ2) {
        TCanvas *c_bin_zPt = new TCanvas(Form("%s c_efficiency4D_xQ2_%d", tag.Data(), ixQ2+1), Form("%s Canvas %d", tag.Data(), ixQ2+1), 800, 800);
        //c_bin_zPt->SetLogz();
        hist_efficiency_xQ2_zPt[ixQ2]->SetStats(0);
        hist_efficiency_xQ2_zPt[ixQ2]->Draw("colz");

        vector<TPolyLine*> gridLines_zp;
        vector<TText*> labels_zp;

        int bin_index_zp = 1;

        // General loop
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

                // rectangle
                TPolyLine *rect_zp = new TPolyLine(5, z, Pt);
                rect_zp->SetLineWidth(2);
                rect_zp->SetLineColor(kBlack);
                rect_zp->Draw("same");
                gridLines_zp.push_back(rect_zp);

                // centre
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
        { 999,  5,  8, 12, 15},
        { 2,  4,  7, 11, 14},
        { 1,  3,  6, 10,  0}
    };


    //--- Efficiency plot for all bins
    TCanvas *c_layout = new TCanvas(Form("%s c_4D_efficiency", tag.Data()), Form("%s All efficiency bins", tag.Data()), 1800, 1200);
    c_layout->cd();

    // === OVERLAY PAD for log axes ===
    TPad *pad_axes = new TPad(Form("%s pad_axes",tag.Data()), Form("%s Global Axes", tag.Data()), 0, 0, 1, 1);
    //pad_axes->SetFillStyle(4000);
    //pad_axes->SetFrameFillStyle(0);
    pad_axes->SetFrameLineWidth(0);
    pad_axes->Draw();
    pad_axes->cd();

    pad_axes->SetTickx(0);
    pad_axes->SetTicky(0);

    // NDC coordinates
    double x_ndc_min = 0.05, x_ndc_max = 0.96;
    double y_ndc_min = 0.05, y_ndc_max = 0.96;

    // x and y axis-range (log)
    double x_min = 1e-5, x_max = 1;
    double y_min = 1, y_max = 1000;

    // draw main lines
    TLine *xAxisLine = new TLine(x_ndc_min, y_ndc_min, x_ndc_max, y_ndc_min);
    xAxisLine->SetLineWidth(2);
    xAxisLine->SetNDC(true);
    xAxisLine->Draw();

    TLine *yAxisLine = new TLine(x_ndc_min, y_ndc_min, x_ndc_min, y_ndc_max);
    yAxisLine->SetLineWidth(2);
    yAxisLine->SetNDC(true);
    yAxisLine->Draw();

    // x-axis
    for (int i = -4; i <= 0; ++i) {
        double val  = pow(10, i);
        double logv = log10(val);

        double norm;
        if (logv <= -3.0) {
            // 10^-4 → 10^-3  (molto largo)
            norm = 2.0 * (logv + 4.0);
        }
        else if (logv <= -2.0) {
            // 10^-3 → 10^-2
            norm = (2.0 + (logv + 3.0));
        }
        else if (logv <= -1.0) {
            // 10^-2 → 10^-1
            norm = (3.0 + 0.7 * (logv + 2.0));
        }
        else {
            // 10^-1 → 1  (molto compresso)
            norm = (3.7 + 0.3 * (logv + 1.0));
        }

        // normalizzazione totale (2 + 1 + 0.7 + 0.3 = 4)
        norm /= 4.0;

        double pos = x_ndc_min + norm * (x_ndc_max - x_ndc_min);
        double dx = 0;
        //if (i == -4) dx = +0.1; // only for 18x275
        //if (i == -3) dx = -0.095; // 18x275
        if (i == -3) dx = -0.27;
        //if (i == -2) dx = 0.03;
        if (i == -2) dx = - 0.145;
        //if (i == -1) dx = -0.03;
        if (i == -1) dx = -0.08;
        if (i == 0) dx = -0.02;
        pos += dx;

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



    // x-axis title
    TLatex *xlabel = new TLatex(0.90, y_ndc_min - 0.025, "x_{B}");
    xlabel->SetTextSize(0.025);
    xlabel->SetTextAlign(22);
    xlabel->SetNDC(true);
    xlabel->Draw();

    // y-axis title
    for (int i = 0; i <= 3; ++i) {
        double val = pow(10, i);
        double norm;
        if (val <= 10.0)
            norm = 2.0 * (log10(val)) / 4.0;
        else if (val <= 100.0)
            norm = (2.0 + (log10(val) - 1.0)) / 4.0;
        else
            norm = (3.0 + (log10(val) - 2.0)) / 4.0;

        double pos = y_ndc_min + norm * (y_ndc_max - y_ndc_min);
        // micro shift
        double dy = 0.0;
        if (val == 10.0)  dy = 0.1;   // higher
        if (val == 100.0) dy = 0.08;    
        pos += dy;

        TLine *tick = new TLine(x_ndc_min, pos, x_ndc_min - 0.01, pos);
        tick->SetNDC(true);
        tick->Draw();

        TLatex *lab = new TLatex(x_ndc_min - 0.015, pos, Form("10^{%d}", i));
        lab->SetTextAlign(32);
        lab->SetTextSize(0.02);
        lab->SetNDC(true);
        lab->Draw();
    }


    // y-axis title
    TLatex *ylabel = new TLatex(x_ndc_min - 0.025, 0.71, "Q^{2} (GeV^{2})");
    ylabel->SetTextSize(0.025);
    ylabel->SetTextAngle(90);
    ylabel->SetTextAlign(22);
    ylabel->SetNDC(true);
    ylabel->Draw();

    pad_axes->Modified();
    pad_axes->Update();

    c_layout->cd();

    // pad info
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
    // Loop
    for (int iRow = 0; iRow < nRows; ++iRow) {
        for (int iCol = 0; iCol < nCols; ++iCol) {
            int idx = layout[iRow][iCol]; // idx is the bin index
            if (idx == 0) continue;

	   
            // normalized coordinates: (x1,y1) bottom-left, (x2,y2) top-right
            double x1 = xMargin + iCol * padW;
            double x2 = x1 + padW;
            // ATTENZIONE: iRow = 0 significa riga in alto nella tua matrice, per cui dobbiamo invertire l'ordine per le y in coordinate normalized
            double y2 = 1.0 - yMargin - iRow * padH;
            double y1 = y2 - padH;

            // pìonly 1 pad
            TString padName = Form("%s_pad_r%d_c%d", tag.Data(), iRow, iCol);
	    if (idx == 999) padName = Form("%s_inlet", tag.Data());
	    
            c_layout->cd();

      	    
            TPad *pad = new TPad(padName, padName, x1, y1, x2, y2);
            pad->SetRightMargin(idx == 999 ? 0.05 : 0.0);
            pad->SetLeftMargin(idx == 999 ? 0.12 : 0.0);
            pad->SetBottomMargin(idx == 999 ? 0.12 : 0.0);
            pad->SetTopMargin(0.0);
            pad->Draw();
            pad->cd();

	

            // Histogram draw (no colorbar)
	    if (idx != 999){
	      hist_efficiency_xQ2_zPt[idx - 1]->SetTitle("");
	      hist_efficiency_xQ2_zPt[idx - 1]->Draw("col"); // idx-1 because 0-based vector
	    }
	    else{
	      hist_efficiency_xQ2_zPt_inlet->SetTitle("");
	      hist_efficiency_xQ2_zPt_inlet->Draw("");
	    }
	      for (auto &rect : grid_zp) rect->DrawClone("same");
	   
            // Come back to the canvas, before moving to the next pad
            c_layout->cd();
        }
    }
    c_layout->Modified();

    TPad *pad_palette = new TPad(Form("%s pad_palette",tag.Data()), "", 0.94, 0.1, 1, 0.9);
    pad_palette->SetRightMargin(0.5);
    pad_palette->Draw();
    pad_palette->cd();

    //
    TH2F *h_ref = (TH2F*) hist_efficiency_xQ2_zPt[2]->Clone(Form("%s h_ref_for_palette",tag.Data()));
    h_ref->SetStats(0);
    h_ref->Draw("COLZ");
    gPad->Update();

    //
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
    TLatex *globalTitle = new TLatex(0.10, 0.95, Form("%s 4D Efficiency Distribution", label.Data()));
    globalTitle->SetNDC(true);
    //globalTitle->SetTextFont(62);
    globalTitle->SetTextSize(0.04);
    globalTitle->SetTextAlign(13);
    globalTitle->Draw();
    TLatex *subTitle = new TLatex(0.15, 0.9, Form("RECO %s / MC %s",label.Data(), label.Data()));
    subTitle->SetNDC(true);
    subTitle->SetTextSize(0.03);
    subTitle->SetTextAlign(13);
    subTitle->Draw();
    TLatex *subsub = new TLatex(0.17, 0.86, "(x_{B}, Q^{2}, z, P_{hT})");
    subsub->SetNDC(true);
    subsub->SetTextSize(0.03);
    subsub->SetTextAlign(13);
    subsub->Draw();
    // End
    c_layout->Update();
    c_layout->Write();



    //--- purity
    double global_min_all = 1e9;
    double global_max_all = -1e9;
    for (int ix = 0; ix < nBin_xQ2; ++ix) {
        double min_tmp = hist_purity_xQ2_zPt[ix]->GetMinimum();
        double max_tmp = hist_purity_xQ2_zPt[ix]->GetMaximum();
        if (min_tmp < global_min_all) global_min_all = min_tmp;
        if (max_tmp > global_max_all) global_max_all = max_tmp;
        if(global_max_all > 1) global_max_all = 1;
    }

    for (int ix = 0; ix < nBin_xQ2; ++ix) {
        hist_purity_xQ2_zPt[ix]->SetMinimum(0);
        hist_purity_xQ2_zPt[ix]->SetMaximum(1);
    }

    // exclude zmin = zmax
    if (global_min_all == global_max_all) global_max_all = global_min_all + 1e-6;


    //--- purity plots
    for (int ixQ2 = 0; ixQ2 < nBin_xQ2; ++ixQ2) {
        TCanvas *c_bin_zPt = new TCanvas(Form("%s c_purity4D_xQ2_%d", tag.Data(), ixQ2+1), Form("%s Canvas %d", tag.Data(), ixQ2+1), 800, 800);
        //c_bin_zPt->SetLogz();
        hist_purity_xQ2_zPt[ixQ2]->SetStats(0);
        hist_purity_xQ2_zPt[ixQ2]->Draw("colz");

        vector<TPolyLine*> gridLines_zp;
        vector<TText*> labels_zp;

        int bin_index_zp = 1;

        // General loop
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

                // rectangle
                TPolyLine *rect_zp = new TPolyLine(5, z, Pt);
                rect_zp->SetLineWidth(2);
                rect_zp->SetLineColor(kBlack);
                rect_zp->Draw("same");
                gridLines_zp.push_back(rect_zp);

                // centre
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


    //--- Purity plot for all bins
    TCanvas *c_layout_all = new TCanvas(Form("%s c_4D_purity", tag.Data()), Form("%s All purity bins", tag.Data()), 1800, 1200);
    c_layout_all->cd();

    pad_axes->Draw();
    pad_axes->cd();

    pad_axes->SetTickx(0);
    pad_axes->SetTicky(0);

    xAxisLine->Draw();
    yAxisLine->Draw();

    

    xlabel->Draw();

    

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

    // Loop and creation of the pads
    for (int iRow = 0; iRow < nRows; ++iRow) {
        for (int iCol = 0; iCol < nCols; ++iCol) {
            int idx = layout[iRow][iCol]; // idx is the bin index
            if (idx == 0) continue;

            // normalized coordinates: (x1,y1) bottom-left, (x2,y2) top-right
            double x1 = xMargin + iCol * padW;
            double x2 = x1 + padW;
            // ATTENZIONE: iRow = 0 significa riga in alto nella tua matrice, per cui dobbiamo invertire l'ordine per le y in coordinate normalized
            double y2 = 1.0 - yMargin - iRow * padH;
            double y1 = y2 - padH;

            // only 1 pad
            TString padName = Form("%s_pad_r%d_c%d", tag.Data(), iRow, iCol);
	    if (idx == 999) padName = Form("%s_inlet", tag.Data());
            c_layout_all->cd();

            TPad *pad = new TPad(padName, padName, x1, y1, x2, y2);
            pad->SetRightMargin(idx == 999 ? 0.05 : 0.0);
            pad->SetLeftMargin(idx == 999 ? 0.12 : 0.0);
            pad->SetBottomMargin(idx == 999 ? 0.12 : 0.0);
            pad->SetTopMargin(0.0);
            pad->Draw();
            pad->cd();

            // Histogram draw (mo colorbar)
	    if (idx != 999){
	      hist_purity_xQ2_zPt[idx - 1]->SetTitle("");
	      hist_purity_xQ2_zPt[idx - 1]->Draw("col"); // idx-1 perché vettore 0-based
	    }
	    else{
	      hist_purity_xQ2_zPt_inlet->SetTitle("");
	      hist_purity_xQ2_zPt_inlet->Draw("");
	    }
	      
            for (auto &rect : grid_zp) rect->DrawClone("same");
            // Come back to the canvas, before moving to the next pad
            c_layout_all->cd();
        }
    }

    c_layout_all->Modified();

    TPad *pad_palette_all = new TPad(Form("%s_pad_palette_all",tag.Data()), "", 0.94, 0.1, 1, 0.9);
    pad_palette_all->SetRightMargin(0.5);  // spazio per la palette
    pad_palette_all->Draw();
    pad_palette_all->cd();

    //
    TH2F *h_ref_all = (TH2F*) hist_purity_xQ2_zPt[2]->Clone(Form("%s h_ref_all_for_palette",tag.Data()));
    h_ref_all->SetStats(0);
    h_ref_all->Draw("COLZ");
    gPad->Update();

    //
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
    TLatex *globalTitle_all = new TLatex(0.12, 0.95, Form("%s 4D Purity Distribution",label.Data()));
    globalTitle_all->SetNDC(true);
    //globalTitle->SetTextFont(62);
    globalTitle_all->SetTextSize(0.04);
    globalTitle_all->SetTextAlign(13);
    globalTitle_all->Draw();
    TLatex *subTitle_all = new TLatex(0.08, 0.9, Form("(MC %s -> RECO %s) / (MC all -> RECO %s)", label.Data(), label.Data(), label.Data()));
    subTitle_all->SetNDC(true);
    subTitle_all->SetTextSize(0.025);
    subTitle_all->SetTextAlign(13);
    subTitle_all->Draw();
    TLatex *subsub_all = new TLatex(0.17, 0.86, "(x_{B}, Q^{2}, z, P_{hT})");
    subsub_all->SetNDC(true);
    subsub_all->SetTextSize(0.03);
    subsub_all->SetTextAlign(13);
    subsub_all->Draw();

    c_layout_all->Update();
    c_layout_all->Write();


    //outFile.Write();
    outFile.Close();
    //chain.Close();

    cout << "ROOT output file: " << outputFile << endl;
}
