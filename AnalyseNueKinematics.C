// Program that loads in the ttree created from the NueKinematics Module and makes some plots as does the Larsoft Module. 
// WARNING: make sure a plots folder exits in the current directory otherwise the histogram will not get saved and get an error. 
// Execute over the NueKinematics.root file which has been run over the larsoft module NueKinematics_module.cc
// This file can be run by simply typing root AnalyseNueKinematics.C

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TExec.h"
#include "TStyle.h"
#include "TColor.h"

#include <iostream>
#include <fstream>

// Function to make a Th1D histogram (Histogram variable name, output file path/name.extension, options for histogram)
void Plot1DHist (TH1D *histogram, const char * print_name, const char *Options){
	auto *c1 = new TCanvas(); // Create a TCanvas
	c1->cd();
	
    histogram->Draw();       // Draw hist and set the options
    histogram->SetOption(Options);
    histogram->SetLineWidth(2);

	c1->Print(print_name);   // Save the histogram
    c1->Close();             // Close the Canvas
}

// Function to make a Th2D histogram (Histogram variable name, output file path/name, options for histogram)
void Plot2DHist (TH2D *histogram, const char * print_name, const char *Options){
	
    auto *c1 = new TCanvas(); // Create a TCanvas
	c1->cd();

   gStyle->SetPalette(kBird,0,1);
   histogram->SetLineWidth(1);
   histogram->SetLineColor(kBlack);


    histogram->Draw(); // Draw hist and set the options
    histogram->SetOption(Options);
    
	
    c1->Print(print_name); // Save the histogram
    c1->Close(); // Close the Canvas
}

// Function that plots 2 th2d histograms with differnt pallete colours on the same plot(Histogram variable name, output file path/name, options for histogram)
void Plot2DHistSAME (TH2D *histogram, TH2D *histogram2, const char * print_name){
	
    auto *c1 = new TCanvas(); // Create a TCanvas
	c1->cd();

    TPaveText *text = new TPaveText(5.,160., 7., 180.);
    TExec *ex1 = new TExec("ex1","gStyle->SetPalette(kBird, 0, 0.1);");         // Colour pallete for BNB
    TExec *ex2 = new TExec("ex2","gStyle->SetPalette(kCherry , 0, 0.1);"); // Colour pallete for NuMI

	histogram->Draw("COL"); // Draw hist and set the options
    ex1->Draw();

    ex2->Draw();
    histogram2->Draw("COL, SAME");

    text->AddText("BNB -- Blue, NuMI -- Red");
    text->Draw();

    c1->Print(print_name); // Save the histogram
    c1->Close(); // Close the Canvas

    gStyle->SetPalette(kBird); // Reset back to the default colour pallete. 
}

// Function to make a Th2D histogram (Histogram variable name, output file path/name, options for histogram)
void Plot3DHist (TH3D *histogram, const char * print_name, const char *Options){
	auto *c1 = new TCanvas(); // Create a TCanvas
	c1->cd();

	histogram->Draw(); // Draw hist and set the options
    histogram->SetOption(Options);
	
    c1->Print(print_name); // Save the histogram
    c1->Close(); // Close the Canvas
}


// Main Function
void AnalyseNueKinematics() {
    gStyle->SetOptStat(0);

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                  Initialize
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    TFile *MyFile = new TFile("plots/Plots_NeuKinematics.root","RECREATE");
    if ( MyFile->IsOpen() ) printf("File opened successfully\n");

    // Create Histograms
    // Nue All BNB
    TH3D*   hNue_E_vs_Theta_vs_Phi = new TH3D("Nue_E_vs_Theta_vs_Phi_All","Nue_E_vs_Theta_All; Energy [GeV]; Theta [degrees]; Phi [degrees]", 20., 0., 10. , 10., 0., 180, 10., -180., 180 );

    TH2D*   hNue_E_vs_Theta = new TH2D("Nue_E_vs_Theta_All","Nue_E_vs_Theta_All; Energy [GeV]; Theta [degrees]",15., 0., 8. , 10., 0., 180);
    TH2D*   hNue_E_vs_Phi   = new TH2D("Nue_E_vs_Phi_All","Nue_E_vs_Phi_All; Energy [GeV]; Phi [degrees]",10., 0., 10. , 10., -180., 180);

    TH1D* 	hNue_Energy = new TH1D("Nue_Energy_All","Nue_Energy_All; E [GeV]; Events",50., 0., 5);
    TH1D* 	hNue_Theta =  new TH1D("Nue_Theta_All","Nue_Theta_All; Theta [Degrees]; Events", 100., 0., 1);
    TH1D* 	hNue_Phi =    new TH1D("Nue_Phi_All","Nue_Phi_All; Phi [Degrees]; Events", 10., -180., 180);

    
    // Electron all BNB
    TH3D*   helectron_E_vs_Theta_vs_Phi = new TH3D("electron_E_vs_Theta_vs_Phi_All","electron_E_vs_Theta_All; Energy [GeV]; Theta [degrees]; Phi [degrees]", 20., 0., 10. , 10., 0., 180, 10., -180., 180 );

    TH2D*   helectron_E_vs_Theta = new TH2D("electron_E_vs_Theta_All","electron_E_vs_Theta_All; Energy [GeV]; Theta [degrees]",20., 0., 10. , 10., 0., 180);
    TH2D*   helectron_E_vs_Phi   = new TH2D("electron_E_vs_Phi_All","electron_E_vs_Phi_All; Energy [GeV]; Phi [degrees]",20., 0., 10. , 10., -180., 180);
  
    TH1D* 	helectron_Energy = new TH1D("electron_Energy_All","electron_Energy_All; E [GeV]; Events",50., 0., 5.);
    TH1D* 	helectron_Theta  = new TH1D("electron_Theta_All","electron_Theta_All; Theta [Degrees]; Events", 10., 0., 180);
    TH1D* 	helectron_Phi    = new TH1D("electron_Phi_All","electron_Phi_All; Phi [Degrees]; Events", 10., -180., 180);


    //-- NuMI Plots --
    // Nue NuMI
    TH2D*   hNue_E_vs_Theta_NuMI = new TH2D("Nue_E_vs_Theta_All_NuMI"," ",15., 0., 8. , 10., 0., 180);
    TH2D*   hNue_E_vs_Phi_NuMI   = new TH2D("Nue_E_vs_Phi_All_NuMI"," ",10., 0., 10. , 10., -180., 180);


    // El NuMI
    TH2D*   helectron_E_vs_Theta_NuMI = new TH2D("electron_E_vs_Theta_All_NuMI",  "  ",20., 0., 10. , 10., 0., 180);
    TH2D*   helectron_E_vs_Phi_NuMI   = new TH2D("electron_E_vs_Phi_All_NuMI","  ",20., 0., 10. , 10., -180., 180);
    // --  -- -- 


    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                  Read In TTree
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
	// Open the NueKinematics File
	TFile *fileIN = TFile::Open("NueKinematics.root");
	if (fileIN == 0) {
	  // if we cannot open the file, print an error message and return immediatly
	  printf("Error: cannot open NueKinematics.root!\n");
	  return;
	}

	// Create the tree reader and its data containers
	TTreeReader myReader("microboonewvtof/EventTree", fileIN);

	TTreeReaderValue<int>    PDGRV(myReader, "PDG"); 					// Load in the PDG variable
	TTreeReaderValue<double> NueEnergyRV(myReader, "NueEnergy"); 	// Load in the NueEnergy variable
    TTreeReaderValue<double> NueThetaRV(myReader, "NueTheta"); 	    // Load in the NueTheta variable
    TTreeReaderValue<double> NuePhiRV(myReader, "NuePhi"); 	        // Load in the NuePhi variable

    TTreeReaderValue<double> ElectronEnergyRV(myReader, "ElectronEnergy"); 	// Load in the ElectronEnergy variable
    TTreeReaderValue<double> ElectronThetaRV(myReader, "ElectronTheta"); 	// Load in the ElectronTheta variable
    TTreeReaderValue<double> ElectronPhiRV(myReader, "ElectronPhi"); 	    // Load in the ElectronPhi variable


    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                         Loop over all entries of the TTree or TChain.
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	while (myReader.Next()) {

        // Fill Nue Histograms
        hNue_E_vs_Theta_vs_Phi ->Fill(*NueEnergyRV, *NueThetaRV, *NuePhiRV );

        hNue_E_vs_Theta->Fill(*NueEnergyRV, *NueThetaRV );
        hNue_E_vs_Phi->Fill(*NueEnergyRV, *NuePhiRV);

		hNue_Energy->Fill(*NueEnergyRV); 
        hNue_Theta->Fill(*NueThetaRV); 
        hNue_Phi->Fill(*NuePhiRV); 
        
        // Fill Electron Histograms
        helectron_E_vs_Theta_vs_Phi ->Fill(*ElectronEnergyRV, *ElectronThetaRV, *ElectronPhiRV );

        helectron_E_vs_Theta->Fill(*ElectronEnergyRV, *ElectronThetaRV );
        helectron_E_vs_Phi->Fill(*ElectronEnergyRV, *ElectronPhiRV);

        helectron_Energy->Fill(*ElectronEnergyRV);
        helectron_Theta->Fill(*ElectronThetaRV);
        helectron_Phi->Fill(*ElectronPhiRV);

	}

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                      Read in the Ttree from Coltons Selection NuMI
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
	// Open the NueKinematics File
	TFile *fileIN2 = TFile::Open("KrishTree.root");
	if (fileIN2 == 0) {
	  // if we cannot open the file, print an error message and return immediatly
	  printf("Error: cannot open KrishTree.root!\n");
	  return;
	}

	// Create the tree reader and its data containers
	TTreeReader myReader2("KrishTree", fileIN2);

	TTreeReaderValue<double> NuMI_NueEnergy_RV(myReader2, "mc_nu_energy"); 	        // Load in the NueEnergy variable
    TTreeReaderValue<double> NuMI_NueTheta_RV(myReader2, "mc_cos_theta"); 	        // Load in the NueTheta variable
    TTreeReaderValue<double> NuMI_NuePhi_RV(myReader2, "mc_phi"); 	                // Load in the NuePhi variable

    TTreeReaderValue<double> NuMI_ElectronEnergy_RV(myReader2, "mc_ele_energy"); 	// Load in the ElectronEnergy variable
    TTreeReaderValue<double> NuMI_ElectronTheta_RV(myReader2, "mc_ele_theta"); 	    // Load in the ElectronTheta variable
    TTreeReaderValue<double> NuMI_ElectronPhi_RV(myReader2, "mc_ele_phi"); 	        // Load in the ElectronPhi variable


    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                         Loop over all entries of the TTree or TChain.
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
    double NuMI_theta{0};

	while (myReader2.Next()) {
        NuMI_theta = acos(*NuMI_NueTheta_RV) * (180 / 3.1415);

        // // Fill Nue Histograms
        // hNue_E_vs_Theta_vs_Phi ->Fill(*NueEnergyRV, *NueThetaRV, *NuePhiRV );

        hNue_E_vs_Theta_NuMI->Fill(*NuMI_NueEnergy_RV,  NuMI_theta );
        hNue_E_vs_Phi_NuMI->Fill(*NuMI_NueEnergy_RV, *NuMI_NuePhi_RV);

		// hNue_Energy->Fill(*NueEnergyRV); 
        // hNue_Theta->Fill(*NueThetaRV); 
        // hNue_Phi->Fill(*NuePhiRV); 
        
        // // Fill Electron Histograms
        // helectron_E_vs_Theta_vs_Phi ->Fill(*ElectronEnergyRV, *ElectronThetaRV, *ElectronPhiRV );

        helectron_E_vs_Theta_NuMI->Fill(*NuMI_ElectronEnergy_RV, *NuMI_ElectronTheta_RV );
        helectron_E_vs_Phi_NuMI->Fill(*NuMI_ElectronEnergy_RV, *NuMI_ElectronPhi_RV);

        // helectron_Energy->Fill(*ElectronEnergyRV);
        // helectron_Theta->Fill(*ElectronThetaRV);
        // helectron_Phi->Fill(*ElectronPhiRV);

	}






    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                  Make Histograms
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // Nue BNB Only
    Plot3DHist(hNue_E_vs_Theta_vs_Phi, "plots/BNB/Nue_Energy_vs_Theta_vs_Phi.png", "LEGO2");

    Plot2DHist(hNue_E_vs_Theta, "plots/BNB/Nue_Energy_vs_Theta.png", "COLZ,TEXT00" );
    Plot2DHist(hNue_E_vs_Phi,   "plots/BNB/Nue_Energy_vs_Phi.png",   "COLZ,TEXT00" );
    Plot1DHist(hNue_Energy,     "plots/BNB/Nue_Energy.png",          "HIST" );
    Plot1DHist(hNue_Theta,      "plots/BNB/Nue_Theta.png",           "HIST" ); 
    Plot1DHist(hNue_Phi,        "plots/BNB/Nue_Phi.png",             "HIST" );

    // NuMi Nue
    Plot2DHist(hNue_E_vs_Theta_NuMI, "plots/NuMI/Nue_Energy_vs_Theta_NuMI.png", "COLZ, TEXT00" );
    Plot2DHist(hNue_E_vs_Phi_NuMI, "plots/NuMI/Nue_Energy_vs_Phi_NuMI.png", "COLZ, TEXT00" );

    // Electron BNB Only
    Plot3DHist(helectron_E_vs_Theta_vs_Phi, "plots/BNB/El_Energy_vs_Theta_vs_Phi.png", "LEGO2" );

    Plot2DHist(helectron_E_vs_Theta,"plots/BNB/El_Energy_vs_Theta.png", "COLZ,TEXT00" );
    Plot2DHist(helectron_E_vs_Phi,  "plots/BNB/El_Energy_vs_Phi.png",   "COLZ,TEXT00" );
    Plot1DHist(helectron_Energy,    "plots/BNB/El_Energy.png",          "HIST" );
    Plot1DHist(helectron_Theta,     "plots/BNB/El_Theta.png",           "HIST" ); 
    Plot1DHist(helectron_Phi,       "plots/BNB/El_Phi.png",             "HIST" );

     // NuMi Electron
    Plot2DHist(helectron_E_vs_Theta_NuMI, "plots/NuMI/El_Energy_vs_Theta_NuMI.png", "COLZ, TEXT00" );
    Plot2DHist(helectron_E_vs_Phi_NuMI, "plots/NuMI/El_Energy_vs_Phi_NuMI.png", "COLZ, TEXT00" );

    // BNB and NuMI Plots for Nue (Superimposed)
    Plot2DHistSAME(hNue_E_vs_Theta , hNue_E_vs_Theta_NuMI, "plots/Nue_Energy_vs_Theta_Overlaid.png" );
    Plot2DHistSAME(hNue_E_vs_Phi , hNue_E_vs_Phi_NuMI, "plots/Nue_Energy_vs_Phi_Overlaid.png" );

    Plot2DHistSAME(helectron_E_vs_Theta , helectron_E_vs_Theta_NuMI, "plots/El_Energy_vs_Theta_Overlaid.png" );
    Plot2DHistSAME(helectron_E_vs_Phi , helectron_E_vs_Phi_NuMI, "plots/El_Energy_vs_Phi_Overlaid.png" );


    MyFile->Write(); // Save to a root file 
    MyFile->Close();
    gSystem->Exit(1); // Quit ROOT

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                     END
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

} // END MAIN FUNCTION


