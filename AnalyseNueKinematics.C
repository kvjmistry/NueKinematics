// Program that loads in the ttree created from the NueKinematics Module and makes some plots as does the Larsoft Module. 
// WARNING: make sure a plots folder exits in the current directory otherwise the histogram will not get saved and get an error. 
// Execute over the NueKinematics.root file which has been run over the larsoft module NueKinematics_module.cc
// This file can be run by simply typing root AnalyseNueKinematics.C

#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TEfficiency.h"

#include <iostream>
#include <fstream>

// Function to make a Th1D histogram (Histogram variable name, output file path/name.extension, options for histogram)
void Plot1DHist (TH1D *histogram, const char * print_name, const char *Options){
	auto *c1 = new TCanvas(); // Create a TCanvas
	c1->cd();
	
    histogram->Draw();       // Draw hist and set the options
    histogram->SetOption(Options);

	c1->Print(print_name);   // Save the histogram
    c1->Close();             // Close the Canvas
}

// Function to make a Th2D histogram (Histogram variable name, output file path/name, options for histogram)
void Plot2DHist (TH2D *histogram, const char * print_name, const char *Options){
	auto *c1 = new TCanvas(); // Create a TCanvas
	c1->cd();

	histogram->Draw(); // Draw hist and set the options
    histogram->SetOption(Options);
	
    c1->Print(print_name); // Save the histogram
    c1->Close(); // Close the Canvas
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

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                  Initialize
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    TFile *MyFile = new TFile("plots/Plots_NeuKinematics.root","RECREATE");
    if ( MyFile->IsOpen() ) printf("File opened successfully\n");

    // Create Histograms
    // Nue All
    TH3D*   hNue_E_vs_Theta_vs_Phi = new TH3D("Nue_E_vs_Theta_vs_Phi_All","Nue_E_vs_Theta_All; Energy [GeV]; Theta [degrees]; Phi [degrees]", 20., 0., 10. , 10., 0., 180, 10., -180., 180 );

    TH2D*   hNue_E_vs_Theta = new TH2D("Nue_E_vs_Theta_All","Nue_E_vs_Theta_All; Energy [GeV]; Theta [degrees]",15., 0., 7. , 10., 0., 180);
    TH2D*   hNue_E_vs_Phi   = new TH2D("Nue_E_vs_Phi_All","Nue_E_vs_Phi_All; Energy [GeV]; Phi [degrees]",10., 0., 10. , 10., -180., 180);

    TH1D* 	hNue_Energy = new TH1D("Nue_Energy_All","Nue_Energy_All; E [GeV]; Events",50., 0., 5);
    TH1D* 	hNue_Theta =  new TH1D("Nue_Theta_All","Nue_Theta_All; Theta [Degrees]; Events", 100., 0., 180);
    TH1D* 	hNue_Phi =    new TH1D("Nue_Phi_All","Nue_Phi_All; Phi [Degrees]; Events", 10., -180., 180);

    

    // Electron all
    TH3D*   helectron_E_vs_Theta_vs_Phi = new TH3D("electron_E_vs_Theta_vs_Phi_All","electron_E_vs_Theta_All; Energy [GeV]; Theta [degrees]; Phi [degrees]", 20., 0., 10. , 10., 0., 180, 10., -180., 180 );

    TH2D*   helectron_E_vs_Theta = new TH2D("electron_E_vs_Theta_All","electron_E_vs_Theta_All; Energy [GeV]; Theta [degrees]",20., 0., 10. , 10., 0., 180);
    TH2D*   helectron_E_vs_Phi   = new TH2D("electron_E_vs_Phi_All","electron_E_vs_Phi_All; Energy [GeV]; Phi [degrees]",20., 0., 10. , 10., -180., 180);
  
    TH1D* 	helectron_Energy = new TH1D("electron_Energy_All","electron_Energy_All; E [GeV]; Events",50., 0., 5.);
    TH1D* 	helectron_Theta  = new TH1D("electron_Theta_All","electron_Theta_All; Theta [Degrees]; Events", 10., 0., 180);
    TH1D* 	helectron_Phi    = new TH1D("electron_Phi_All","electron_Phi_All; Phi [Degrees]; Events", 10., -180., 180);

    

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
    //                                                  Make Histograms
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // Nue
    Plot3DHist(hNue_E_vs_Theta_vs_Phi, "plots/Nue_Energy_vs_Theta_vs_Phi.png", "LEGO2");

    Plot2DHist(hNue_E_vs_Theta, "plots/Nue_Energy_vs_Theta.png", "COLZ,TEXT00" );
    Plot2DHist(hNue_E_vs_Phi,   "plots/Nue_Energy_vs_Phi.png",   "COLZ,TEXT00" );
    Plot1DHist(hNue_Energy,     "plots/Nue_Energy.png",          "HIST,TEXT00" );
    Plot1DHist(hNue_Theta,      "plots/Nue_Theta.png",           "HIST,TEXT00" ); 
    Plot1DHist(hNue_Phi,        "plots/Nue_Phi.png",             "HIST,TEXT00" );

    
    // Electron
    Plot3DHist(helectron_E_vs_Theta_vs_Phi, "plots/El_Energy_vs_Theta_vs_Phi.png", "LEGO2" );

    Plot2DHist(helectron_E_vs_Theta,"plots/El_Energy_vs_Theta.png", "COLZ,TEXT00" );
    Plot2DHist(helectron_E_vs_Phi,  "plots/El_Energy_vs_Phi.png",   "COLZ,TEXT00" );
    Plot1DHist(helectron_Energy,    "plots/El_Energy.png",          "HIST,TEXT00" );
    Plot1DHist(helectron_Theta,     "plots/El_Theta.png",           "HIST,TEXT00" ); 
    Plot1DHist(helectron_Phi,       "plots/El_Phi.png",             "HIST,TEXT00" );


    MyFile->Write(); // Save to a root file 

    gSystem->Exit(1); // Quit ROOT
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                     END
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

} // END MAIN FUNCTION


