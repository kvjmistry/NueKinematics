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

// Function that plots 2 th2d histograms with differnt pallete colours on the same plot(Histogram variable name, second histogram output file path/name, vecotor of TLines for phase space. )
void Plot2DHistSAME (TH2D *histogram, TH2D *histogram2, const char * print_name, std::vector<TLine*> vLine ){
	
    auto *c1 = new TCanvas(); // Create a TCanvas
	c1->cd();

	histogram->Draw("COLZ"); // Draw hist and set the options

    // Loop over all lines in the TLine vecotr and draw on top of graph for phase space comparisons. 
    for (unsigned int i =0; i < vLine.size(); i++){
       
        vLine[i]->SetLineColor(kBlack); // Line specifiers
        vLine[i]->SetLineWidth(2);
        vLine[i]->Draw("SAME");
    }

    // Define a text box for the numi phase space label
    TLatex latex;
    latex.SetTextSize(0.05);
    latex.DrawLatex(6,120,"#Box NuMI Phase Space");

    c1->Print(print_name);  // Save the histogram
    c1->Close();            // Close the Canvas

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
    TH3D*   hNue_E_vs_Theta_vs_Phi = new TH3D("Nue_E_vs_Theta_vs_Phi_All","Nue_E_vs_Theta_vs_Phi_All; Energy [GeV]; Theta [degrees]; Phi [degrees]", 20., 0., 10. , 10., 0., 180, 10., -180., 180 );

    TH2D*   hNue_E_vs_Theta = new TH2D("Nue_E_vs_Theta_All","Nue_E_vs_Theta_All; Energy [GeV]; Theta [degrees]",15., 0., 10. , 10., 0., 180);
    TH2D*   hNue_E_vs_Phi   = new TH2D("Nue_E_vs_Phi_All","Nue_E_vs_Phi_All; Energy [GeV]; Phi [degrees]",10., 0., 10. , 10., -180., 180);

    TH1D* 	hNue_Energy = new TH1D("Nue_Energy_All","Nue_Energy_All; E [GeV]; Events",50., 0., 5);
    TH1D* 	hNue_Theta =  new TH1D("Nue_Theta_All","Nue_Theta_All; Theta [Degrees]; Events", 100., 0., 1);
    TH1D* 	hNue_Phi =    new TH1D("Nue_Phi_All","Nue_Phi_All; Phi [Degrees]; Events", 10., -180., 180);

    
    // Electron all BNB
    TH3D*   helectron_E_vs_Theta_vs_Phi = new TH3D("electron_E_vs_Theta_vs_Phi_All","electron_E_vs_Theta_vs_Phi_All; Energy [GeV]; Theta [degrees]; Phi [degrees]", 20., 0., 10. , 10., 0., 180, 10., -180., 180 );

    TH2D*   helectron_E_vs_Theta = new TH2D("electron_E_vs_Theta_All","electron_E_vs_Theta_All; Energy [GeV]; Theta [degrees]",20., 0., 10. , 10., 0., 180);
    TH2D*   helectron_E_vs_Phi   = new TH2D("electron_E_vs_Phi_All","electron_E_vs_Phi_All; Energy [GeV]; Phi [Rad]",20., 0., 10. , 10., -3.145, 3.145);
  
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
        helectron_Phi->Fill(*ElectronPhiRV * ( 3.1415 / 180. ));

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
        helectron_E_vs_Phi_NuMI->Fill(*NuMI_ElectronEnergy_RV, *NuMI_ElectronPhi_RV * (180 / 3.1415));

        // helectron_Energy->Fill(*ElectronEnergyRV);
        // helectron_Theta->Fill(*ElectronThetaRV);
        // helectron_Phi->Fill(*ElectronPhiRV);

	}

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                        Define the lines for coltons numi phase space
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // 
    
    // Nue E vs Theta NuMI
    std::vector<TLine*> vLine_Nue_E_vs_Theta;

    TLine *l_theta_1  = new TLine(0., 18.0, 6.67, 18.0);     vLine_Nue_E_vs_Theta.push_back(l_theta_1);
    TLine *l_theta_2  = new TLine(6.67, 18.0, 6.67, 35.9);   vLine_Nue_E_vs_Theta.push_back(l_theta_2);
    TLine *l_theta_3  = new TLine(2.1, 35.9, 6.67, 35.9);    vLine_Nue_E_vs_Theta.push_back(l_theta_3);
    TLine *l_theta_4  = new TLine(2.1, 35.9, 2.1, 54.1);     vLine_Nue_E_vs_Theta.push_back(l_theta_4);
    TLine *l_theta_5  = new TLine(1.1, 54.1, 2.1, 54.1);     vLine_Nue_E_vs_Theta.push_back(l_theta_5);
    TLine *l_theta_6  = new TLine(1.1, 54.1, 1.1, 72.1);     vLine_Nue_E_vs_Theta.push_back(l_theta_6);
    TLine *l_theta_7  = new TLine(0.53, 72.1, 1.1, 72.1);    vLine_Nue_E_vs_Theta.push_back(l_theta_7);
    TLine *l_theta_8  = new TLine(0.53, 72.1, 0.53, 125.9);  vLine_Nue_E_vs_Theta.push_back(l_theta_8);
    TLine *l_theta_9  = new TLine(0.0, 125.9, 0.53, 125.9);  vLine_Nue_E_vs_Theta.push_back(l_theta_9);
    TLine *l_theta_10 = new TLine(7.4, 18.0, 8.0, 18.0);     vLine_Nue_E_vs_Theta.push_back(l_theta_10);
    TLine *l_theta_11 = new TLine(8.0, 18.0, 8.0, 35.9);     vLine_Nue_E_vs_Theta.push_back(l_theta_11);
    TLine *l_theta_12 = new TLine(7.4, 35.9, 8.0, 35.9);     vLine_Nue_E_vs_Theta.push_back(l_theta_12);
    TLine *l_theta_13 = new TLine(7.4,18.0, 7.4, 35.9);      vLine_Nue_E_vs_Theta.push_back(l_theta_13);

    // Nue E vs Phi NuMI
    std::vector<TLine*> vLine_Nue_E_vs_Phi;
    TLine *l_phi_1  = new TLine(0., 0., 8.0, 0.0);        vLine_Nue_E_vs_Phi.push_back(l_phi_1);
    TLine *l_phi_2  = new TLine(8.0, 0.0, 8.0, 36.1);     vLine_Nue_E_vs_Phi.push_back(l_phi_2);
    TLine *l_phi_3  = new TLine(0.0, 36.1, 8.0, 36.1);    vLine_Nue_E_vs_Phi.push_back(l_phi_3);

    // El E vs Theta
    std::vector<TLine*> vLine_El_E_vs_Theta;

    TLine *l_el_theta_1  = new TLine(0, 0, 4, 0);               vLine_El_E_vs_Theta.push_back(l_el_theta_1);
    TLine *l_el_theta_2  = new TLine(4, 0, 4 ,18);              vLine_El_E_vs_Theta.push_back(l_el_theta_2);
    TLine *l_el_theta_3  = new TLine(4, 18, 5.5, 18);           vLine_El_E_vs_Theta.push_back(l_el_theta_3);
    TLine *l_el_theta_4  = new TLine(5.5, 18, 5.5, 36.0);       vLine_El_E_vs_Theta.push_back(l_el_theta_4);
    TLine *l_el_theta_5  = new TLine(5.5, 36, 3.5, 36.0);       vLine_El_E_vs_Theta.push_back(l_el_theta_5);
    TLine *l_el_theta_6  = new TLine(3.5, 36.0, 3.5, 53.8);     vLine_El_E_vs_Theta.push_back(l_el_theta_6);
    TLine *l_el_theta_7  = new TLine(3.5, 53.8, 2.0, 53.8);     vLine_El_E_vs_Theta.push_back(l_el_theta_7);
    TLine *l_el_theta_8  = new TLine(2.0, 53.8, 2.0, 71.9);     vLine_El_E_vs_Theta.push_back(l_el_theta_8);
    TLine *l_el_theta_9  = new TLine(2.0, 71.9, 1.5, 71.9);     vLine_El_E_vs_Theta.push_back(l_el_theta_9);
    TLine *l_el_theta_10 = new TLine(1.5, 71.9, 1.5, 90.0);     vLine_El_E_vs_Theta.push_back(l_el_theta_10);
    TLine *l_el_theta_11 = new TLine(1.5, 90.0, 1 ,90);         vLine_El_E_vs_Theta.push_back(l_el_theta_11);
    TLine *l_el_theta_12 = new TLine(1, 90, 1, 126);            vLine_El_E_vs_Theta.push_back(l_el_theta_12);
    TLine *l_el_theta_13 = new TLine(1, 126, 0.5, 126);         vLine_El_E_vs_Theta.push_back(l_el_theta_13);
    TLine *l_el_theta_14 = new TLine(0.5, 126, 0.5, 180);       vLine_El_E_vs_Theta.push_back(l_el_theta_14);
    TLine *l_el_theta_15 = new TLine(6.0, 18, 7, 18);           vLine_El_E_vs_Theta.push_back(l_el_theta_15);
    TLine *l_el_theta_16 = new TLine(7, 18, 7, 36);             vLine_El_E_vs_Theta.push_back(l_el_theta_16);
    TLine *l_el_theta_17 = new TLine(7, 36, 6, 36);             vLine_El_E_vs_Theta.push_back(l_el_theta_17);
    TLine *l_el_theta_18 = new TLine(6,36,6,18);                vLine_El_E_vs_Theta.push_back(l_el_theta_18);

    // Nue E vs Phi NuMI
    std::vector<TLine*> vLine_El_E_vs_Phi;
    TLine *l_el_phi_1  = new TLine(0, -35.6, 5.5, -35.6);        vLine_El_E_vs_Phi.push_back(l_el_phi_1);
    TLine *l_el_phi_2  = new TLine(5.5, -35.6, 5.5, 0);     vLine_El_E_vs_Phi.push_back(l_el_phi_2);
    TLine *l_el_phi_3  = new TLine(5.5, 0 ,5, 0);    vLine_El_E_vs_Phi.push_back(l_el_phi_3);

    TLine *l_el_phi_4  = new TLine(5, 0, 5, 35.6);    vLine_El_E_vs_Phi.push_back(l_el_phi_4);
    TLine *l_el_phi_5  = new TLine(5, 35.6, 0 ,35.6);    vLine_El_E_vs_Phi.push_back(l_el_phi_5);
    TLine *l_el_phi_6  = new TLine(6, 0 , 7, 0);    vLine_El_E_vs_Phi.push_back(l_el_phi_6);
    TLine *l_el_phi_7  = new TLine(6, 35.6, 7, 35.6);    vLine_El_E_vs_Phi.push_back(l_el_phi_7);
    TLine *l_el_phi_8  = new TLine(6, 0, 6, 35.6);    vLine_El_E_vs_Phi.push_back(l_el_phi_8);
    TLine *l_el_phi_9  = new TLine(7, 0, 7, 35.6);    vLine_El_E_vs_Phi.push_back(l_el_phi_9);


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
    Plot2DHistSAME(hNue_E_vs_Theta , hNue_E_vs_Theta_NuMI, "plots/Nue_Energy_vs_Theta_Overlaid.png", vLine_Nue_E_vs_Theta );
    Plot2DHistSAME(hNue_E_vs_Phi , hNue_E_vs_Phi_NuMI, "plots/Nue_Energy_vs_Phi_Overlaid.png", vLine_Nue_E_vs_Phi );

    Plot2DHistSAME(helectron_E_vs_Theta , helectron_E_vs_Theta_NuMI, "plots/El_Energy_vs_Theta_Overlaid.png", vLine_El_E_vs_Theta );
    Plot2DHistSAME(helectron_E_vs_Phi , helectron_E_vs_Phi_NuMI, "plots/El_Energy_vs_Phi_Overlaid.png", vLine_El_E_vs_Phi );


    MyFile->Write(); // Save to a root file 
    MyFile->Close();
    gSystem->Exit(1); // Quit ROOT

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                     END
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

} // END MAIN FUNCTION


