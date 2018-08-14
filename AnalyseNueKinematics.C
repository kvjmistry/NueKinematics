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
void Plot1DHist (TH1D *histogram, const char * print_name, const char *Options, double scalefactor){
	auto *c1 = new TCanvas(); // Create a TCanvas
	c1->cd();
	
    histogram->Draw();       // Draw hist and set the options
    histogram->SetOption(Options);
    histogram->SetLineWidth(2);

    histogram->Scale(scalefactor); // Scale the histogram to the number of NuMI events

	c1->Print(print_name);   // Save the histogram
    c1->Close();             // Close the Canvas
}

// Function to make a Th2D histogram (Histogram variable name, output file path/name, options for histogram)
void Plot2DHist (TH2D *histogram, const char * print_name, const char *Options, double scalefactor){
	
    auto *c1 = new TCanvas(); // Create a TCanvas
	c1->cd();

    gStyle->SetPalette(kBird,0,1);
    histogram->SetLineWidth(1);
    histogram->SetLineColor(kBlack);


    histogram->Draw(); // Draw hist and set the options
    histogram->SetOption(Options);

    histogram->Scale(scalefactor); // Scale the histogram to the number of NuMI events
    
	
    c1->Print(print_name); // Save the histogram
    c1->Close(); // Close the Canvas
}

// Function that plots 2 th2d histograms with differnt pallete colours on the same plot(Histogram variable name, second histogram output file path/name, vecotor of TLines for phase space. )
void Plot2DHistSAME (TH2D *histogram, TH2D *histogram2, const char * print_name, std::vector<TLine*> vLine, double scalefactor ){
	
    auto *c1 = new TCanvas(); // Create a TCanvas
	c1->cd();

	histogram->Draw("COLZ TEXT00"); // Draw hist and set the options

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

    histogram->SetMarkerSize(0.6);

    //histogram->Scale(scalefactor); // Scale the histogram to the number of NuMI events

    c1->Print(print_name);  // Save the histogram
    c1->Close();            // Close the Canvas

    gStyle->SetPalette(kBird); // Reset back to the default colour pallete. 
}

// Function to make a Th2D histogram (Histogram variable name, output file path/name, options for histogram)
void Plot3DHist (TH3D *histogram, const char * print_name, const char *Options, double scalefactor){
	auto *c1 = new TCanvas(); // Create a TCanvas
	c1->cd();

	histogram->Draw(); // Draw hist and set the options
    histogram->SetOption(Options);
	
    histogram->Scale(scalefactor); // Scale the histogram to the number of NuMI events

    c1->Print(print_name); // Save the histogram
    c1->Close(); // Close the Canvas
}


// Main Function
void AnalyseNueKinematics() {
    // General formatting
    gStyle->SetOptStat(0); // Stats Box
    gStyle->SetPaintTextFormat("1.1f"); 

    // Counters for Number of events
    double BNB_Counter{0}, NuMI_Counter{0};


    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                       Define a variable bin widths for the histograms.
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    

    // ---- Nue E vs Phi -----
    Float_t bins_Phi_Nue[] = { -180, -140, -100, -45, 0, 45, 90, 140, 180 }; // The bins for Phi in the Nue E vs Phi Plot
    Int_t   binnum_Phi_Nue = sizeof(bins_Phi_Nue)/sizeof(Float_t) - 1; 

    Float_t bins_EPhi_Nue[] = { 0, 0.67,  1, 1.4, 2, 6, 10}; // The bins for E in the Nue E vs Phi Plot
    Int_t   binnum_EPhi_Nue = sizeof(bins_EPhi_Nue)/sizeof(Float_t) - 1; // As above, but for the bin numbers. 
    // ---- ------------ -----

    // ---- Nue E vs Theta -----
    Float_t bins_Theta_Nue[] = { 0,  45,  180 }; // The bins for Theta in the Nue E vs Theta Plot
    Int_t   binnum_Theta_Nue = sizeof(bins_Theta_Nue)/sizeof(Float_t) - 1; 

    Float_t bins_ETheta_Nue[] = { 0,  0.3, 0.4, 0.48, 0.52, 0.56, 0.6, 0.64, 0.66, 0.7, 0.725,  0.75, 0.8, 0.84, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98,  1, 1.05,  1.1, 1.15, 1.2, 1.22, 1.24, 1.26, 1.28, 1.3,1.34, 1.38, 1.42, 1.48, 1.54, 1.6,1.64, 1.68, 1.74, 1.8, 1.85, 1.9, 2, 2.1, 2.25, 2.3, 2.35, 2.5, 2.6, 2.85, 3.2, 3.5, 6, 10}; // The bins for E in the Nue E vs Theta Plot
    Int_t   binnum_ETheta_Nue = sizeof(bins_ETheta_Nue)/sizeof(Float_t) - 1; // As above, but for the bin numbers. 
    // ---- ------------ -----

    // ---- Nue E vs Phi -----
    Float_t bins_Phi_El[] = { -180, -135, -90, -45, 0, 45, 90, 135, 180 }; // The bins for Phi in the El E vs Phi Plot
    Int_t   binnum_Phi_El = sizeof(bins_Phi_El)/sizeof(Float_t) - 1; 

    Float_t bins_EPhi_El[] = { 0, 0.3,  0.53, 0.82, 1.3, 6, 10}; // The bins for E in the El E vs Phi Plot
    Int_t   binnum_EPhi_El = sizeof(bins_EPhi_El)/sizeof(Float_t) - 1; // As above, but for the bin numbers. 
    // ---- ------------ -----

    // ---- Nue E vs Theta -----
    Float_t bins_Theta_El[] = { 0, 30, 50, 65, 100,  180 }; // The bins for Theta in the Nue E vs Theta Plot
    Int_t   binnum_Theta_El = sizeof(bins_Theta_El)/sizeof(Float_t) - 1; 

    Float_t bins_ETheta_El[] = { 0, 0.5, 0.65, 0.77, 0.9, 1, 1.1, 1.2, 1.3,1.4, 1.5, 1.65,  1.83,2.15, 2.6, 5.5, 10 }; // The bins for E in the El E vs Theta Plot
    Int_t   binnum_ETheta_El = sizeof(bins_ETheta_El)/sizeof(Float_t) - 1; // As above, but for the bin numbers. 
    // ---- ------------ -----


    // ++++++++++++++++++++++++++++++++++s++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                  Initialize
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    TFile *MyFile = new TFile("plots/Plots_NeuKinematics.root","RECREATE");
    if ( MyFile->IsOpen() ) printf("File opened successfully\n");

    // Create Histograms
    // Nue All BNB
    TH3D*   hNue_E_vs_Theta_vs_Phi = new TH3D("Nue_E_vs_Theta_vs_Phi_All","Nue_E_vs_Theta_vs_Phi_All; Energy [GeV]; Theta [degrees]; Phi [degrees]", 20., 0., 10. , 10., 0., 180, 10., -180., 180 );

    TH2D*   hNue_E_vs_Theta = new TH2D("Nue_E_vs_Theta_All","Nue_E_vs_Theta_All; Energy [GeV]; Theta [degrees]",15., 0., 10. , 10., 0., 180);  // Fixed Vals 15., 0., 10. , 10., 0., 180
    TH2D*   hNue_E_vs_Phi   = new TH2D("Nue_E_vs_Phi_All","Nue_E_vs_Phi_All; Energy [GeV]; Phi [degrees]", 10., 0., 10. , 10., -180., 180); // Fixed Values 10., 0., 10. , 10., -180., 180
    TH2D*   hNue_E_vs_Theta_EB = new TH2D("Nue_E_vs_Theta_All_EB","Nue_E_vs_Theta_All; Energy [GeV]; Theta [degrees]",binnum_ETheta_Nue ,bins_ETheta_Nue, binnum_Theta_Nue ,bins_Theta_Nue);  // Even binned histograms
    TH2D*   hNue_E_vs_Phi_EB    = new TH2D("Nue_E_vs_Phi_All_EB","Nue_E_vs_Phi_All; Energy [GeV]; Phi [degrees]", binnum_EPhi_Nue, bins_EPhi_Nue , binnum_Phi_Nue, bins_Phi_Nue); // Even binned histograms
    
    TH1D* 	hNue_Energy = new TH1D("Nue_Energy_All","Nue_Energy_All; E [GeV]; Events",50., 0., 5);
    TH1D* 	hNue_Theta =  new TH1D("Nue_Theta_All","Nue_Theta_All; Theta [Degrees]; Events", 100., 0., 1);
    TH1D* 	hNue_Phi =    new TH1D("Nue_Phi_All","Nue_Phi_All; Phi [Degrees]; Events", 10., -180., 180);

    
    // Electron all BNB
    TH3D*   helectron_E_vs_Theta_vs_Phi = new TH3D("electron_E_vs_Theta_vs_Phi_All","electron_E_vs_Theta_vs_Phi_All; Energy [GeV]; Theta [degrees]; Phi [degrees]", 20., 0., 10. , 10., 0., 180, 10., -180., 180 );

    TH2D*   helectron_E_vs_Theta = new TH2D("electron_E_vs_Theta_All","electron_E_vs_Theta_All; Energy [GeV]; Theta [degrees]",20., 0., 10. , 10., 0., 180);
    TH2D*   helectron_E_vs_Phi   = new TH2D("electron_E_vs_Phi_All","electron_E_vs_Phi_All; Energy [GeV]; Phi [degrees]",20., 0., 10. , 10., -180., 180  ); // Fixed Values 20., 0., 10. , 10., -180., 180
    TH2D*   helectron_E_vs_Theta_EB = new TH2D("electron_E_vs_Theta_All_EB","electron_E_vs_Theta_All; Energy [GeV]; Theta [degrees]", binnum_ETheta_El, bins_ETheta_El, binnum_Theta_El, bins_Theta_El );  // Even binned histograms
    TH2D*   helectron_E_vs_Phi_EB   = new TH2D("electron_E_vs_Phi_All_EB","electron_E_vs_Phi_All; Energy [GeV]; Phi [degrees]", binnum_EPhi_El,bins_EPhi_El, binnum_Phi_El, bins_Phi_El  );                // Even binned histograms

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
        BNB_Counter++;

        // Fill Nue Histograms
        hNue_E_vs_Theta_vs_Phi ->Fill(*NueEnergyRV, *NueThetaRV, *NuePhiRV );

        hNue_E_vs_Theta->Fill(*NueEnergyRV, *NueThetaRV );
        hNue_E_vs_Phi->Fill(*NueEnergyRV, *NuePhiRV);

        hNue_E_vs_Theta_EB->Fill(*NueEnergyRV, *NueThetaRV ); // Also fill the even binned versions of the histograms. 
        hNue_E_vs_Phi_EB->Fill(*NueEnergyRV, *NuePhiRV);

		hNue_Energy->Fill(*NueEnergyRV); 
        hNue_Theta->Fill(*NueThetaRV); 
        hNue_Phi->Fill(*NuePhiRV); 
        
        // Fill Electron Histograms
        helectron_E_vs_Theta_vs_Phi ->Fill(*ElectronEnergyRV, *ElectronThetaRV, *ElectronPhiRV );

        helectron_E_vs_Theta->Fill(*ElectronEnergyRV, *ElectronThetaRV );
        helectron_E_vs_Phi->Fill(*ElectronEnergyRV, *ElectronPhiRV);

        helectron_E_vs_Theta_EB->Fill(*ElectronEnergyRV, *ElectronThetaRV ); // Also fill the even binned versions of the histograms. 
        helectron_E_vs_Phi_EB->Fill(*ElectronEnergyRV, *ElectronPhiRV);

        helectron_Energy->Fill(*ElectronEnergyRV);
        helectron_Theta->Fill(*ElectronThetaRV);
        helectron_Phi->Fill(*ElectronPhiRV );

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
        NuMI_Counter++;

        NuMI_theta = acos(*NuMI_NueTheta_RV) * (180 / 3.1415);

        // // Fill Nue Histograms
        // hNue_E_vs_Theta_vs_Phi ->Fill(*NueEnergyRV, *NueThetaRV, *NuePhiRV );

        hNue_E_vs_Theta_NuMI->Fill(*NuMI_NueEnergy_RV,  NuMI_theta );
        hNue_E_vs_Phi_NuMI->Fill(*NuMI_NueEnergy_RV, *NuMI_NuePhi_RV * (180 / 3.1415));

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

    TLine *l_theta_1  = new TLine(0., 18.0, 7.4, 18.0);     vLine_Nue_E_vs_Theta.push_back(l_theta_1);
    //TLine *l_theta_2  = new TLine(6.67, 18.0, 6.67, 35.9);   vLine_Nue_E_vs_Theta.push_back(l_theta_2);
    TLine *l_theta_3  = new TLine(2.1, 35.9, 7.4, 35.9);    vLine_Nue_E_vs_Theta.push_back(l_theta_3);
    TLine *l_theta_4  = new TLine(2.1, 35.9, 2.1, 54.1);     vLine_Nue_E_vs_Theta.push_back(l_theta_4);
    TLine *l_theta_5  = new TLine(1.1, 54.1, 2.1, 54.1);     vLine_Nue_E_vs_Theta.push_back(l_theta_5);
    TLine *l_theta_6  = new TLine(1.1, 54.1, 1.1, 72.1);     vLine_Nue_E_vs_Theta.push_back(l_theta_6);
    TLine *l_theta_7  = new TLine(0.53, 72.1, 1.1, 72.1);    vLine_Nue_E_vs_Theta.push_back(l_theta_7);
    TLine *l_theta_8  = new TLine(0.53, 72.1, 0.53, 125.9);  vLine_Nue_E_vs_Theta.push_back(l_theta_8);
    TLine *l_theta_9  = new TLine(0.0, 125.9, 0.53, 125.9);  vLine_Nue_E_vs_Theta.push_back(l_theta_9);
    TLine *l_theta_10 = new TLine(7.4, 18.0, 8.0, 18.0);     vLine_Nue_E_vs_Theta.push_back(l_theta_10);
    TLine *l_theta_11 = new TLine(8.0, 18.0, 8.0, 35.9);     vLine_Nue_E_vs_Theta.push_back(l_theta_11);
    TLine *l_theta_12 = new TLine(7.4, 35.9, 8.0, 35.9);     vLine_Nue_E_vs_Theta.push_back(l_theta_12);
    //TLine *l_theta_13 = new TLine(7.4,18.0, 7.4, 35.9);      vLine_Nue_E_vs_Theta.push_back(l_theta_13);

    // Nue E vs Phi NuMI
    std::vector<TLine*> vLine_Nue_E_vs_Phi;
    TLine *l_phi_1  = new TLine(0., 0., 8.0, 0.0);        vLine_Nue_E_vs_Phi.push_back(l_phi_1);
    TLine *l_phi_2  = new TLine(8.0, 0.0, 8.0, 36.1);     vLine_Nue_E_vs_Phi.push_back(l_phi_2);
    TLine *l_phi_3  = new TLine(1.0, 36.1, 8.0, 36.1);    vLine_Nue_E_vs_Phi.push_back(l_phi_3);
    TLine *l_phi_4  = new TLine(1.0, 36.1, 1.0, 72.6);    vLine_Nue_E_vs_Phi.push_back(l_phi_4);
    TLine *l_phi_5  = new TLine(0.0, 72.6, 1.0, 72.6);    vLine_Nue_E_vs_Phi.push_back(l_phi_5);

    // El E vs Theta
    std::vector<TLine*> vLine_El_E_vs_Theta;

    TLine *l_el_theta_1  = new TLine(0, 0, 4, 0);               vLine_El_E_vs_Theta.push_back(l_el_theta_1);
    TLine *l_el_theta_2  = new TLine(4, 0, 4 ,18);              vLine_El_E_vs_Theta.push_back(l_el_theta_2);
    TLine *l_el_theta_3  = new TLine(4, 18, 5.5, 18);           vLine_El_E_vs_Theta.push_back(l_el_theta_3);
    TLine *l_el_theta_4  = new TLine(5.5, 18, 6.0, 18.0);       vLine_El_E_vs_Theta.push_back(l_el_theta_4);
    TLine *l_el_theta_5  = new TLine(6.0, 36, 3.5, 36.0);       vLine_El_E_vs_Theta.push_back(l_el_theta_5);
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
    //TLine *l_el_theta_18 = new TLine(6,36,6,18);                vLine_El_E_vs_Theta.push_back(l_el_theta_18);

    // El E vs Phi NuMI
    std::vector<TLine*> vLine_El_E_vs_Phi;
    TLine *l_el_phi_1  = new TLine(0, -180, 1.5, -180);         vLine_El_E_vs_Phi.push_back(l_el_phi_1);
    TLine *l_el_phi_2  = new TLine(1.5, -180, 1.5 ,-144);       vLine_El_E_vs_Phi.push_back(l_el_phi_2);
    TLine *l_el_phi_3  = new TLine(1.5, -144, 2.0, -144);       vLine_El_E_vs_Phi.push_back(l_el_phi_3);
    TLine *l_el_phi_4  = new TLine(2.0, -144, 2.0, -108);       vLine_El_E_vs_Phi.push_back(l_el_phi_4);
    TLine *l_el_phi_5  = new TLine(2.0, -108, 2.5, -108);       vLine_El_E_vs_Phi.push_back(l_el_phi_5);
    TLine *l_el_phi_6  = new TLine(2.5, -108 , 2.5, -72);       vLine_El_E_vs_Phi.push_back(l_el_phi_6);
    TLine *l_el_phi_7  = new TLine(2.5, -72, 3.0, -72);         vLine_El_E_vs_Phi.push_back(l_el_phi_7);
    TLine *l_el_phi_8  = new TLine(3.0, -72, 3.0, -36);         vLine_El_E_vs_Phi.push_back(l_el_phi_8);
    TLine *l_el_phi_9  = new TLine(3.0, -36, 5.5, -36);         vLine_El_E_vs_Phi.push_back(l_el_phi_9);
    TLine *l_el_phi_10  = new TLine(5.5, -36, 5.5, 0);          vLine_El_E_vs_Phi.push_back(l_el_phi_10);
    TLine *l_el_phi_11  = new TLine(5.5, 0, 6.0, 0);            vLine_El_E_vs_Phi.push_back(l_el_phi_11);
    //TLine *l_el_phi_12  = new TLine(5.0, 0, 5.0, 36);           vLine_El_E_vs_Phi.push_back(l_el_phi_12);
    TLine *l_el_phi_13  = new TLine(6.0, 36, 3.0, 36);          vLine_El_E_vs_Phi.push_back(l_el_phi_13);
    TLine *l_el_phi_14  = new TLine(3.0, 36, 3.0, 72);          vLine_El_E_vs_Phi.push_back(l_el_phi_14);
    TLine *l_el_phi_15  = new TLine(3.0, 72, 2.0, 72);          vLine_El_E_vs_Phi.push_back(l_el_phi_15);
    TLine *l_el_phi_16  = new TLine(2.0, 72, 2.0, 108);         vLine_El_E_vs_Phi.push_back(l_el_phi_16);
    TLine *l_el_phi_17  = new TLine(2.0, 108, 1.5, 108);        vLine_El_E_vs_Phi.push_back(l_el_phi_17);
    TLine *l_el_phi_18  = new TLine(1.5, 108, 1.5, 180);        vLine_El_E_vs_Phi.push_back(l_el_phi_18);
    TLine *l_el_phi_19  = new TLine(6, 0, 7, 0);                vLine_El_E_vs_Phi.push_back(l_el_phi_19);
    TLine *l_el_phi_20  = new TLine(6, 36, 7, 36);              vLine_El_E_vs_Phi.push_back(l_el_phi_20);
    //TLine *l_el_phi_21  = new TLine(6, 0, 6, 36);               vLine_El_E_vs_Phi.push_back(l_el_phi_21);
    TLine *l_el_phi_22  = new TLine(7, 0, 7, 36);               vLine_El_E_vs_Phi.push_back(l_el_phi_22);



    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                  Make Histograms
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Define a scale factor to normalise histograms
    double scalefactor{NuMI_Counter / BNB_Counter};

    std::cout <<" ++++++++++++++++====================+++++++++++++++++++++++++++++++ \n" << std::endl;
    std::cout << "Scale Factor for BNB events:\t" << scalefactor  << std::endl;
    std::cout <<"\n ++++++++++++++=====================+++++++++++++++++++++++++++++++ " << std::endl;

    std::cout << "\nCreating Histograms....\n" << std::endl;

    // Nue BNB Only
    //Plot3DHist(hNue_E_vs_Theta_vs_Phi, "plots/BNB/Nue_Energy_vs_Theta_vs_Phi.png", "LEGO2", scalefactor);

    Plot2DHist(hNue_E_vs_Theta, "plots/BNB/Nue_Energy_vs_Theta.png", "COLZ,TEXT00", scalefactor );
    Plot2DHist(hNue_E_vs_Phi,   "plots/BNB/Nue_Energy_vs_Phi.png",   "COLZ,TEXT00", scalefactor );
    Plot2DHist(hNue_E_vs_Theta_EB, "plots/BNB/Nue_Energy_vs_Theta_Even_Binning.png", "COLZ,TEXT00", scalefactor );
    Plot2DHist(hNue_E_vs_Phi_EB,   "plots/BNB/Nue_Energy_vs_Phi_Even_Binning.png",   "COLZ,TEXT00", scalefactor );

    Plot1DHist(hNue_Energy,     "plots/BNB/Nue_Energy.png",          "HIST",        scalefactor );
    Plot1DHist(hNue_Theta,      "plots/BNB/Nue_Theta.png",           "HIST",        scalefactor ); 
    Plot1DHist(hNue_Phi,        "plots/BNB/Nue_Phi.png",             "HIST",        scalefactor );

    // NuMi Nue
    Plot2DHist(hNue_E_vs_Theta_NuMI, "plots/NuMI/Nue_Energy_vs_Theta_NuMI.png", "COLZ, TEXT00",     scalefactor );
    Plot2DHist(hNue_E_vs_Phi_NuMI, "plots/NuMI/Nue_Energy_vs_Phi_NuMI.png", "COLZ, TEXT00",         scalefactor );

    // Electron BNB Only
    Plot3DHist(helectron_E_vs_Theta_vs_Phi, "plots/BNB/El_Energy_vs_Theta_vs_Phi.png", "LEGO2",     scalefactor );

    Plot2DHist(helectron_E_vs_Theta,"plots/BNB/El_Energy_vs_Theta.png", "COLZ,TEXT00", scalefactor );
    Plot2DHist(helectron_E_vs_Phi,  "plots/BNB/El_Energy_vs_Phi.png",   "COLZ,TEXT00", scalefactor );
    Plot2DHist(helectron_E_vs_Theta_EB,"plots/BNB/El_Energy_vs_Theta_Even_Binning.png", "COLZ,TEXT00", scalefactor );
    Plot2DHist(helectron_E_vs_Phi_EB,  "plots/BNB/El_Energy_vs_Phi_Even_Binning.png",   "COLZ,TEXT00", scalefactor );

    Plot1DHist(helectron_Energy,    "plots/BNB/El_Energy.png",          "HIST",        scalefactor );
    Plot1DHist(helectron_Theta,     "plots/BNB/El_Theta.png",           "HIST",        scalefactor ); 
    Plot1DHist(helectron_Phi,       "plots/BNB/El_Phi.png",             "HIST",        scalefactor );

     // NuMi Electron
    Plot2DHist(helectron_E_vs_Theta_NuMI, "plots/NuMI/El_Energy_vs_Theta_NuMI.png", "COLZ, TEXT00", scalefactor );
    Plot2DHist(helectron_E_vs_Phi_NuMI, "plots/NuMI/El_Energy_vs_Phi_NuMI.png", "COLZ, TEXT00",     scalefactor );

    // BNB and NuMI Plots for Nue (Superimposed)
    Plot2DHistSAME(hNue_E_vs_Theta , hNue_E_vs_Theta_NuMI, "plots/Nue_Energy_vs_Theta_Overlaid.eps", vLine_Nue_E_vs_Theta , scalefactor);
    Plot2DHistSAME(hNue_E_vs_Phi , hNue_E_vs_Phi_NuMI, "plots/Nue_Energy_vs_Phi_Overlaid.eps", vLine_Nue_E_vs_Phi,          scalefactor );

    Plot2DHistSAME(helectron_E_vs_Theta , helectron_E_vs_Theta_NuMI, "plots/El_Energy_vs_Theta_Overlaid.eps", vLine_El_E_vs_Theta, scalefactor );
    Plot2DHistSAME(helectron_E_vs_Phi , helectron_E_vs_Phi_NuMI, "plots/El_Energy_vs_Phi_Overlaid.eps", vLine_El_E_vs_Phi,         scalefactor );

    std::cout <<" ++++++++++++++++====================+++++++++++++++++++++++++++++++ \n" << std::endl;
    std::cout <<"BNB Events:\t" << BNB_Counter << "\tNuMI Events:\t" << NuMI_Counter << std::endl;
    std::cout <<"\n ++++++++++++++=====================+++++++++++++++++++++++++++++++ " << std::endl;

    MyFile->Write(); // Save to a root file 
    MyFile->Close();

    

    gSystem->Exit(1); // Quit ROOT

    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                     END
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

} // END MAIN FUNCTION


