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

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                                                       Function Definitions
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// Function that will get a histogram and turn its edges into an outline for plotting a phase space comparison. 
std::vector<TLine*> MakeTLineVector(TH2D *histogram){

    std::vector<TLine*> vLine; // The final vector containing the phase space. 
    
    // Get the bin max values to loop over. 
    int YRange{histogram->GetNbinsY()};
    int XRange{histogram->GetNbinsX()};

    // Positions to draw the line. // Initializing. 
    double x0{histogram->GetXaxis()->GetBinLowEdge(1)};
    double x1{histogram->GetYaxis()->GetBinLowEdge(1)}; 
    double y0{histogram->GetYaxis()->GetBinLowEdge(1)};
    double y1{histogram->GetYaxis()->GetBinLowEdge(1)}; 

    // Loop over the histogram elements
    for (int y = 1; y < YRange + 2 ; y++){ // Loop over each row I havent figured out yet why the plus 2, but it works ?!

        // Create a Tline with current x1, y0 and new y. Skip the very first go as need initial x1 to be set. 
        if (y != 1 ) {
            y1 = histogram->GetYaxis()->GetBinLowEdge(y);
            TLine *l  = new TLine(x0, y0, x1, y1 );  vLine.push_back(l); // Add vertical Tline
            y0 = y1;
        }

        for (int x = 1; x < XRange; x++){ // loop over each element in a row

            if (histogram->GetBinContent(x,y) == 0){
                x1 = histogram->GetXaxis()->GetBinLowEdge(x);
                TLine *l2  = new TLine(x0, y0, x1, y1 );  vLine.push_back(l2); // Add horizontal TLine. 
                x0 = x1;
                break;
            }

        } // End loop over elements in a row
    } // End loop over the rows

    return vLine;

}

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

// Function that plots 2 th2d histograms with same plot with the second histogram being outlined on top to map a phase space. (Histogram variable name, second histogram output file path/name to plot on top,scale factor, label for legend )
void Plot2DHistSAME (TH2D *histogram, TH2D *histogram2, const char * print_name, double scalefactor, std::string label ){
	
    auto *c1 = new TCanvas(); // Create a TCanvas
	c1->cd();

	histogram->Draw("COLZ TEXT00"); // Draw hist and set the options

    std::vector<TLine*> vLine = MakeTLineVector(histogram2); // Make the second histogram into an outline to plot on top

    // Loop over all lines in the TLine vecotr and draw on top of graph for phase space comparisons. 
    for (unsigned int i =0; i < vLine.size(); i++){
       
        vLine[i]->SetLineColor(kBlack); // Line specifiers
        vLine[i]->SetLineWidth(2);
        vLine[i]->Draw("SAME");
    }

    // Define a text box for the numi phase space label
    TLatex latex;
    latex.SetTextSize(0.05);
    if (label == "NuMI") latex.DrawLatex(6,120,"#Box NuMI Phase Space");
    else latex.DrawLatex(6,120,"#Box BNB Phase Space");

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

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                                                       Main Function
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
    TH2D*   hNue_E_vs_Theta_NuMI = new TH2D("Nue_E_vs_Theta_All_NuMI"," ",15., 0., 10. , 10., 0., 180);
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
	TTreeReaderValue<double> NueEnergyRV(myReader, "NueEnergy"); 	    // Load in the NueEnergy variable
    TTreeReaderValue<double> NueThetaRV(myReader, "NueTheta"); 	        // Load in the NueTheta variable
    TTreeReaderValue<double> NuePhiRV(myReader, "NuePhi"); 	            // Load in the NuePhi variable

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

    Plot2DHist(hNue_E_vs_Theta,    "plots/BNB/Nue_Energy_vs_Theta.png",              "COLZ,TEXT00", scalefactor );
    Plot2DHist(hNue_E_vs_Phi,      "plots/BNB/Nue_Energy_vs_Phi.png",                "COLZ,TEXT00", scalefactor );
    Plot2DHist(hNue_E_vs_Theta_EB, "plots/BNB/Nue_Energy_vs_Theta_Even_Binning.png", "COLZ,TEXT00", scalefactor );
    Plot2DHist(hNue_E_vs_Phi_EB,   "plots/BNB/Nue_Energy_vs_Phi_Even_Binning.png",   "COLZ,TEXT00", scalefactor );

    Plot1DHist(hNue_Energy, "plots/BNB/Nue_Energy.png", "HIST", scalefactor );
    Plot1DHist(hNue_Theta,  "plots/BNB/Nue_Theta.png",  "HIST", scalefactor ); 
    Plot1DHist(hNue_Phi,    "plots/BNB/Nue_Phi.png",    "HIST", scalefactor );

    // NuMi Nue
    Plot2DHist(hNue_E_vs_Theta_NuMI, "plots/NuMI/Nue_Energy_vs_Theta_NuMI.png", "COLZ, TEXT00", scalefactor );
    Plot2DHist(hNue_E_vs_Phi_NuMI, "plots/NuMI/Nue_Energy_vs_Phi_NuMI.png", "COLZ, TEXT00",     scalefactor );

    // Electron BNB Only
    Plot3DHist(helectron_E_vs_Theta_vs_Phi, "plots/BNB/El_Energy_vs_Theta_vs_Phi.png", "LEGO2",  scalefactor );

    Plot2DHist(helectron_E_vs_Theta,   "plots/BNB/El_Energy_vs_Theta.png",              "COLZ,TEXT00", scalefactor );
    Plot2DHist(helectron_E_vs_Phi,     "plots/BNB/El_Energy_vs_Phi.png",                "COLZ,TEXT00", scalefactor );
    Plot2DHist(helectron_E_vs_Theta_EB,"plots/BNB/El_Energy_vs_Theta_Even_Binning.png", "COLZ,TEXT00", scalefactor );
    Plot2DHist(helectron_E_vs_Phi_EB,  "plots/BNB/El_Energy_vs_Phi_Even_Binning.png",   "COLZ,TEXT00", scalefactor );

    Plot1DHist(helectron_Energy, "plots/BNB/El_Energy.png", "HIST", scalefactor );
    Plot1DHist(helectron_Theta,  "plots/BNB/El_Theta.png",  "HIST", scalefactor ); 
    Plot1DHist(helectron_Phi,    "plots/BNB/El_Phi.png",    "HIST", scalefactor );

     // NuMi Electron
    Plot2DHist(helectron_E_vs_Theta_NuMI, "plots/NuMI/El_Energy_vs_Theta_NuMI.png", "COLZ, TEXT00", scalefactor );
    Plot2DHist(helectron_E_vs_Phi_NuMI, "plots/NuMI/El_Energy_vs_Phi_NuMI.png", "COLZ, TEXT00",     scalefactor );

    // BNB and NuMI Plots for Nue (Superimposed)
    
    // NuMI on BNB
    Plot2DHistSAME(hNue_E_vs_Theta, hNue_E_vs_Theta_NuMI, "plots/BNB/Nue_Energy_vs_Theta_Overlaid.eps", scalefactor, "NuMI");
    Plot2DHistSAME(hNue_E_vs_Phi,   hNue_E_vs_Phi_NuMI,    "plots/BNB/Nue_Energy_vs_Phi_Overlaid.eps",  scalefactor, "NuMI" );

    Plot2DHistSAME(helectron_E_vs_Theta, helectron_E_vs_Theta_NuMI, "plots/BNB/El_Energy_vs_Theta_Overlaid.eps", scalefactor, "NuMI" );
    Plot2DHistSAME(helectron_E_vs_Phi,   helectron_E_vs_Phi_NuMI,    "plots/BNB/El_Energy_vs_Phi_Overlaid.eps",  scalefactor, "NuMI" );

    // BNB on NuMI
    Plot2DHistSAME( hNue_E_vs_Theta_NuMI, hNue_E_vs_Theta, "plots/NuMI/Nue_Energy_vs_Theta_Overlaid.eps", scalefactor, "BNB");
    Plot2DHistSAME( hNue_E_vs_Phi_NuMI,   hNue_E_vs_Phi,    "plots/NuMI/Nue_Energy_vs_Phi_Overlaid.eps",  scalefactor, "BNB" );

    Plot2DHistSAME( helectron_E_vs_Theta_NuMI, helectron_E_vs_Theta, "plots/NuMI/El_Energy_vs_Theta_Overlaid.eps", scalefactor, "BNB" );
    Plot2DHistSAME( helectron_E_vs_Phi_NuMI,   helectron_E_vs_Phi,   "plots/NuMI/El_Energy_vs_Phi_Overlaid.eps",   scalefactor, "BNB" );

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


