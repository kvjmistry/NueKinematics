// Program that loads in the ttree created from the NueKinematics Module and makes some plots as does the Larsoft Module. 
// WARNING: make sure a plots folder exits in the current directory otherwise the histogram will not get saved and get an error. 

#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TEfficiency.h"

#include <iostream>
#include <fstream>


void Plot1DHist (TH1D *histogram, const char * print_name){
	auto *c1 = new TCanvas();
	c1->cd();
	histogram->Draw();
    histogram->SetOption("HIST,TEXT00");

	c1->Print(print_name);
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
    TH2D*   hNue_E_vs_Theta = new TH2D("Nue_E_vs_Theta_All","Nue_E_vs_Theta_All; Energy [GeV]; Theta [degrees]",15., 0., 7. , 10., 0., 180);
    TH2D*   hNue_E_vs_Phi   = new   TH2D("Nue_E_vs_Phi_All","Nue_E_vs_Phi_All; Energy [GeV]; Phi [degrees]",10., 0., 10. , 10., -180., 180);

    TH1D* 	hNue_Energy = new TH1D("Nue_Energy_All","Nue_Energy_All; E [GeV]; Events",50., 0., 5);
    TH1D* 	hNue_Theta = new TH1D("Nue_Theta_All","Nue_Theta_All; Theta [Degrees]; Events", 100., 0., 180);
    TH1D* 	hNue_Phi = new TH1D("Nue_Phi_All","Nue_Phi_All; Phi [Degrees]; Events", 10., -180., 180);

    // Electron all
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

	TTreeReaderValue<int> PDGRV(myReader, "PDG"); 					// Load in the PDG variable
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

        hNue_E_vs_Theta->Fill(*NueEnergyRV, *NueThetaRV );
        hNue_E_vs_Phi->Fill(*NueEnergyRV, *NuePhiRV);


		hNue_Energy->Fill(*NueEnergyRV); 
        hNue_Theta->Fill(*NueThetaRV); 
        hNue_Phi->Fill(*NuePhiRV); 
        

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
    auto cNue_E_vs_Theta = new TCanvas();
    hNue_E_vs_Theta->Draw();
    hNue_E_vs_Theta->SetOption("COLZ,TEXT00");
    cNue_E_vs_Theta->Print("plots/Nue_Energy_vs_Theta.png");
    cNue_E_vs_Theta->Close();
    
    auto cNue_E_vs_Phi = new TCanvas();
    hNue_E_vs_Phi->Draw();
    hNue_E_vs_Phi->SetOption("COLZ,TEXT00");
    cNue_E_vs_Phi->Print("plots/Nue_Evergy_vs_Phi.png");
    
    auto cNue_Energy = new TCanvas();
    hNue_Energy->Draw();
    hNue_Energy->SetOption("HIST,TEXT00");
    cNue_Energy->Print("plots/Nue_Energy.png");

    auto cNue_Theta = new TCanvas();
    hNue_Theta->Draw();
    hNue_Theta->SetOption("HIST,TEXT00");
    cNue_Theta->Print("plots/Nue_Theta.png");

    auto cNue_Phi = new TCanvas();
    hNue_Phi->Draw();
    hNue_Phi->SetOption("HIST,TEXT00");
    cNue_Phi->Print("plots/Nue_Phi.png");

    
    // Electron
    auto celectron_E_vs_Theta = new TCanvas();
    helectron_E_vs_Theta->Draw();
    helectron_E_vs_Theta->SetOption("COLZ,TEXT00");
    celectron_E_vs_Theta->Print("plots/El_Energy_vs_Theta.png");
    
    auto celectron_E_vs_Phi = new TCanvas();
    helectron_E_vs_Phi->Draw();
    helectron_E_vs_Phi->SetOption("COLZ,TEXT00");
    celectron_E_vs_Phi->Print("plots/El_Energy_vs_Phi.png");
    
    auto celectron_Energy = new TCanvas();
    helectron_Energy->Draw();
    helectron_Energy->SetOption("HIST,TEXT00");
    celectron_Energy->Print("plots/El_Energy.png");

    //auto celectron_Theta = new TCanvas();
    // helectron_Theta->Draw();
    // helectron_Theta->SetOption("HIST,TEXT00");
    // celectron_Theta->Print("plots/El_Theta.png");

    //Plot1DHist(helectron_Theta, "plots/El_Theta.png" ); 
    //Plot1DHist(helectron_Phi, "plots/El_Phi.png" );

    //auto celectron_Phi = new TCanvas();
    //helectron_Phi->Draw();
    //helectron_Phi->SetOption("HIST,TEXT00");
    //celectron_Phi->Print("plots/El_Phi.png");
    

    MyFile->Write();
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //                                                     END
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

} // END MAIN FUNCTION


