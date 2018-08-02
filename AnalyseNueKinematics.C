// Program that loads in the ttree created from the NueKinematics Module and makes some plots as does the Larsoft Module. 

#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TEfficiency.h"

#include <iostream>
#include <fstream>

// Main Function
void AnalyseNueKinematics() {
    		
	// Open the NueKinematics File
	TFile *fileIN = TFile::Open("NueKinematics.root");
	if (fileIN == 0) {
	  // if we cannot open the file, print an error message and return immediatly
	  printf("Error: cannot open NueKinematics.root!\n");
	  return;
	}

    double Total_events = 0;
	double Passed_Events = 0;

	// Create the tree reader and its data containers
	TTreeReader myReader("microboonewvtof/EventTree", fileIN);

	TTreeReaderValue<std::vector<int>> PDGRV(myReader, "PDG"); 					// Load in the PDG variable
	TTreeReaderValue<std::vector<double>> NueEnergyRV(myReader, "NueEnergy"); 	// Load in the NueEnergy variable

	// Loop over all entries of the TTree or TChain.
	while (myReader.Next()) {
		if ( (*PDGRV)[0] >= 1  ) Total_events++; 	

		if ( (*NueEnergyRV)[0] == 1   ) Passed_Events++; 

	}

    std::cout << "\nTest Output : " << Total_events << std::endl;
    std::cout << "\nTest Output2 : " << Passed_Events << std::endl;

}


