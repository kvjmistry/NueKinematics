////////////////////////////////////////////////////////////////////////
// Class:       NueKinematics
// Plugin Type: analyzer (art v2_05_01)
// File:        NueKinematics_module.cc
//
// This module intends to loop over a data set and get the nue and el 
// kinematics: energy, phi and theta. 
//
// Generated at Mon Jul 23 07:00:45 2018 by Krishan Mistry using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

// Default art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"

// LARSOFT includes
#include "lardataobj/Simulation/SimChannel.h" 
#include "lardataobj/Simulation/SimPhotons.h" 
#include "lardataobj/Simulation/AuxDetSimChannel.h" 
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h" 
#include "larsim/Simulation/LArG4Parameters.h"                          
#include "larcore/Geometry/Geometry.h" 
#include "lardataobj/RecoBase/Hit.h" 
#include "nutools/ParticleNavigation/ParticleList.h" 
#include "nutools/ParticleNavigation/EmEveIdCalculator.h" 
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" 
//#include "nusimdata/SimulationBase/MCTrajectory.h" // for MCtrajectroy
//#include "lardata/DetectorInfo/DetectorPropertiesStandard.h"          // Commented out as causing a error when building//#include "larreco/Calorimetry/CalorimetryAlg.h" 			               
//#include "larreco/RecoAlg/ClusterRecoUtil/ClusterParamsAlg.h"		      // for cluster params alg function 
//#include "larevt/CalibrationDBI/Interface/ElectronLifetimeService.h"    // For dEdx conversion
//#include "larevt/CalibrationDBI/Interface/ElectronLifetimeProvider.h"   // For dEdx conversion
//#include "larreco/Calorimetry/CalorimetryAlg.h" 


// ROOT includes
#include "TH1.h" 
#include "TH2.h" 
#include "TH3.h" 
#include "TTree.h" 
#include "TDatabasePDG.h" 
#include "TParticlePDG.h" 
#include "TCanvas.h" 
#include "TVectorT.h" 
#include "TMatrixT.h" 
#include "TMatrixDUtils.h" 
#include "TGraph.h" 
#include "TF1.h" 
#include "TMath.h"

// C++ Includes
#include <cmath>
#include <map>
#include <iostream>
#include <algorithm>

#define PI 3.14159265


class NueKinematics;
class NueKinematics : public art::EDAnalyzer {
public:
  explicit NueKinematics(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NueKinematics(NueKinematics const &) = delete;
  NueKinematics(NueKinematics &&) = delete;
  NueKinematics & operator = (NueKinematics const &) = delete;
  NueKinematics & operator = (NueKinematics &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  
  // Optional Member functions
  double Calc_Theta(double Pz, double P);
  double Calc_Phi  (double Px, double Py, double P);


private:

  // Declare member data here.
  std::string fSimulationProducerLabel;			// The name of the producer that tracked simulated particles through the detector

  // Histograms
  // Nue all
  TH1D* 	hNue_Energy_All;	 
  TH1D* 	hNue_Theta_All;	  
  TH1D* 	hNue_Phi_All;	

  TH2D* 	hNue_E_vs_Theta_All;	 
  TH2D* 	hNue_E_vs_Phi_All;   

  // Nue
  TH1D* 	hNue_Energy;	
	TH1D*	  hNue_Theta;		 
	TH1D* 	hNue_Phi;			  
  
  // Nue bar
  TH1D* 	hNue_bar_Energy;	
	TH1D*	  hNue_bar_Theta;		 
	TH1D* 	hNue_bar_Phi;			 

  // Electron all
  TH1D* 	helectron_Energy_All;	 
  TH1D* 	helectron_Theta_All;	  
  TH1D* 	helectron_Phi_All;	    

  TH2D* 	helectron_E_vs_Theta_All;	 
  TH2D* 	helectron_E_vs_Phi_All;  
 
  // EMinus
  TH1D* 	heminus_Energy;	
	TH1D*	  heminus_Theta;		 
	TH1D* 	heminus_Phi;	

    // EPlus
  TH1D* 	heplus_Energy;	
	TH1D*	  heplus_Theta;		 
	TH1D* 	heplus_Phi;	

};


NueKinematics::NueKinematics(fhicl::ParameterSet const & p) : EDAnalyzer(p) {
   // More initializers here.
  fSimulationProducerLabel = p.get< std::string >("SimulationLabel");
}


double NueKinematics::Calc_Theta(double Pz, double P){

  double Z_dir{ Pz / P }; // Z direction
  double Theta{ std::acos(Z_dir) * (180. / PI) }; // Theta

  return Theta; 
}

double NueKinematics::Calc_Phi(double Px, double Py, double P){

  double Y_dir {Py / P}; // Y direction 
  double X_dir {Px / P}; // X direction
  
  double Phi{std::atan2(Y_dir, X_dir) * (180. / PI)}; // Phi

  return Phi; 
}


void NueKinematics::beginJob() {
  // Implementation of optional member function here.

  // Access ART's TFileService, which will handle histograms/trees/etc.
	art::ServiceHandle<art::TFileService> tfs;
	
  // Make directories
  art::TFileDirectory Nue_All_dir = tfs->mkdir( "Nue_All_dir" );
  art::TFileDirectory Nue_dir     = tfs->mkdir( "Nue_dir"     );
  art::TFileDirectory Nue_bar_dir = tfs->mkdir( "Nue_bar_dir" );

  art::TFileDirectory el_All_dir = tfs->mkdir( "Nue_All_dir" );
  art::TFileDirectory eMinus_dir = tfs->mkdir( "eMinus_dir"  );
  art::TFileDirectory ePlus_dir  = tfs->mkdir( "ePlus_dir"   );

  // Histograms 

  // Nue All
  hNue_E_vs_Theta_All = Nue_All_dir.make<TH2D>("hNue_E_vs_Theta_All","hNue_E_vs_Theta_All; Energy [GeV]; Theta [degrees]",30, 0., 20 , 30, -180., 180);
  hNue_E_vs_Phi_All   = Nue_All_dir.make<TH2D>("hNue_E_vs_Phi_All","hNue_E_vs_Phi_All; Energy [GeV]; Phi [degrees]",30, 0., 20 , 30, -180., 180);
  
  hNue_Energy_All = Nue_All_dir.make<TH1D>("Nue_Energy_All","Nue_Energy_All; E [GeV]; Events",30, 0., 20);
  hNue_Theta_All  = Nue_All_dir.make<TH1D>("Nue_Theta_All","Nue_Theta_All; Theta [Degrees]; Events", 30, -180., 180);
  hNue_Phi_All    = Nue_All_dir.make<TH1D>("Nue_Phi_All","Nue_Phi_All; Theta [Degrees]; Events", 30, -180., 180);

  // Nue
  hNue_Energy = Nue_dir.make<TH1D>("Nue_Energy","Nue_Energy; E [GeV]; Events",30, 0., 20);
  hNue_Theta  = Nue_dir.make<TH1D>("Nue_Theta","Nue_Theta; Theta [Degrees]; Events", 30, -180., 180);
  hNue_Phi    = Nue_dir.make<TH1D>("Nue_Phi","Nue_Phi; Theta [Degrees]; Events", 30, -180., 180);

  // Nue bar
  hNue_bar_Energy = Nue_bar_dir.make<TH1D>("Nue_bar_Energy","Nue_bar_Energy; E [GeV]; Events",30, 0., 20);
  hNue_bar_Theta  = Nue_bar_dir.make<TH1D>("Nue_bar_Theta","Nue_bar_Theta; Theta [Degrees]; Events", 30, -180., 180);
  hNue_bar_Phi    = Nue_bar_dir.make<TH1D>("Nue_bar_Phi","Nue_bar_Phi; Theta [Degrees]; Events", 30, -180., 180);

  // Electron all
  helectron_E_vs_Theta_All = el_All_dir.make<TH2D>("helectron_E_vs_Theta_All","helectron_E_vs_Theta_All; Energy [GeV]; Theta [degrees]",30, 0., 20 , 30, -180., 180);
  helectron_E_vs_Phi_All   = el_All_dir.make<TH2D>("helectron_E_vs_Phi_All","helectron_E_vs_Phi_All; Energy [GeV]; Phi [degrees]",30, 0., 20 , 30, -180., 180);
  
  helectron_Energy_All = el_All_dir.make<TH1D>("electron_Energy_All","electron_Energy_All; E [GeV]; Events",30, 0., 20);
  helectron_Theta_All  = el_All_dir.make<TH1D>("electron_Theta_All","electron_Theta_All; Theta [Degrees]; Events", 30, -180., 180);
  helectron_Phi_All    = el_All_dir.make<TH1D>("electron_Phi_All","electron_Phi_All; Theta [Degrees]; Events", 30, -180., 180);

  // Eminus
  heminus_Energy = eMinus_dir.make<TH1D>("eminus_Energy","eminus_Energy; E [GeV]; Events",30, 0., 20);
  heminus_Theta  = eMinus_dir.make<TH1D>("eminus_Theta","eminus_Theta; Theta [Degrees]; Events", 30, -180., 180);
  heminus_Phi    = eMinus_dir.make<TH1D>("eminus_Phi","eminus_Phi; Theta [Degrees]; Events", 30, -180., 180);

  // EPlus
  heplus_Energy = ePlus_dir.make<TH1D>("eplus_Energy","eplus_Energy; E [GeV]; Events",30, 0., 20);
  heplus_Theta  = ePlus_dir.make<TH1D>("eplus_Theta","eplus_Theta; Theta [Degrees]; Events", 30, -180., 180);
  heplus_Phi    = ePlus_dir.make<TH1D>("eplus_Phi","eplus_Phi; Theta [Degrees]; Events", 30, -180., 180);

}

void NueKinematics::analyze(art::Event const & e) {
  // Implementation of required member function here.
	
	// Create vector of particles involved in event
	art::Handle <std::vector <simb::MCParticle>> particleHandle;
	e.getByLabel(fSimulationProducerLabel, particleHandle);
	
  double ParticleE{ 0. };
  double Theta{ 0. }; 
  double Phi{ 0. }; 

  for(auto particlePtr = particleHandle->begin(); particlePtr != particleHandle->end(); ++particlePtr) { // Loop over PDG particles in the event
	  const simb::MCParticle& particle = (*particlePtr); // De-reference the particle pointer

    double Theta =  Calc_Theta{ particle.Pz(), particle.P()  };                          // Calculate Theta
    //double Phi   =  Calc_Phi{   particle.Momentum(0).Px(), particle.Momentum(0).Py(),  particle.Momentum(0).P()  }; // Calculate Phi
    

    ParticleE = particle.Momentum(0).E() ;

    // BEGIN SELECTING PARTICLE BLOCK
		if (particle.PdgCode() == 12){ // nue in the event
      
      // Fill a histogram with the pdg code, energy, theta, phi
      hNue_Energy_All ->Fill( ParticleE );
      hNue_Energy     ->Fill( ParticleE );

      
    } 
    else if (particle.PdgCode() == -12){ // nue bar in the event
      
      // Fill a histogram with the pdg code, energy, theta, phi
      hNue_Energy_All ->Fill( ParticleE );
      hNue_bar_Energy ->Fill( ParticleE );

    } 
    else if (particle.PdgCode() == 11){ // e- in the event
      // Fill a histogram with the pdg code, energy, theta, phi
      helectron_Energy_All ->Fill( ParticleE );
      heminus_Energy   ->Fill( ParticleE );

    } 
    else if (particle.PdgCode() == -11){ // e+ in the event
      // Fill a histogram with the pdg code, energy, theta, phi
      helectron_Energy_All ->Fill( ParticleE );
      heplus_Energy        ->Fill( ParticleE );

    }// END IF CONDITION BLOCK


	} // End loop over PDG particles in the event

}


void NueKinematics::endJob() {
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(NueKinematics)