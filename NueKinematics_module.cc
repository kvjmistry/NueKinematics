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
#include "lardataobj/Simulation/SimPhotons.h" 
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h" 
#include "larsim/Simulation/LArG4Parameters.h"                          
#include "larcore/Geometry/Geometry.h" 
#include "lardataobj/RecoBase/Hit.h" 
#include "nutools/ParticleNavigation/ParticleList.h" 
#include "nutools/ParticleNavigation/EmEveIdCalculator.h" 
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" 
#include "nusimdata/SimulationBase/MCTruth.h"


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
  std::string fMCTruth_label; 

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

  int NC_event{0};

  // Add Ttree variables 
  TTree *DataTree;
  int run, subrun, evt;
  double NueEnergy,       NueTheta,      NuePhi;
  double ElectronEnergy, ElectronTheta, ElectronPhi;
  int PDG; 


};


NueKinematics::NueKinematics(fhicl::ParameterSet const & p) : EDAnalyzer(p) {
   // More initializers here.
  fSimulationProducerLabel = p.get< std::string >("GeantModule");
  fMCTruth_label = p.get<std::string>("MCTruthProduct", "generator");
  
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


  // Create the TTree and add relavent branches
  DataTree = tfs->make<TTree>("EventTree","EventTree");
	
  // Make directories
  art::TFileDirectory Nue_All_dir = tfs->mkdir( "Nue_All_dir" );
  art::TFileDirectory Nue_dir     = tfs->mkdir( "Nue_dir"     );
  art::TFileDirectory Nue_bar_dir = tfs->mkdir( "Nue_bar_dir" );

  art::TFileDirectory el_All_dir = tfs->mkdir( "el_All_dir" );
  art::TFileDirectory eMinus_dir = tfs->mkdir( "eMinus_dir"  );
  art::TFileDirectory ePlus_dir  = tfs->mkdir( "ePlus_dir"   );

  // Histograms 

  // Nue All
  hNue_E_vs_Theta_All = Nue_All_dir.make<TH2D>("Nue_E_vs_Theta_All","Nue_E_vs_Theta_All; Energy [GeV]; Theta [degrees]",15., 0., 7. , 10., 0., 180);
  hNue_E_vs_Phi_All   = Nue_All_dir.make<TH2D>("Nue_E_vs_Phi_All","Nue_E_vs_Phi_All; Energy [GeV]; Phi [degrees]",10., 0., 10. , 10., -180., 180);
  hNue_E_vs_Theta_All ->SetOption("COLZ,TEXT");
  hNue_E_vs_Phi_All   ->SetOption("COLZ,TEXT");
  
  hNue_Energy_All = Nue_All_dir.make<TH1D>("Nue_Energy_All","Nue_Energy_All; E [GeV]; Events",50., 0., 5);
  hNue_Theta_All  = Nue_All_dir.make<TH1D>("Nue_Theta_All","Nue_Theta_All; Theta [Degrees]; Events", 100., 0., 180);
  hNue_Phi_All    = Nue_All_dir.make<TH1D>("Nue_Phi_All","Nue_Phi_All; Phi [Degrees]; Events", 10., -180., 180);
  hNue_Energy_All ->SetOption("HIST,TEXT00");
  hNue_Theta_All  ->SetOption("HIST,TEXT00");
  hNue_Phi_All    ->SetOption("HIST,TEXT00");

  // Nue
  hNue_Energy = Nue_dir.make<TH1D>("Nue_Energy","Nue_Energy; E [GeV]; Events",50., 0., 5.);
  hNue_Theta  = Nue_dir.make<TH1D>("Nue_Theta","Nue_Theta; Theta [Degrees]; Events", 200., 0., 180);
  hNue_Phi    = Nue_dir.make<TH1D>("Nue_Phi","Nue_Phi; Phi [Degrees]; Events", 10., -180., 180);
  hNue_Energy ->SetOption("HIST,TEXT00");
  hNue_Theta  ->SetOption("HIST,TEXT00");
  hNue_Phi    ->SetOption("HIST,TEXT00");

  // Nue bar
  hNue_bar_Energy = Nue_bar_dir.make<TH1D>("Nue_bar_Energy","Nue_bar_Energy; E [GeV]; Events",50., 0., 5.);
  hNue_bar_Theta  = Nue_bar_dir.make<TH1D>("Nue_bar_Theta","Nue_bar_Theta; Theta [Degrees]; Events", 10., 0., 180);
  hNue_bar_Phi    = Nue_bar_dir.make<TH1D>("Nue_bar_Phi","Nue_bar_Phi; Phi [Degrees]; Events", 10., -180., 180);
  hNue_bar_Energy ->SetOption("HIST,TEXT00");
  hNue_bar_Theta  ->SetOption("HIST,TEXT00");
  hNue_bar_Phi    ->SetOption("HIST,TEXT00");  

  // Electron all
  helectron_E_vs_Theta_All = el_All_dir.make<TH2D>("electron_E_vs_Theta_All","electron_E_vs_Theta_All; Energy [GeV]; Theta [degrees]",20., 0., 10. , 10., 0., 180);
  helectron_E_vs_Phi_All   = el_All_dir.make<TH2D>("electron_E_vs_Phi_All","electron_E_vs_Phi_All; Energy [GeV]; Phi [degrees]",20., 0., 10. , 10., -180., 180);
  helectron_E_vs_Theta_All->SetOption("COLZ,TEXT");
  helectron_E_vs_Phi_All->SetOption("COLZ,TEXT");
  
  helectron_Energy_All = el_All_dir.make<TH1D>("electron_Energy_All","electron_Energy_All; E [GeV]; Events",50., 0., 5.);
  helectron_Theta_All  = el_All_dir.make<TH1D>("electron_Theta_All","electron_Theta_All; Theta [Degrees]; Events", 10., 0., 180);
  helectron_Phi_All    = el_All_dir.make<TH1D>("electron_Phi_All","electron_Phi_All; Phi [Degrees]; Events", 10., -180., 180);
  helectron_Energy_All  ->SetOption("HIST,TEXT00");
  helectron_Theta_All   ->SetOption("HIST,TEXT00");
  helectron_Phi_All     ->SetOption("HIST,TEXT00");

  // Eminus
  heminus_Energy = eMinus_dir.make<TH1D>("eminus_Energy","eminus_Energy; E [GeV]; Events",50., 0., 5.);
  heminus_Theta  = eMinus_dir.make<TH1D>("eminus_Theta","eminus_Theta; Theta [Degrees]; Events", 10., 0., 180);
  heminus_Phi    = eMinus_dir.make<TH1D>("eminus_Phi","eminus_Phi; Phi [Degrees]; Events", 10., -180., 180);
  heminus_Energy ->SetOption("HIST,TEXT00");
  heminus_Theta ->SetOption("HIST,TEXT00");
  heminus_Phi   ->SetOption("HIST,TEXT00");

  // EPlus
  heplus_Energy = ePlus_dir.make<TH1D>("eplus_Energy","eplus_Energy; E [GeV]; Events",50., 0., 5.);
  heplus_Theta  = ePlus_dir.make<TH1D>("eplus_Theta","eplus_Theta; Theta [Degrees]; Events", 10., 0., 180);
  heplus_Phi    = ePlus_dir.make<TH1D>("eplus_Phi","eplus_Phi; Phi [Degrees]; Events", 10., -180., 180);
  heplus_Energy->SetOption("HIST,TEXT00");
  heplus_Theta->SetOption("HIST,TEXT00");
  heplus_Phi->SetOption("HIST,TEXT00");

  // Add Tree Information
	DataTree->Branch("run",   &run);
	DataTree->Branch("subrun",&subrun);
  DataTree->Branch("event", &evt);
  DataTree->Branch("PDG",    &PDG);

  DataTree->Branch("NueEnergy", &NueEnergy);
  DataTree->Branch("NueTheta",  &NueTheta);
  DataTree->Branch("NuePhi",    &NueEPhi);
  
  DataTree->Branch("ElectronEnergy", &ElectronEnergy);
  DataTree->Branch("ElectronTheta",  &ElectronTheta;
  DataTree->Branch("ElectronPhi",    &ElectronPhi);

}

void NueKinematics::analyze(art::Event const & e) {
  // Implementation of required member function here.

  std::cout << "Running over Event:\t" << e.event() << std::endl;

  // Determine event ID
  run = e.id().run();
  subrun = e.id().subRun();
  evt = e.id().event();

	
  double ParticleE{ 0. };
  double Theta{ 0. }; 
  double Phi{ 0. }; 
  NC_event = 0; // reset NC counter

  // create vector of GENIE MCTruth involved in event
	art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  
  if (e.getByLabel("generator",mctruthListHandle)) art::fill_ptr_vector(mclist, mctruthListHandle);

  int iList{ 0 }; // 1 nu int per spill

  for (int p = 0; p < mclist[iList]->NParticles(); p++) { // Loop over GENIE MCTruth Particles

    //if (mclist[iList]->GetNeutrino().CCNC() == 1) NC_event ++; // The event was NC and we dont want to include the additional scattered neutrino
    if (mclist[iList]->GetNeutrino().CCNC() == 1) continue; // The event was NC then skip the event


    simb::MCParticle particle{mclist[iList]->GetParticle(p)}; // Get a MC Particle

    Theta =  Calc_Theta( particle.Pz(), particle.P()  );                 // Calculate Theta
    Phi   =  Calc_Phi(   particle.Px(), particle.Py(),  particle.P()  ); // Calculate Phi
    ParticleE = particle.E() ;                                           // Energy

    if ( NC_event >= 2 ) continue;// Skip fill for the second neutrino if NC

    if (mclist[iList]->Origin() == simb::kBeamNeutrino){ // Require a beam neutrino
      //DEBUG: std::cout << "PDG:\t"<< particle.PdgCode()<< std::endl;
      
      // BEGIN SELECTING PARTICLE BLOCK
      if (particle.PdgCode() == 12){ // nue in the event
        
        // Fill a histogram with the pdg code, energy, theta, phi
        hNue_Energy_All             ->Fill( ParticleE );
        hNue_Theta_All              ->Fill(Theta); 
        hNue_Phi_All                ->Fill(Phi); 
        hNue_E_vs_Theta_All         ->Fill(ParticleE, Theta);
        hNue_E_vs_Phi_All           ->Fill(ParticleE, Phi);
        
        hNue_Energy ->Fill( ParticleE );
        hNue_Theta  ->Fill(Theta); 
        hNue_Phi    ->Fill(Phi);

        // Set laues for filling TTree
        NueEnergy =  ParticleE;
        NueTheta  =  Theta; 
        NuePhi    =  Phi;
        PDG       =  particle.PdgCode();

        
      } 
      else if (particle.PdgCode() == -12){ // nue bar in the event

        // Fill a histogram with the pdg code, energy, theta, phi
        hNue_Energy_All       ->Fill( ParticleE );
        hNue_Theta_All        ->Fill( Theta ); 
        hNue_Phi_All          ->Fill( Phi ); 
        hNue_E_vs_Theta_All   ->Fill( ParticleE, Theta );
        hNue_E_vs_Phi_All     ->Fill( ParticleE, Phi );

        hNue_bar_Energy ->Fill( ParticleE );
        hNue_bar_Theta  ->Fill( Theta ); 
        hNue_bar_Phi    ->Fill( Phi );

        // Set laues for filling TTree
        NueEnergy =  ParticleE;
        NueTheta  =  Theta; 
        NuePhi    =  Phi;
        PDG       =  particle.PdgCode();

      } 
      else if (particle.PdgCode() == 11){ // e- in the event

        // Fill a histogram with the pdg code, energy, theta, phi
        helectron_Energy_All     ->Fill( ParticleE );
        helectron_Theta_All      ->Fill( Theta );
        helectron_Phi_All        ->Fill( Phi );
        helectron_E_vs_Theta_All ->Fill( ParticleE, Theta );
        helectron_E_vs_Phi_All   ->Fill( ParticleE, Phi );

        heminus_Energy   ->Fill( ParticleE );
        heminus_Theta    ->Fill( Theta );
        heminus_Phi      ->Fill( Phi );

        // Set laues for filling TTree
        ElectronEnergy =  ParticleE;
        ElectronTheta  =  Theta; 
        ElectronPhi    =  Phi;
        PDG            =  particle.PdgCode();

      } 
      else if (particle.PdgCode() == -11){ // e+ in the event

        // Fill a histogram with the pdg code, energy, theta, phi
        helectron_Energy_All     ->Fill( ParticleE );
        helectron_Theta_All      ->Fill( Theta );
        helectron_Phi_All        ->Fill( Phi );
        helectron_E_vs_Theta_All ->Fill( ParticleE, Theta );
        helectron_E_vs_Phi_All   ->Fill( ParticleE, Phi );

        heplus_Energy   ->Fill( ParticleE );
        heplus_Theta    ->Fill( Theta );
        heplus_Phi      ->Fill( Phi );

        // Set laues for filling TTree
        ElectronEnergy =  ParticleE;
        ElectronTheta  =  Theta; 
        ElectronPhi    =  Phi;
        PDG            =  particle.PdgCode();

      }// END IF CONDITION BLOCK  
    
    }// END if a beam neutrino

  }// END loop over mclist

  DataTree->Fill();

}



void NueKinematics::endJob() {
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(NueKinematics)
