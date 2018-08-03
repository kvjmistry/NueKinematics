A module that gets the MC True Kinematics of Nue and electrons from a Genie MC Sample.

This repository consists of a larsoft module NueKinetmatics_module.cc which can be exectued by running  

`lar -c NueKinematics.fcl`

The module will create a .root file containing some histograms of the key variables of interest and will also spit most of the interesting variables into a TTree. 

The AnalyseNueKinematics.C script is a c script written for analysing the output n-tuples created from the larsoft module. This scripts allows for a quick analysis of the variables rather than having to re-run the larsoft module which is time consuming. 
