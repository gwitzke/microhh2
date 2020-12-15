# CONSTRAIN Cold Air Outbreak case

## Case Description
This case is a Large Eddy Model intercomparison study of the CONSTRAIN cold air outbreak case. 
The CONSTRAIN CAO is a cold weather variant of Stratocumlus to Cumulus Transitions (SCTs) 
in which a cold stratocumumlus cloud deck is advected over an increasingly warmer sea surface temperature
(SST) resulting in a deepening of the boundry layer and the subsequent breaking up of the stratocumulus 
deck into a open-cell structured cumulus cloud field.

Research Article: "De Roode et al.(2019): Turbulent Transport in the Gray Zone: 
A Large Eddy Model Intercomparison Study of the CONSTRAIN Cold Air Outbreak Case"

## Running a Simulation with Condor
Copy necessary files into directory in which you wish the data to be saved to.
Necessary files:
constrain.ini
constrain_input.py
constrain_setup_forcing.nc
constrain.sh
constrain.sub
microhh

Then its a simple command of : condor_submit constrain.sub
