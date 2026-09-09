# Deriving novel signatures

In this folder, we include the code for deriving CIN signatures that can be mapped to individual events.

This folder needs to be run in the following order:

1.  Run the `00_WORKFLOW_Derive_GoldSignatures_newCHP.R` script for quantifying gold-standard exposures used for signature derivation

2.  Run the `01_WORKFLOW_Derive_NewSignatures.R` script for deriving the novel signatures used for mapping each individual event to specific CIN type


We suggest running the MC simulations in a HPC cluster. The required dependencies to run the scripts are stored in env_mcsimulations.yml. 
Activate as follows: 

	conda update --all
	conda env create -f env_mcsimulations.yml
	conda activate focas_mcsims
