# Yeast lipid Metabolism model
Model of the lipid metabolism in the yeast *Saccharomyces cerevisiae*    
Written by Vera Schützhold, Jens Hahn - 2015    
Written in Python 2.7

## Files 
### Model files
* components: stores multiple component classes for initialisation of lipid objects. Called by **reactions**.
	* lipids (phospholipids, TAG, CL)
	* fatty acids
	* sterol
	* sterylester
	* sphingolipid

* reactions: reaction class, performs all reactions in one time step. Called by **model**
    * run_step (permute order of reactions)
    * glycerol 3-p synthesis
    * inositol synthesis
    * ceramide synthesis
    * acetyl-coa synthase
    * acyl synthase
    * PA synthesis
    * lyso PA synthase
    * PA synthase
    * CDP-DG synthase
    * TAG synthesis
    * DAG synthase
    * TAG synthase
    * TAG lipase
    * DAG kinase
    * PS synthase
    * PI synthase
    * PE synthase
    * PC synthase
    * CL synthase
    * ergosterol synthase
    * sterylester synthase
    * sphingolipid synthase
    * transport reaction

* model: model class
    * __init__ (initialise parameters, lists, precursors, molecules and simulation)
    * run (start simulation, return 3 Python dicts - saturation dist., membrane lengths and composition)
    * start (generate initial cell state)
    * cell_cycle (cell cycle tracker)
    * numbers (counts and saves membrane lengths)
    * membrane_compositions (saves membrane compositions)
    * saturation_counter (saves FA saturations)

## Examples
* comparison: comparison between stochastic and deterministic MM (supplementary Figure S4)

## Simulation
* simulate_print: perform 1000 simulations and save results in csv files
* test_ergosterol: disable ergosterol synthase, perform 1000 simulations and save results in csv files
* test_inositol: increase inositol addition, perform 1000 simulations and save results in csv files

## Sensitivity analysis
* rates_SA: change rate of reaction, perform 1000 simulations and save results in csv files
This script has to be called with reaction name as argument:
e.g.: ```python rates_SA.py glycerol_3_p_synthesis```

Names of reactions and change of rates are shown here, initial values shown in brackets:    
glycerol_3_p_synthesis:  + 1 (8)    
inositol_synthesis: + 1 (5)    
ceramide_synthesis: + + 1 (2)    
acetyl_coa_synthase: + 65 (650)    
acyl_synthase: + 45 (450)    
PA_synthesis: + 2 (17)    
CDP_DG_synthase: + 2 (20)     
TAG_synthesis: + 3 (34)    
TAG_lipase: + 2 (23)    
DAG_kinase: + 4 (40)    
PS_synthase: + 2 (18)    
PI_synthase: + 1 (6)    
PE_synthase: + 2 (12)    
PC_synthase: + 1 (5)    
CL_synthase: + 1 (2)    
ergosterol_synthase: + 3 (25)    
sterylester_synthase: + 3 (27)    
sphingolipid_synthase: + 1 (2)
