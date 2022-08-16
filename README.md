# MUFT-model
Metabolism estimation in Unsteady Flow conditions and Transient storage zones. 

MUFT_model general:

MUFT_model_general includes scripts for general model application in rivers. The general solute model presented here includes equations for the ADZ formulation (including transient storage), but the equations can be modified for a well-mixed reach (without transient storage) or ADV (advetion dominated solute transport) formulation if required (equations provided in the manuscript).
The general model scripts use artificial input data that were generated based on the availablabe data in the River Otra.
Please see input files (.csv) for the input data requirements and format your input data accordingly. Please change parameter values and other information (e.g., number of reaches, model time-steps, etc.) in the script (at the top) according to your case study.
A sample model application is shown here for 2 river reaches using input data of 3 days.
To accomodate the model to look for values in the past (for transient storage time-lag), the model simulates data for a period of n-1 days at each reach, where n is the total number of simulation days for the preceding reach.

MUFT_Otra:

MUFT_Otra contains model scripts and input files used for flow and solute routing in the River Otra. 
The model for the River Otra uses a different approach for parameter estimation compared to the general model and thus, the MUFT model equations are modified according to the data availability in the River Otra (explained in the manuscript).
MUFT_flow_param.py optimised the *Fr* parameter, which is then used for routing along the river network in the MUFT_flow_routing.py.
The estimates of *P<sub>GPPmax</sub>*, *k<sub>PAR</sub>* and *R<sub>ER</sub>* parameters in the two-station models were used to constrain inverse model parameters.
The scripts provided in the MUFT_model_general can be directly used for your case study. 