# MUFT-model
### Metabolism estimation in Unsteady Flow conditions and Transient storage zones (MUFT). 

The MUFT model presents a modelling approach to estimate whole-stream metabolism in hydropeaking rivers with/without significant transient storage. The model is developed by coupling an unsteady flow routing model (adapted from Sincock and Lees, 2002) with the two-station stream metabolism method (Odum, 1956). 
The solute model includes two formulations, ADV (advection; Beck and Young, 1975) and ADZ (aggregated dead-zone; Wallis et al., 1989), that can be selected based on the dominant solute transport mechanism in the river of interest.
The MUFT model can be applied using an accounting method and/or inverse modelling approach. 

## MUFT_model general:

MUFT_model_general includes scripts for general model application in rivers. 

The general solute model presented here includes equations for the ADZ formulation (including transient storage), but the equations can be modified for a well-mixed reach (without transient storage) or ADV (advetion dominated solute transport) formulation if required (equations provided in the manuscript).

<p>The general model scripts here use artificial input data that were generated based on the availablabe data in the River Otra, Norway.
A sample model application is shown here for 2 river reaches using input data of 3 days.
The scripts provided in the MUFT_model_general can be directly used for your case study.</p> 
<ul class="no-bullets">
  <li>See input files (.csv) for the input data requirements and format your input data accordingly.</li>
  <li>Change parameter values and other information (e.g., number of reaches, model time-steps, etc.) in the script (at the top) according to your case study.</li>
  <li>To accomodate the model to look for values in the past (for transient storage time-lag), the model simulates data for a period of n-1 days at each reach, where n is the total number of simulation days for the preceding reach.
</li>
</ul>


## MUFT_Otra:

MUFT_Otra contains model scripts and input files used for flow and solute modelling in the River Otra. 
The model for the River Otra uses a different approach for parameter estimation compared to the general model and thus, the MUFT model equations are modified according to the data availability in the River Otra (explained in the manuscript).

The flow routing model is adapted from the unsteady QUASAR model. The general flow routing model uses the same equations as Sincock and Lees, 2002 (3 parameter model). For the River Otra, we have modified the flow model equations (1 parameter model) to reduce the number of paraneters and to make use of the available data in the river (discussed in the manuscript).
MUFT_flow_param.py optimises the *Fr* parameter, which is then used for routing along the river network in the MUFT_flow_routing.py.

For metabolism estimation, Otra_twostation demonstrates the accounting approach and Otra_inverse demonstrates the inverse modelling approach.
In the case study presented here, the estimates of *P<sub>GPPmax</sub>*, *k<sub>PAR</sub>* and *R<sub>ER</sub>* parameters derived in the two-station models were used to constrain the inverse model parameters.

## <p>References:</p>
<ul class="no-bullets">
  <li>Beck, M., & Young, P. C. (1975). A dynamic model for do-bod relationships in a non-tidal stream. Water Research, 9 (9), 769–776.</li>
  <li>Odum, H. T. (1956). Primary production in flowing waters. Limnology and Oceanography, 1 (2), 102–117.
</li>
<li>Sincock, A. M., & Lees, M. (2002). Extension of the quasar river-water quality model to unsteady flow conditions. Water and Environment Journal, 16 (1), 12–17.
  <li>Wallis, S. G., Young, P., & Beven, K. (1989). Experimental investigation of the ag935 gregated dead zone model for longitudinal solute transport in stream channels. Proceedings of the Institution of Civil Engineers, 87 (1), 1-22.
</li>
</ul>