# Network models of spatiotemporal On-Off dynamics


## What is contained here

This repository contains the code for the model simulations presented in :  

 **Yan-Liang Shi,  Nicholas A. Steinmetz,  Tirin Moore,  Kwabena Boahen,  and Tatiana A. Engel, Cortical state dynamics and selective attention define the spatial pattern of correlated variability in neocortex. Nat Commun (2021)**.

**Link:** [https://doi.org/10.1101/2020.09.02.279893](https://doi.org/10.1101/2020.09.02.279893).





* Simulation code for a dynamical system network model of population On-Off dynamics 
* Simulation code for a binary-unit reduced network model of population On-Off dynamics  
* Computation of noise correlations based on simulations


## Code description 


### 1. Dynamical system network model

* Directory: /Dynamical_system_model
* Simulation_dynamical_system_model.m is the simulation code for the model, with input from random initial conditions in FHN_init.mat. After each simulation, On/Off population phases of 121 units in attentional and control condition are generated.
* After running Simulation_dynamical_system_model.m, noise correlations are computed based on default uniform On/Off firing rates (r_{on}=130 Hz, r_{off}=50 Hz).
* To compute noise correlations with On/Off firing rates sampled from data (fitted with a 2-state HMM across 31 recording sessions), we can use NC_with_input_firing_rate.m, with input 1) On/Off population phases generated by Simulation_dynamical_system_model.m, 2) distribution of On/Off firing rates sampled from data (distribution_firing_rate.mat). To obtain sufficient samples of On/Off firing rates, the sampling procedure is repeated 10 times and the average noise correlation is computed.



### 2. Binary-unit network model

* Directory: /Binary_units_model
* Simulation_binary_units_model.m is the simulation code for the model, with input from random initial conditions in SPIN_init.mat. After each simulation, On/Off population phases of 121 units in attentional and control condition are generated.
* After running Simulation_binary_units_model.m, noise correlations are computed based on default uniform On/Off firing rates (r_{on}=125 Hz, r_{off}=25 Hz).

