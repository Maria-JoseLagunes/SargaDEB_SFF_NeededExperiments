%%% Writing the protocol for the identification procecedure
1) **Define reference parameters + model structure**
   Functions : pars_init_Sargassum_fluitans + get_mECENV
2) **Simulate existent data**
   Functions :   Identifiability/get_prdData.m\
   To create file with the predicted data based on reference parameters and see the plots of these data
3) **Create pseudo experiments**
   Functions : Identifiability/ExperimentalCreation/creationExperiments.m\
   *Note* I made a loop for experiments to run each time the code is ran. If an new experiment is created, one has to do it manually and modified manually the mydata, predict and list ofdatasets to include this new dataset.
  - To modify pseudo experiments forcing conditions, time of experience or observable initial conditions go to the specific binder file. when running a whole simulation modifiy the init.m file and if modifying N uptake, modify the run_Nuptake.m.
  - To modify the measurement times, in the creationExperiments.m modify which data points are going to be saved for each experiment using the *setDays* variable for each pseudo experiment. 

5) **Run the identification procedure**
   Function : Identifiability/run_ID_analysis.m\
   First in *all_which_data* variable, select those datasets to be used in the estimation procedure. A loop can be done by including different combinations of datasets\
   Exemple:  all_which_data = {
    [1 3 4 5 ]
    [1 3 4 5 6]};\
Results will be saved in SSE_Analysis binder as SSE_results_"nummber of parameter estimated"_datasets used
 
   Calls :
   - IdentificationProcedure_NE.m which calls the estimation, calculates the SSE for each parameters and saves the results in a structure. Here we select the variation factor used
        - Calls:
             -  Estimate_parameters : changes weights of datasets not selected and does the estimation according to the variation factor selected.
