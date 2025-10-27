%% Estimation procedure function
% Estimate parameters based on selected data and new seed of reference
% parameters

function resultsStructure = estimate_parameters(varFactor, label, data_setweight_Estimate)

% Called by Identification_procedure_NE function
% Calls :
% my_data_pets :  to obtain data of experiments created and those available
% in the literature
% pars_init_pets : to obtain parameter values of species and identifies
% those with free = 1 to conduct the estimation
% predict_pets : to conduct the estimation based on the experiments
% selected
% petregr_f : to conduct the estimation

global pets


% set all of the estimation options:
estim_options('default'); 
nmregr_options('report', 0);  % does not report to screen to save time

% Load data and initial parameters
[data, auxData, metaData, txtData, weights] = feval(['mydata_', pets{1}]); 
fieldsData = fieldnames(data); 

% Set weight of datasets not selected == 0 so they do not have an influence in
% the parameter estimation

for i=1:length(fieldsData)
    if ~contains(fieldsData(i), data_setweight_Estimate)
    weights.(fieldsData{i}) = weights.(fieldsData{i}) * 0; 
    end
end


[par, metaPar, txtPar] = feval(['pars_init_',  pets{1}], metaData);
[prdData0, info] = feval(['predict_', pets{1}], par, data, auxData); 
filternm = ['filter_', metaPar.model];


% Identify free = 1  parameters to compute estimation of those parameters
% and deviate the seed with the selected variation
paramNames = fieldnames(par.free);
freeParams = paramNames(cellfun(@(p) par.free.(p) == 1, paramNames));

% Store original parameters
originalPar = par; 

% Vary parameters by variation factor 
for i = 1:length(freeParams)
    pName = freeParams{i};
    par.(pName) = originalPar.(pName) + (originalPar.(pName) * varFactor);
end


% Start estimation
% tStart = tic;
info = 0; % To start the while loop
nCont = 1; % To initialize number of continuations
max_nCont = 5; % Max number of continuations 

while info == 0 && nCont <= max_nCont  
    [par, info, nsteps, fval] = petregr_f(@predict_data_psd, par, data, auxData, weights, filternm);   
    disp(info)
    if info == 1
        break; % Exit the loop if info becomes 1, so if estimation was convergent it exits the loop
    end
    nCont =  nCont + 1; 
end 

% tEnd = toc(tStart);

% Obtain predictions based on new estimation
[prdData, info] = predict_Sargassum_fluitans(par, data, auxData); 
% to change to predict_data_psd
% as it should no include the pets name

% Store results
resultsStructure.par = par;
resultsStructure.prdData = prdData;
resultsStructure.variation_factor = varFactor;
resultsStructure.originalPar = originalPar; 
resultsStructure.freeParams = freeParams; 
resultsStructure.nCont = nCont; 
end

% auxiliary functions :
function [prdData, info] = predict_data_psd(par, data, auxData)
% Predictions, using parameters and data
% Adds pseudodata predictions into predictions structure 

global pets

[prdData, info] = feval(['predict_',pets{1}], par, data, auxData);


end
