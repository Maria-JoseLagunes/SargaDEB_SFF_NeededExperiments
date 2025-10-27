%Work
addpath(genpath('/home/LAGUNES/Documents/GitHub/SargaDEB_working/DEBtool_M-master'))
cd /home/LAGUNES/Documents/GitHub/SargaDEB_working/matlab/multi_DEB/Paper_1
addpath(genpath("/home/LAGUNES/Documents/GitHub/SargaDEB_working/matlab/multi_DEB/Paper_1/"))

%Laure
  % addpath(genpath('/Users/laurepecquerie/Documents/GitHub/SargaDEB_working/DEBtool_M-master'))
  % cd /Users/laurepecquerie/Documents/GitHub/SargaDEB_working/matlab/multi_DEB/Paper_1/
  % addpath(genpath("/Users/laurepecquerie/Documents/GitHub/SargaDEB_working/matlab/multi_DEB/Paper_1"))
 % 
%Home
% addpath(genpath('C:\Users\Majo\Documents\PhD\SargaDEB_home\DEBtool_M-master'))
% cd C:\Users\Majo\Documents\PhD\SargaDEB_home\matlab\multi_DEB\Cleaningcodes
% addpath("C:\Users\Majo\Documents\PhD\SargaDEB_home\matlab\multi_DEB\Cleaningcodes")


cd("Identifiability/ExperimentalCreation")
%% Selection of datasets not to estimate
clc; clear; close all;
%% Datasets available for estimation
% The PE refers to : pseudo numerical experiment
%
% 1 : HT2023
% 2 : HS1987
% 3 : MG2023a
% 4 : MG2023b
% 5 : VE2023
% 6 : PE : TC_Ww temperatures (10:3:40) with datapoints at [0 2 4 6 8 10] days 
% 7 : PE : N uptake at three temperatures
% 8 : PE : TC_Ww temperatures (10:3:40) with datapoints at [0 10] days 
% 9 : PE : Starvation with datapoints at [0 10 20 30] days Ww and CNtotalratio at 25°C

%% Estimation
% Selection of which data sets to use for the estimation procedure. Those
% not included will be set with weight = 0 in the EstimationProcedure_NE
% function
global pets

list_datasets = ["HT2023", "HS1987", "MG2023a","MG2023b","VE2023", "TC_Ww",...
    "N_uptake", "TC_WW_Twodays", 'Starvation']; %List of existent datasets on mydata_pets



pets={'Sargassum_fluitans'};
%% RUN ONLY IF NEW ANALYSIS BEING MADE
varFactor_list = [0.05 0.1 0.2] ; % Variation factor, 10%
varFactor_list = 0.1; 

for j = 1:length(varFactor_list)
    varFactor = varFactor_list(j);

    all_which_data = {
    [3 4]
    [1 3 4]
    [1 3 4 5]
    [3 4 5]
    [1 3 4 5 7]
    [1 3 4 5 9]
    [1 3 4 5 7 9]
    [3 4 5 7 9]
    [1 3 4 7 9]
    };

    % all_which_data = {
    % [3 4]
    % [3 4 5]
    % [3 4 5 7 9]
    % };
    % 

    
    
    for i = 1:length(all_which_data)
        which_data = all_which_data{i};
       % Calls identification procedure according to each which data selected
       [originalPar, positiveParameters, negativeParameters, SSE_structure] = IdentificationProcedure_NE(varFactor, which_data, list_datasets); 
    end
end

% Sound when identification procedure for each dataset is finished 
load handel
sound(y,Fs)
%% Plots if direct analysis

clf


for i = 1:length(all_which_data)
    which_data = all_which_data{i};
    config_label = strjoin(string(which_data), '_');
    config_label = string(length(positiveParameters.freeParams)) + '_' + config_label;
    filename = strjoin(['SSE_Analysis/SSE_results', config_label , '.mat']);
    if isfile(filename)
        load(filename); % Should load SSE_results
        plot_SSE_variation(SSE_results.SSE_structure, which_data, list_datasets, i);
    else
    end
end

