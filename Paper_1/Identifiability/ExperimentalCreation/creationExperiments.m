%% Create experimental data based on conditions specified for each pseudoexperiment
% Created by : Maria José Lagunes


% Called by : Not function, output is saved to a .mat file with the pseudoexperiments observables, forcing conditions and intial conditions
% Calls :
% my_data_pets :  to obtain data of experiments created and those available
% in the literature
% pars_init_pets : to obtain parameter values of species and identifies
% those with free = 1 to conduct the estimation
% predict_pets : to conduct the estimation based on the experiments
% selected
% petregr_f : to conduct the estimation


clc; clear; close all;

addpath(genpath('/home/LAGUNES/Documents/GitHub/SargaDEB_working/DEBtool_M-master'))
myPath = '/home/LAGUNES/Documents/GitHub/SargaDEB_working/matlab/multi_DEB/Paper_1';
addpath(genpath(myPath))
addpath(genpath(strcat(myPath,'/Identifiability')));
workPath = strcat(myPath, '/Identifiability/ExperimentalCreation/Experiments'); 



%% 1. Create experimental data 
% !!! ATTENTION !!! you must be inside the directory of experience being
% created
%
% Write the experiment you want to create by writing
%
% * TC_Ww : if modification of data to evaluate wet weight evolution at
% different temperatures
% * TC_Ww_twoDays : if modification of data to evaluate wet weight evolution at
% different temperatures with only initial and final weight
% * Nutrient_uptake : if simulation of nutrient uptake SU for nitrogen

% prompt = "What experimental data are you creating? ";
% experimentName = input(prompt); 

% !!! ATTENTION !! Now codes are automated to perform the for experiments
% each time something changes 
global pets
pets = {'Sargassum_fluitans'};

listExperiments = {'TC_Ww', 'N_uptake', 'TC_Ww_twoDays', 'Starvation'};


for i = 1:length(listExperiments)
    cd(workPath)
    experimentName = listExperiments{i};
    if ~exist(experimentName, 'dir')
       mkdir(experimentName)
    end
    addpath(experimentName)
    cd(experimentName)
    

    
    if contains(experimentName, 'TC_Ww') || contains(experimentName, 'Starvation') ...
            || contains(experimentName, 'TC_Ww_twoDays')
        [simu, obs] = run_simulation(); %Function called to perform the complete DEB multireserves ...
        % model based on a preivously edited init file
    elseif contains(experimentName, 'N_uptake')
        [simu, obs] = run_Nuptake(0:7:50,[22 26 28]); %Function called to perform assimilation box of the DEB model for N uptake
    end
    save_simulation_data(simu, obs, pets); %Function created to save simu and observation for later use as mydata

end

%% 2. Saving experimental data
% %Save to use in anaylsis
% save('simu.mat', 'simu'); 
% save('obs.mat', 'obs');

% Get today's date and time in MATLAB's serial date number format
% Save to replicate base on data generated



%% 3. Importing files from date in format 'ddmmyyyy'
% Choose date to import, now or write date to import
cd(workPath)

formattedDate = datestr(now, 'ddmmyyyy'); %To import last experiments created, change this manually if another date is to evaluate


data = struct(); 
simulationConditions = struct();
% listExperiments = {'TC_Ww'};
% listExperiments = {'TC_Ww', 'N_uptake', 'TC_Ww_twoDays'};
listExperiments = {'TC_Ww', 'N_uptake', 'TC_Ww_twoDays', 'Starvation'};
dataExperiment = struct();

% Creating a data structure with the needed information for each dataset
% to be included in the mydata file and to be used in the predict file

for expIdx = 1:length(listExperiments)
    expName = listExperiments{expIdx};
    load(fullfile(expName, strcat("obs_", pets,".mat")));
    load(fullfile(expName, strcat("simu_", pets,".mat")));
    data.(expName)= obs;
    simulationConditions.(expName) = simu;
    
    if strcmp(expName, 'TC_Ww')
       setDays = [0 2 4 6 8 10];
       for i = 1:length(data.(expName))
            dataTemp = data.(expName)(i);
            simuTemp = simulationConditions.(expName)(i); 
            
            tempName = sprintf("%s_%d", expName, i);
            nline = [1 ((setDays(2:length(setDays)))* 24)+1]; 
            nline = nline(nline <= length(dataTemp.Ww));  % remove indices that are out of bounds
            dataExperiment.(tempName).timePoints = setDays; % in hours
            dataExperiment.(tempName).Ww = dataTemp.Ww(nline); 
            dataExperiment.(tempName).Ww0 = dataTemp.Ww_0;
            dataExperiment.(tempName).time = round(simuTemp.tFinal/24); 

            dataExperiment.(tempName).temp = simuTemp.temp; 
            dataExperiment.(tempName).light = simuTemp.lightIntensity; 
            dataExperiment.(tempName).nitrogen = simuTemp.N; 
            dataExperiment.(tempName).CO2 = simuTemp.CO2; 
            dataExperiment.(tempName).phosphorous = simuTemp.P; 
            dataExperiment.(tempName).photoPeriod = simuTemp.photoPeriod; 

        end
    

    elseif strcmp(expName, 'TC_Ww_twoDays')
       setDays = [0 10];
       for i = 1:length(data.(expName))
               dataTemp = data.(expName)(i);
            simuTemp = simulationConditions.(expName)(i); 
            
            tempName = sprintf("%s_%d", expName, i);
            nline = [1 ((setDays(2:length(setDays)))* 24)+1]; 
            nline = nline(nline <= length(dataTemp.Ww));  % remove indices that are out of bounds
            dataExperiment.(tempName).timePoints = setDays; % in hours
            dataExperiment.(tempName).Ww = dataTemp.Ww(nline); 
            dataExperiment.(tempName).Ww0 = dataTemp.Ww_0;
            dataExperiment.(tempName).time = round(simuTemp.tFinal/24); 

            dataExperiment.(tempName).temp = simuTemp.temp; 
            dataExperiment.(tempName).light = simuTemp.lightIntensity; 
            dataExperiment.(tempName).nitrogen = simuTemp.N; 
            dataExperiment.(tempName).CO2 = simuTemp.CO2; 
            dataExperiment.(tempName).phosphorous = simuTemp.P; 
            dataExperiment.(tempName).photoPeriod = simuTemp.photoPeriod; 
        end

    elseif strcmp(expName, 'N_uptake')
        for i = 1:length(data.(expName))
            dataTemp = data.(expName)(i);
            simuTemp = simulationConditions.(expName)(i); 
            tempName = sprintf("%s_%d", expName, i);
            dataExperiment.(tempName).NitrateAssimilation = dataTemp.j_EN_A_inmicroMol;  % micro mol N g dW-1 h-1
            dataExperiment.(tempName).NitrateConcentration = simuTemp.N_microConcentration;
            dataExperiment.(tempName).temp = simuTemp.temp;
            dataExperiment.(tempName).Ww0 = simuTemp.Ww0;
        end


     elseif strcmp(expName, 'Starvation')
       setDays = [0 10 20 30];
       for i = 1:length(data.(expName))
            dataTemp = data.(expName)(i);
            simuTemp = simulationConditions.(expName)(i); 
            
            tempName = sprintf("%s_%d", expName, i);
            nline = [1 ((setDays(2:length(setDays)))* 24)+1]; 
            nline = nline(nline <= length(dataTemp.Ww));  % remove indices that are out of bounds
            dataExperiment.(tempName).timePoints = setDays; % in hours
            dataExperiment.(tempName).Ww = dataTemp.Ww(nline); 
            dataExperiment.(tempName).Ww0 = dataTemp.Ww_0;
            dataExperiment.(tempName).CNtotalRatio = dataTemp.CNtotalRatio(nline); 
            dataExperiment.(tempName).CPtotalRatio = dataTemp.CPtotalRatio(nline); 
            dataExperiment.(tempName).NPtotalRatio = dataTemp.NPtotalRatio(nline); 
            dataExperiment.(tempName).Ww0 = dataTemp.Ww_0;
            dataExperiment.(tempName).time = round(simuTemp.tFinal/24); 


            dataExperiment.(tempName).temp = simuTemp.temp; 
            dataExperiment.(tempName).light = simuTemp.lightIntensity; 
            dataExperiment.(tempName).nitrogen = simuTemp.N; 
            dataExperiment.(tempName).CO2 = simuTemp.CO2; 
            dataExperiment.(tempName).phosphorous = simuTemp.P; 
            dataExperiment.(tempName).photoPeriod = simuTemp.photoPeriod; 
       end
    end
end

%% Saving mat of experiment data


save(fullfile(workPath, 'dataExperiment'), 'dataExperiment')

