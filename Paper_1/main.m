%% Main file to do simulations
% Call for solver 
% Maria José Lagunes
% Last modified : 11/12/2024
% ---------------------------------------------------------------------

clear all

% Lab

addpath(genpath('/home/LAGUNES/Documents/GitHub/SargaDEB_working/DEBtool_M-master'))
addpath(genpath('/home/LAGUNES/Documents/GitHub/SargaDEB_working/matlab/multi_DEB/Paper_1'))
cd /home/LAGUNES/Documents/GitHub/SargaDEB_working/matlab/multi_DEB/Paper_1/

%Laure

% addpath(genpath('/Users/laurepecquerie/Documents/GitHub/DEBtool_M'))
% addpath(genpath('/Users/laurepecquerie/Documents/GitHub/SargaDEB_working/matlab/multi_DEB/Paper_1'))
% cd /Users/laurepecquerie/Documents/GitHub/SargaDEB_working/matlab/multi_DEB/Paper_1/

% Home 
% addpath(genpath('C:\Users\Majo\Documents\PhD\SargaDEB_home\DEBtool_M-master'))
% addpath(genpath('C:\Users\Majo\Documents\PhD\SargaDEB_home\matlab\multi_DEB\Cleaningcodes'))
% cd C:\Users\Majo\Documents\PhD\SargaDEB_home\matlab\multi_DEB\Cleaningcodes

%% 1 Call integration and perform it to exit integration values
% Initialize iteration counter
global pets
% close all; clear; 

% Initialize time, parameters, initial conditions, ...
% Initialize time, parameters, initial conditions, ...
simu = init;

[t, mECENV, J_struct] = indiv(simu); 
simu.t = t; 
simu.mECENV = mECENV; 
simu.CIs = mECENV(1,:); 
simu.J_struct = J_struct; 


obs = get_obs(mECENV,J_struct, simu);
get_plots_2(t,mECENV,J_struct,obs, simu);



%Calculate state variables
for i = 1 : length(simu)

    [t, mECENV, J_struct] = indiv(simu(i));
    simu(i).t = t; 
    simu(i).mECENV = mECENV; 
    simu(i).CIs = mECENV(1,:); 
    simu(i).J_struct = J_struct; 


    % Observables
    obs(i) = get_obs(mECENV,J_struct, simu(i));


    %Plots
    
    get_plots_2(t,mECENV,J_struct,obs(i), simu(i));



end 



%% Saving the plots
%Get today's date and time in MATLAB's serial date number format
currentDate = now;

% Convert the serial date number to a formatted date string
formattedDate = datestr(currentDate, 'ddmmyyyy');
directoryPath =  pwd; 

folderName = strcat("data",'_',formattedDate); 
mkdir(fullfile(directoryPath,folderName))

folder_path = fullfile(directoryPath,folderName); 


simulationData_fileName = strcat('simulationData_','.mat'); 
observablesData_fileName = strcat('obsData_','.mat'); 


save(fullfile(folder_path, simulationData_fileName), 'simu');
save(fullfile(folder_path, observablesData_fileName), 'obs');

save(fullfile(pwd,'simu'), 'simu');
save(fullfile(pwd,'obs'), 'obs');


