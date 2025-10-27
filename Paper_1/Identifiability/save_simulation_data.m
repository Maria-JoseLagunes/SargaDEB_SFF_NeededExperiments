function save_simulation_data(simu, obs, pets)

% Create timestamped folder
    formattedDate = datestr(now, 'ddmmyyyy');
    folderName = strcat("data_", formattedDate);
    mkdir(folderName);
    folderPath = fullfile(pwd, folderName);
    folderPath_savelast = fullfile(pwd); 

    % Save files with experiment name
    simulationData_fileName = strcat('simu_', pets{1}, '.mat'); 
    observablesData_fileName = strcat('obs_', pets{1}, '.mat'); 

    save(fullfile(folderPath_savelast, simulationData_fileName), 'simu');
    save(fullfile(folderPath_savelast, observablesData_fileName), 'obs');

     save(fullfile(folderPath, simulationData_fileName), 'simu');
    save(fullfile(folderPath, observablesData_fileName), 'obs');


    fprintf("Simulation data saved in %s\n", folderPath);

end