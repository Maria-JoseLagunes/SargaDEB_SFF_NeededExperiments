function forcingVariables = get_forc(forcingConditions)
global pets

nitrogen = [];
light = [];
CO2 = []; 
photoPeriod = [];


if contains(forcingConditions,"SLL2024") || contains(forcingConditions,"MG2023a") || contains(forcingConditions,"MG2023b")
[data, auxData, metaData, txtData, weights] = mydata_pets;
% Nitrogen concentration NH4 + NO3
NFields = fieldnames(auxData.(pets{1}).nitrogen);

for l = 1:length(NFields)
    fieldName = NFields{l};
    
    % Check if the field name contains the current dataset
if contains(fieldName, forcingConditions)
        % Extract the initial wet weight 
        nitrogen = [nitrogen; auxData.(pets{1}).nitrogen.(fieldName)]; 
    end
end


% CO2
CO2Fields = fieldnames(auxData.(pets{1}).CO2);

for l = 1:length(CO2Fields)
    fieldName = CO2Fields{l};
    
    % Check if the field name contains the current dataset
    if contains(fieldName, forcingConditions)
        % Extract the initial wet weight 
        CO2 = [CO2; auxData.(pets{1}).CO2.(fieldName)]; 
    end
end

% Light
lightFields = fieldnames(auxData.(pets{1}).light);

for l = 1:length(lightFields)
    fieldName = lightFields{l};
    
    % Check if the field name contains the current dataset
    if contains(fieldName, forcingConditions)
        % Extract the initial wet weight 
        light = [light; auxData.(pets{1}).light.(fieldName)]; 
    end
end


% Photoperiod
photoPeriodFields = fieldnames(auxData.(pets{1}).photoPeriod);

for l = 1:length(photoPeriodFields)
    fieldName = photoPeriodFields{l};
    
    % Check if the field name contains the current dataset
    if contains(fieldName, forcingConditions)
        % Extract the initial wet weight 
        photoPeriod = [photoPeriod; auxData.(pets{1}).photoPeriod.(fieldName)]; 
    end
end

forcingVariables.CO2 = CO2 ; % mol DIC
forcingVariables.nitrogen = nitrogen ; % mol NO3 and NH4 L-1
forcingVariables.lightIntensity = light; % mol gamma m-2 h-1
forcingVariables.photoPeriod = photoPeriod; % -
forcingVariables.P =  0.10 * 1e-6; % mol P043- L-1


elseif contains(forcingConditions, "temperatureSimulations")
    lightIntensity = 500;
    forcingVariables.lightIntensity = lightIntensity * 1e-6 * 3600; % mol gamma m-2 h-1
    forcingVariables.nitrogen = 3 * 1e-6 ; % mol N L-1
    forcingVariables.CO2 = 0.002; % mol C L-1 
    forcingVariables.P = 0.2 * 1e-6 ; % mol PO43- L-1 
    forcingVariables.photoPeriod = 1; % -

 elseif forcingConditions == "nutrientSimulations"    
    lightIntensity = 500;
    forcingVariables.lightIntensity = lightIntensity * 1e-6 * 3600; % mol gamma m-2 h-1
    forcingVariables.nitrogen = 1.5 * 1e-6 ; % mol N L-1
    forcingVariables.CO2 = 0.002; % mol C L-1 
    forcingVariables.P = 0.10 * 1e-6 ; % mol PO43- L-1 
    forcingVariables.photoPeriod = 1; % -

end

