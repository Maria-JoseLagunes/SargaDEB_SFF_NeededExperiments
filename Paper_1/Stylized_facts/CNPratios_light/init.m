% ---------------------------------------------------------------------
%% Initital conditions for iteration and time vector
% Maria José Lagunes
% First created : 2024/03/12
%Define initial conditions for simulation and the time of simulation
% Input
% modelParameters --> calling set_pars function
%
% ---------------------------------------------------------------------
function simulationConditions = init

global pets

% colorMatrix= [0    0.4471    0.7412;
%                    0.4667    0.6745    0.1882;
%                   1.0000    0.4118    0.1608;
%                   0 0 1 ] ;

dataset_list = {'SLL2024', 'MG2023a', 'MG2023b', ...
    'temperatureSimulations', 'nutrientSimulations'};

pets = {'Sargassum_fluitans'};
pars_init_method = 1; %1 = from initial parameters %2 from .mat estimation of parameters


dataset = dataset_list{5}; % Choose dataset from list of datasets

pars_initnm = ['pars_init_', pets{1}];
resultsnm   = ['results_', pets{1}, '.mat'];
auxDatanm = ['auxData.', pets{1}];


if pars_init_method == 1
   [par, metaPar, txtPar] = feval(pars_initnm, []);
elseif pars_init_method == 2
   load(resultsnm, 'par');
end


vars_pull(par);
pars_T = [par.T_A; par.T_L; par.T_H; par.T_AL; par.T_AH]; % K



%Modifying elemental composition of structure
par.n_PV = 0.0012;     free.n_PV = 0;    units.n_PV = '-';    label.n_PV = 'Chemical index of phosphorous in structure from Redfield ratio (106:263:110:16:1)';



if contains(dataset,"SLL2024") || contains(dataset,"MG2023a") || contains(dataset,"MG2023b")
[data, auxData, metaData, txtData, weights] = mydata_pets;

tFinal = [];
Ww0 = [];
temperature = []; 
N = [];
t0 = 1; 

% Final time
timeFields = fieldnames(auxData.(pets{1}).time);

% Loop through each field in timeFields for time data
for j = 1:length(timeFields)
    fieldName = timeFields{j};
    
    % Check if the field name contains the current dataset
    if contains(fieldName, dataset)
        % Extract the time data and convert it to hours
        tFinal = [tFinal; auxData.(pets{1}).time.(fieldName) * 24]; 
    end
end

%Initial conditions wet weight
Ww0Fields = fieldnames(auxData.(pets{1}).Ww0);

for k = 1:length(Ww0Fields)
    fieldName = Ww0Fields{k};
    
    % Check if the field name contains the current dataset
    if contains(fieldName, dataset)
        % Extract the initial wet weight 
        Ww0= [Ww0; auxData.(pets{1}).Ww0.(fieldName)]; 
    end
end

dW0 = Ww0 * dwratio;  %g dW, initial dry weight biomass
M_V_0 = dW0 / (w_V+   (m_EN_0*w_EN) +   (m_EC_0*w_EC)); % mol V, structural intial mass


%Temperature
tempFields = fieldnames(auxData.(pets{1}).temp);

for l = 1:length(tempFields)
    fieldName = tempFields{l};
    
    % Check if the field name contains the current dataset
    if contains(fieldName, dataset)
        % Extract the initial wet weight 
        temperature = [temperature; auxData.(pets{1}).temp.(fieldName)]; 
    end
end

TC = tempcorr(temperature, par.T_ref, pars_T);  %-
forcingConditions = get_forc(dataset); 

for i = 1:length(TC)
    simulationConditions(i).data = data;
    simulationConditions(i).auxData = auxData;
    simulationConditions(i).CO2 = forcingConditions.CO2(i) ;
    simulationConditions(i).N = forcingConditions.nitrogen(i) ;
    simulationConditions(i).lightIntensity = forcingConditions.lightIntensity(i) ;
    simulationConditions(i).P = forcingConditions.P;
    simulationConditions(i).ct = TC(i);
    simulationConditions(i).temp = temperature(i);
    simulationConditions(i).initialConditions = [m_EC_0, m_EN_0, M_V_0(i)];
    simulationConditions(i).t0 = t0;
    simulationConditions(i).tFinal = tFinal(i) ; 
    simulationConditions(i).dW0 = dW0(i); 
    simulationConditions(i).par= par; 
    simulationConditions(i).col = colorMatrix(i,:); 
end


elseif contains(dataset, 'temperatureSimulations')

       t0 = 1; % h, initial time
       tFinal_days = 10; %d, time of experience
       tFinal = ((tFinal_days + 1) * 24) ; % hours
     
      
     
      
       Ww0 = 2.5 ; %g Ww, initial wet biomass 
       dW0 = Ww0 * dwratio;  %g dW, initial dry weight biomass
       M_V_0 = dW0 / (w_V+   (m_EN_0*w_EN) +   (m_EC_0*w_EC)); % mol V, structural intial mass

       
       
      
    
        % Define the temperature values
        T1 = C2K(22); % K
        T2 = C2K(24); % K
        T3 = C2K(26); % K
        T4 = C2K(28); % K
        % T5 = C2K(30); % K
        % T6 = C2K(32); % K
        

        % Initialize an empty array for selected temperatures
        temperature = [];
        
        % Dynamically add temperatures based on conditions
        if exist('T1', 'var')
            temperature(end+1) = T1; % Add T1
        end
        if exist('T2', 'var')
            temperature(end+1) = T2; % Add T2
        end
        if exist('T3', 'var')
            temperature(end+1) = T3; % Add T3
        end
        if exist('T4', 'var')
            temperature(end+1) = T4; % Add T4
        end
         if exist('T5', 'var')
            temperature(end+1) = T5; % Add T5
        end

         if exist('T6', 'var')
            temperature(end+1) = T6; % Add T6
        end



        
    TC = tempcorr(temperature, par.T_ref, pars_T);  %-
    forcingConditions = get_forc(dataset);
    colorMatrix = colormap(abyss(length(TC)));

    for i = 1:length(TC)
 

    simulationConditions(i).CO2 = forcingConditions.CO2 ;
    simulationConditions(i).N = forcingConditions.nitrogen ;
    simulationConditions(i).lightIntensity = forcingConditions.lightIntensity ;
    simulationConditions(i).photoPeriod = forcingConditions.photoPeriod ;
    simulationConditions(i).P = forcingConditions.P;
    simulationConditions(i).ct = TC(i);
    simulationConditions(i).temp = temperature(i);
    simulationConditions(i).initialConditions = [m_EC_0, m_EN_0, M_V_0];
    simulationConditions(i).t0 = t0;
    simulationConditions(i).tFinal = tFinal ; 
    simulationConditions(i).dW0 = dW0; 
    simulationConditions(i).par= par; 
    simulationConditions(i).col = colorMatrix(i,:); 
    end
     
    %    t0 = 1; % h, initial time
    %    tFinal_days = 10 + 1; %d, time of experience
    %    tFinal = tFinal_days  * 24; % hours
    % 
    % 
    %    Ww0 = 2.5 ; %g Ww, initial wet biomass 
    %    dW0 = Ww0 * dwratio;  %g dW, initial dry weight biomass
    %    M_V_0 = dW0 / (w_V+   (m_EN_0*w_EN) +   (m_EC_0*w_EC)); % mol V, structural intial mass
    % 
    % 
    % 
    % 
    % 
    %     % Define the temperature values
    %     T1 = C2K(23); % K
    %     T2 = C2K(25); % K
    %     T3 = C2K(27); % K
    %     T4 = C2K(29); % K
    % 
    %     % Initialize an empty array for selected temperatures
    %     temperature = [];
    % 
    %     % Dynamically add temperatures based on conditions
    %     if exist('T1', 'var')
    %         temperature(end+1) = T1; % Add T1
    %     end
    %     if exist('T2', 'var')
    %         temperature(end+1) = T2; % Add T2
    %     end
    %     if exist('T3', 'var')
    %         temperature(end+1) = T3; % Add T3
    %     end
    %     if exist('T4', 'var')
    %         temperature(end+1) = T4; % Add T4
    %     end
    % 
    % 
    % 
    % TC = tempcorr(temperature, par.T_ref, pars_T);  %-
    % forcingConditions = get_forc(dataset);
    % 
    % for i = 1:length(TC)
    % 
    % 
    % simulationConditions(i).CO2 = forcingConditions.CO2 ;
    % simulationConditions(i).N = forcingConditions.nitrogen ;
    % simulationConditions(i).lightIntensity = forcingConditions.lightIntensity ;
    % simulationConditions(i).P = forcingConditions.P;
    % simulationConditions(i).ct = TC(i);
    % simulationConditions(i).temp = temperature(i);
    % simulationConditions(i).initialConditions = [m_EC_0, m_EN_0, M_V_0];
    % simulationConditions(i).t0 = t0;
    % simulationConditions(i).tFinal = tFinal ; 
    % simulationConditions(i).dW0 = dW0; 
    % simulationConditions(i).par= par; 
    % simulationConditions(i).col = colorMatrix(i,:); 
    % end
    
    
elseif contains(dataset, 'nutrientSimulations')
      t0 = 1; % h, initial time
       tFinal_days = 10; %d, time of experience
       tFinal = ((tFinal_days + 1) * 24) ; % hours
     
      
     
      
       Ww0 = 2.5 ; %g Ww, initial wet biomass 
       dW0 = Ww0 * dwratio;  %g dW, initial dry weight biomass
       M_V_0 = dW0 / (w_V+   (m_EN_0*w_EN) +   (m_EC_0*w_EC)); % mol V, structural intial mass
       
       
       temperature = C2K(23); % K
       TC = tempcorr(temperature, par.T_ref, pars_T);  %-


       forcingConditions = get_forc(dataset);
       [nm nst] =  fieldnmnst_st(forcingConditions);
    

       % To modify nutrient or light conditions create a vector with values
       % to test. For exemple : 
       forcingConditions.lightIntensity = [200 500 700 1200] * 1e-6 * 3600; 
       % forcingConditions.nitrogen = [1.5 6.5] * 1e-6; 
       % forcingConditions.CO2 = [0.002 0.005] ; 

      
      
       numElements = arrayfun(@(i) numel(forcingConditions.(nm{i})), 1:nst);
       [maxVal, maxIdx] = max(numElements);

       if maxVal == 1 %No maximum value
           maxField = NaN; 
       else
            maxField = nm{maxIdx};
       end


      simulationConditions = struct();
     
       colorMatrix = colormap(abyss(maxVal));
       
        for i = 1:maxVal
             simulationConditions(i).CO2 = getDynamicValue(forcingConditions.CO2, i);
             simulationConditions(i).N = getDynamicValue(forcingConditions.nitrogen, i);
             simulationConditions(i).lightIntensity = getDynamicValue(forcingConditions.lightIntensity, i);
             simulationConditions(i).P = getDynamicValue(forcingConditions.P, i);
              simulationConditions(i).photoPeriod = forcingConditions.photoPeriod ;
             simulationConditions(i).ct = getDynamicValue(TC, i);
             simulationConditions(i).temp = getDynamicValue(temperature, i);
             simulationConditions(i).initialConditions = [m_EC_0, m_EN_0, getDynamicValue(M_V_0, i)];
             simulationConditions(i).t0 = 1;
             simulationConditions(i).tFinal = getDynamicValue(tFinal, i);
             simulationConditions(i).dW0 = getDynamicValue(dW0, i);
             simulationConditions(i).par = par;
             simulationConditions(i).col = colorMatrix(i,:); 

        end
        




end






end



%%




