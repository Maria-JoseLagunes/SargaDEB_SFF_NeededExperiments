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

colorMatrix= [0    0.4471    0.7412;
                   0.4667    0.6745    0.1882;
                  1.0000    0.4118    0.1608;
                  0 0 1 ] ;

dataset_list = {'SLL2024', 'MG2023a', 'MG2023b', ...
    'temperatureSimulations', 'nutrientSimulations'};

pets = {'Sargassum_fluitans'};
pars_init_method = 2; %1 = from initial parameters %2 from .mat estimation of parameters


dataset = dataset_list{4}; % Choose dataset from list of datasets

pars_initnm = ['pars_init_', pets{1}];
resultsnm   = ['results_', pets{1}, '.mat'];
auxDatanm = ['auxData.', pets{1}];


if pars_init_method == 1
   [par, metaPar, txtPar] = feval(pars_initnm, []);
elseif pars_init_method == 2
   load(resultsnm, 'par');
end


% vars_pull(par);
pars_T = [par.T_A; par.T_L; par.T_H; par.T_AL; par.T_AH]; % K



% %Modifying elemental composition of structure
% par.n_PV = 0.0012;     free.n_PV = 0;    units.n_PV = '-';    label.n_PV = 'Chemical index of phosphorous in structure from Redfield ratio (106:263:110:16:1)';
% 
% 
% % %Initial conditions of algae
% par.m_EC_0 = 0.434941555120368;   free.m_EC_0 = 0; units.m_EC_0 = 'mol C / mol V'; label.m_EC_0 = 'initial carbon reserve density'; 
% par.m_EN_0 = 0.00212465563703726; free.m_EN_0 = 0; units.m_EN_0 = 'mol N / mol V'; label.m_EN_0 = 'initial nitrogen reserve density'; 

% par.m_EC_0 = 0.002 * 0.0001;   free.m_EC_0 = 0; units.m_EC_0 = 'mol C / mol V'; label.m_EC_0 = 'initial carbon reserve density'; 
% par.m_EN_0 = 0.01 * 0.0001 ; free.m_EN_0 = 0; units.m_EN_0 = 'mol N / mol V'; label.m_EN_0 = 'initial nitrogen reserve density'; 
% % 

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

dW0 = Ww0 * par.dwratio;  %g dW, initial dry weight biomass
M_V_0 = dW0 / (par.w_V +   (par.m_EN_0 * par.w_EN) +   (par.m_EC_0 * par.w_EC)); % mol V, structural intial mass

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
colorMatrix = colormap(pastel_jet(length(TC)));




for i = 1:length(TC)
    simulationConditions(i).data = data;
    simulationConditions(i).auxData = auxData;
    simulationConditions(i).CO2 = forcingConditions.CO2(i) ;
    simulationConditions(i).N = forcingConditions.nitrogen(i) ;
    simulationConditions(i).lightIntensity = forcingConditions.lightIntensity(i) ;
    simulationConditions(i).P = forcingConditions.P;
    simulationConditions(i).ct = TC(i);
    simulationConditions(i).temp = temperature(i);
    simulationConditions(i).initialConditions = [par.m_EC_0, par.m_EN_0, M_V_0(i)];
    simulationConditions(i).t0 = t0;
    simulationConditions(i).tFinal = tFinal(i) ; 
    simulationConditions(i).dW0 = dW0(i); 
    simulationConditions(i).par= par; 
    simulationConditions(i).col = colorMatrix(i,:); 
    simulationConditions(i).photoPeriod = forcingConditions.photoPeriod; 
end


elseif contains(dataset, 'temperatureSimulations')

  
     
       t0 = 0; % h, initial time
       tFinal_days = 10; %d, time of experience
       tSpinTune_days= 100; 
       tFinal = ( tFinal_days  * 24) ; % hours
       tSpinTune  = ( tSpinTune_days   * 24);
      
       Ww0 = 6; %g Ww, initial wet biomass 
       dW0 = Ww0 * par.dwratio;  %g dW, initial dry weight biomass
       M_V_0 = dW0 / (par.w_V +   (par.m_EN_0 * par.w_EN) +   (par.m_EC_0 * par.w_EC)); % mol V, structural intial mass
       
      
    
            
        % Define the temperature values in Celsius
        % T_values = [22 24 28 30]; % Just update this array!
        T_values = [25]; % Just update this array!
        
        % Convert all to Kelvin in one step
        temperature = C2K(T_values);
        
       
    

    
             
        TC = tempcorr(temperature, par.T_ref, pars_T);  %-
        forcingConditions = get_forc(dataset);
        colorMatrix = colormap(get_idf_palette(length(TC)));


    for i = 1:length(TC)
 

    simulationConditions(i).CO2 = forcingConditions.CO2 ;
    simulationConditions(i).N = forcingConditions.nitrogen ;
    simulationConditions(i).lightIntensity = forcingConditions.lightIntensity ;
    simulationConditions(i).P = forcingConditions.P;
    simulationConditions(i).ct = TC(i);
    simulationConditions(i).temp = temperature(i);
    simulationConditions(i).initialConditions = [par.m_EC_0, par.m_EN_0, M_V_0];
    simulationConditions(i).t0 = t0;
    simulationConditions(i).tFinal = tFinal ; 
    simulationConditions(i).tSpinTune = tSpinTune  ;
    simulationConditions(i).dW0 = dW0; 
    simulationConditions(i).par= par; 
    simulationConditions(i).col = colorMatrix(i,:); 
    simulationConditions(i).photoPeriod = forcingConditions.photoPeriod; 
    end
     
    
elseif contains(dataset, 'nutrientSimulations')


       t0 = 0; % h, initial time
       tFinal_days = 20; %d, time of experience
       tSpinTune_days= 5; 
       tFinal = ( (tFinal_days )  * 24) ; % hours
       tSpinTune  = ( (tSpinTune_days)  * 24);
      
     
      
       Ww0 = 6 ; %g Ww, initial wet biomass 
       dW0 = Ww0 * par.dwratio;  %g dW, initial dry weight biomass
       M_V_0 = dW0 / (par.w_V +   (par.m_EN_0 * par.w_EN) +   (par.m_EC_0 * par.w_EC)); % mol V, structural intial mass
       
       temperature = C2K(25); % K
       TC = tempcorr(temperature, par.T_ref, pars_T);  %-


       forcingConditions = get_forc(dataset);
       [nm nst] =  fieldnmnst_st(forcingConditions);
    

       % To modify nutrient or light conditions create a vector with values
       % to test. For exemple : 
       forcingConditions.lightIntensity = [300] * 1e-6 * 3600; 
       forcingConditions.nitrogen = [1] * 1e-6; 
       forcingConditions.CO2 = [0.0002] ; 

      

       numElements = arrayfun(@(i) numel(forcingConditions.(nm{i})), 1:nst);
       [maxVal, maxIdx] = max(numElements);

       if maxVal == 1 %No maximum value
           maxField = NaN; 
       else
            maxField = nm{maxIdx};
       end


      simulationConditions = struct();
     
        
        for i = 1:maxVal
             simulationConditions(i).CO2 = getDynamicValue(forcingConditions.CO2, i);
             simulationConditions(i).N = getDynamicValue(forcingConditions.nitrogen, i);
             simulationConditions(i).lightIntensity = getDynamicValue(forcingConditions.lightIntensity, i);
             simulationConditions(i).P = getDynamicValue(forcingConditions.P, i);
             simulationConditions(i).ct = getDynamicValue(TC, i);
             simulationConditions(i).temp = getDynamicValue(temperature, i);
             simulationConditions(i).initialConditions = [par.m_EC_0, par.m_EN_0, getDynamicValue(M_V_0, i)];
             simulationConditions(i).t0 = 1;
             simulationConditions(i).tFinal = getDynamicValue(tFinal, i);
             simulationConditions(i).tSpinTune = tSpinTune  ;
             simulationConditions(i).dW0 = getDynamicValue(dW0, i);
             simulationConditions(i).par = par;
             simulationConditions(i).col = colorMatrix(i,:); 
             simulationConditions(i).photoPeriod = forcingConditions.photoPeriod; 

        end
        




end






end



%%




