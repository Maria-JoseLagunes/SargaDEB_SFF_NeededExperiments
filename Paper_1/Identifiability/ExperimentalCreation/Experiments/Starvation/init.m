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


global pets pars_init_method

dataset_list = {'temperatureSimulations', 'nutrientSimulations'};

pets = {'Sargassum_fluitans'};
% pars_init_method = 1; %1 = from initial parameters %2 from .mat estimation of parameters


dataset = dataset_list{2}; % Choose dataset from list of datasets

pars_initnm = ['pars_init_', pets{1}];
resultsnm   = ['results_', pets{1}, '.mat'];
auxDatanm = ['auxData.', pets{1}];


if pars_init_method == 1
   [par, ~, ~] = feval(pars_initnm, []);
elseif pars_init_method == 2
   load(resultsnm, 'par');
end



% vars_pull(par);
pars_T = [par.T_A; par.T_L; par.T_H; par.T_AL; par.T_AH]; % K




if contains(dataset, 'temperatureSimulations')

       t0 = 0; % h, initial time
       tFinal_days = 10; %d, time of experience
       % tSpinTune_days = 5000; %d, time of experience
       tFinal = (tFinal_days  * 24) ; % hours
       % tSpinTune  = ( tSpinTune_days   * 24); % hours
      
       Ww0 = 2.5 ; %g Ww, initial wet biomass 
       dW0_ash = Ww0 * par.phi_dw;  %g dW, initial dry weight biomass
       dW0 = dW0_ash * (1 - par.x_moist - par.x_ash); 
       M_V_0 = dW0 / (par.w_V +   (par.m_EN_0 * par.w_EN) +   (par.m_EC_0 * par.w_EC)); % mol V, structural intial mass
              
        % Define the temperature values in Celsius
        T_values = [10:3:40]; % Just update this array
        
        % Convert all to Kelvin in one step
        temperature = C2K(T_values);

        
    TC = tempcorr(temperature, par.T_ref, pars_T);  %-
    forcingConditions = get_forc(dataset);
    colorMatrix = colormap(pastel_jet(length(TC)));


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
        % simulationConditions(i).tSpinTune = tSpinTune  ;
        simulationConditions(i).dW0 = dW0; 
        simulationConditions(i).par= par; 
        simulationConditions(i).col = colorMatrix(i,:); 
        simulationConditions(i).photoPeriod = forcingConditions.photoPeriod; 
    end

elseif contains(dataset, 'nutrientSimulations')
     
       t0 = 0; % h, initial time
       tFinal_days = 30; %d, time of experience
       % tSpinTune_days = 5000; %d, time of experience
       tFinal = ( tFinal_days  * 24) ; % hours
       % tSpinTune  = ( tSpinTune_days   * 24); % hours
      
      
     
      
       Ww0 = 6 ; %g Ww, initial wet biomass 
       % dW0 = Ww0 * par.dwratio;  %g dW, initial dry weight biomass
       % M_V_0 = dW0 / (par.w_V +   (par.m_EN_0 * par.w_EN) +   (par.m_EC_0 * par.w_EC)); % mol V, structural intial mass

       dW0_ash = Ww0 * par.phi_dw;  %g dW, initial dry weight biomass
       dW0 = dW0_ash * (1 - par.x_moist - par.x_ash); 
       M_V_0 = dW0 / (par.w_V +   (par.m_EN_0 * par.w_EN) +   (par.m_EC_0 * par.w_EC)); % mol V, structural intial mass
       
       temperature = C2K(25); % K
       TC = tempcorr(temperature, par.T_ref, pars_T);  %-


       forcingConditions = get_forc(dataset);
       [nm nst] =  fieldnmnst_st(forcingConditions);
    

       % To modify nutrient or light conditions create a vector with values
       % to test. For exemple : 
       % forcingConditions.lightIntensity = [100 350 700] * 1e-6 * 3600; 
       % forcingConditions.nitrogen = [1.5 6.5] * 1e-6; 
       % forcingConditions.CO2 = [0.002 0.005] ; 


       forcingConditions.lightIntensity = 200 * 1e-6 * 3600; 
       forcingConditions.nitrogen = 1 * 1e-10; 
       forcingConditions.CO2 = 1 * 1e-10; 


       
      

       numElements = arrayfun(@(i) numel(forcingConditions.(nm{i})), 1:nst);
       [maxVal, maxIdx] = max(numElements);

       if maxVal == 1 %No maximum value
           maxField = NaN; 
       else
            maxField = nm{maxIdx};
       end


      simulationConditions = struct();
      colorMatrix = colormap(pastel_jet(maxVal));
     

        for i = 1:maxVal
             simulationConditions(i).CO2 = getDynamicValue(forcingConditions.CO2, i);
             simulationConditions(i).N = getDynamicValue(forcingConditions.nitrogen, i);
             simulationConditions(i).lightIntensity = getDynamicValue(forcingConditions.lightIntensity, i);
             simulationConditions(i).P = getDynamicValue(forcingConditions.P, i);
             simulationConditions(i).ct = getDynamicValue(TC, i);
             simulationConditions(i).temp = getDynamicValue(temperature, i);
             simulationConditions(i).initialConditions = [par.m_EC_0, par.m_EN_0, getDynamicValue(M_V_0, i)];
             simulationConditions(i).t0 =  t0;
             simulationConditions(i).tFinal = getDynamicValue(tFinal, i);
             % simulationConditions(i).tSpinTune = tSpinTune  ;
             simulationConditions(i).dW0 = getDynamicValue(dW0, i);
             simulationConditions(i).par = par;
             simulationConditions(i).col = colorMatrix(i,:); 
             simulationConditions(i).photoPeriod = forcingConditions.photoPeriod; 

        end
        
        




end






end



%%




