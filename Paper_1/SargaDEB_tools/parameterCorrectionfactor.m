%% VENOLIA %%

pars.T_A = 6314.3; 
pars.T_AH = 18702; 
pars.T_AL = 4391.9; 
pars.T_H = 286.536; 
pars.T_L = 273.15; 
pars.T_ref = 293.15; 


% T = -5:0.05:35; %°C, temperatures tested, 2.6 to have maximum 12 points for estimation procedure of DEBtoolM
T = [20 25]; 
% pars_T = [T_A; T_L; T_H; T_AL; T_AH]'; % K
temperatureEffect =  zeros(length(T),1);

for l = 1:length(T)
        temp = T(l); %°C
        ct = Arrhenius_equation(pars, temp, 'simple');
        temperatureEffect(l,1) = ct; % -

end


% Create plot of function cT 
plot(T,temperatureEffect, 'b*');
xlabel("Temperature (°C)"); ylabel("cT (-)")


ct = temperatureEffect(2)  / temperatureEffect(1) ; 


%% Load data and initial parameters
pars.j_CO2m =  0.0075; 
pars.k_I = 0.075; 
pars.j_ECAm = 0.282; 
pars.j_ENAm =  0.00015; 
pars.k_EC = 0.02; 
pars.k_EN = 0.04; 
pars.j_ECM = 1e-06; 
pars.j_ENM = 4e-06;

j_CO2m_CT = pars.j_CO2m * ct; % mol CO2 mol V-1 h-1, max volume specific CO2 uptake rate
k_I_CT = pars.k_I * ct; % mol gamma mol PSU-1 h-1, dissociationrate of photosynthetic products
j_ECAm_CT = pars.j_ECAm * ct; %mol C mol V-1 h-1, max volume specific carbon assimilation
j_ENAm_CT = pars.j_ENAm * ct; %mol N mol V-1 h-1, max volume specific nitrogen assimilation
k_EC_CT = pars.k_EC * ct; %h-1, C reserve turnover
k_EN_CT = pars.k_EN* ct; %h-1, N reserve turnover
j_ECM_CT = pars.j_ECM * ct; %mol C mol V-1 h-1, volume specific maitenance rate paid by C reserve
j_ENM_CT = pars.j_ENM * ct; %mol N mol V-1 h-1, volume specific maitenance rate paid by N reserve
% 
%%
%% LAVAUD %%

pars.T_A = 4323; 
pars.T_AH = 68642; 
pars.T_AL = 60306; 
pars.T_H = 296.5; 
pars.T_L = 278.8; 
pars.T_ref = 293.15; 


% T = -5:0.05:35; %°C, temperatures tested, 2.6 to have maximum 12 points for estimation procedure of DEBtoolM
T = [20 25]; 
% pars_T = [T_A; T_L; T_H; T_AL; T_AH]'; % K
temperatureEffect =  zeros(length(T),1);

for l = 1:length(T)
        temp = T(l); %°C
        ct = Arrhenius_equation(pars, temp, 'simple');
        temperatureEffect(l,1) = ct; % -

end


% Create plot of function cT 
plot(T,temperatureEffect, 'b*');
xlabel("Temperature (°C)"); ylabel("cT (-)")


ct = temperatureEffect(2)  / temperatureEffect(1) ;

pars.j_CO2m =  0.0060; 
pars.k_I = 0.15; 
pars.j_ECAm = 0.006; 
pars.j_ENAm =  0.0005; 
pars.k_EC = 0.02; 
pars.k_EN = 0.04; 
pars.j_ECM = 1e-06; 
pars.j_ENM = 6.7e-04;

j_CO2m_CT = pars.j_CO2m * ct; % mol CO2 mol V-1 h-1, max volume specific CO2 uptake rate
k_I_CT = pars.k_I * ct; % mol gamma mol PSU-1 h-1, dissociationrate of photosynthetic products
j_ECAm_CT = pars.j_ECAm * ct; %mol C mol V-1 h-1, max volume specific carbon assimilation
j_ENAm_CT = pars.j_ENAm * ct; %mol N mol V-1 h-1, max volume specific nitrogen assimilation
k_EC_CT = pars.k_EC * ct; %h-1, C reserve turnover
k_EN_CT = pars.k_EN* ct; %h-1, N reserve turnover
j_ECM_CT = pars.j_ECM * ct; %mol C mol V-1 h-1, volume specific maitenance rate paid by C reserve
j_ENM_CT = pars.j_ENM * ct; %mol N mol V-1 h-1, volume specific maitenance rate paid by N reserve

%% GROSSOWICZ %%

%% LAVAUD %%

pars.T_A = 7964.5000; 
pars.T_ref = 293.15; 


% T = -5:0.05:35; %°C, temperatures tested, 2.6 to have maximum 12 points for estimation procedure of DEBtoolM
T = [20 25]; 
% pars_T = [T_A; T_L; T_H; T_AL; T_AH]'; % K
temperatureEffect =  zeros(length(T),1);

for l = 1:length(T)
        temp = T(l); %°C
        ct = Arrhenius_equation(pars, temp, 'simple');
        temperatureEffect(l,1) = ct; % -

end


% Create plot of function cT 
plot(T,temperatureEffect, 'b*');
xlabel("Temperature (°C)"); ylabel("cT (-)")


ct = temperatureEffect(2)  / temperatureEffect(1) ;

pars.j_CO2m =  0.3225; 
pars.k_I = 0; 
pars.j_ECAm = 0.0698; 
pars.j_ENAm =  0.0054; 
pars.j_EPAm =  0.0081/24; 
pars.k_EC = 0.2167; 
pars.k_EN = 0.2167; 
pars.k_EP = 0.2167; 
pars.j_ECM = 2.25e-3; 
pars.j_ENM = 5e-4;
pars.j_EPM = 0.00075/24;

j_CO2m_CT = pars.j_CO2m * ct; % mol CO2 mol V-1 h-1, max volume specific CO2 uptake rate
k_I_CT = pars.k_I * ct; % mol gamma mol PSU-1 h-1, dissociationrate of photosynthetic products

j_ECAm_CT = pars.j_ECAm * ct; %mol C mol V-1 h-1, max volume specific carbon assimilation
j_ENAm_CT = pars.j_ENAm * ct; %mol N mol V-1 h-1, max volume specific nitrogen assimilation
j_EPAm_CT = pars.j_EPAm * ct; %mol N mol V-1 h-1, max volume specific nitrogen assimilation

k_EC_CT = pars.k_EC * ct; %h-1, C reserve turnover
k_EN_CT = pars.k_EN * ct; %h-1, N reserve turnover
k_EP_CT = pars.k_EP * ct; %h-1, N reserve turnover

j_ECM_CT = pars.j_ECM * ct; %mol C mol V-1 h-1, volume specific maitenance rate paid by C reserve
j_ENM_CT = pars.j_ENM * ct; %mol N mol V-1 h-1, volume specific maitenance rate paid by N reserve
j_EPM_CT = pars.j_EPM * ct; %mol N mol V-1 h-1, volume specific maitenance rate paid by N reserve