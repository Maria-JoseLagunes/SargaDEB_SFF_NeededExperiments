%From estimated at 25°C to 20°C for Table in paper
%% Estimated parameters %%
% Transform estimated parametres at 25°C to 20°c T_ref

pars.T_A = 6603; 
pars.T_AH = 58010; 
pars.T_AL = 54620; 
pars.T_H = 304.3; 
pars.T_L = 288.9; 
pars.T_ref = 298.15; 


% T = -5:0.05:35; %°C, temperatures tested, 2.6 to have maximum 12 points for estimation procedure of DEBtoolM
T = [25 20]; 
% pars_T = [T_A; T_L; T_H; T_AL; T_AH]'; % K
temperatureEffect =  zeros(length(T),1);

for l = 1:length(T)
        temp = T(l); %°C
        % ct = Arrhenius_equation(pars, temp, 'simple');
        ct = tempcorr (C2K(temp), pars.T_ref, pars.T_A);
        temperatureEffect(l,1) = ct; % -

end


% Create plot of function cT 
plot(T,temperatureEffect, 'b*');
xlabel("Temperature (°C)"); ylabel("cT (-)")


ct = temperatureEffect(2)  / temperatureEffect(1) ; 


%% Load data and initial parameters
pars.j_CO2m =   	0.01076; 
pars.k_I = 0.1047; 
pars.j_ECAm = 0.4047 ;
pars.j_ENAm = 0.0001228 	; 
pars.k_EC =  	0.02679; 
pars.k_EN =  	0.3132; 
pars.j_ECM = 1.414e-06; 
pars.j_ENM = 3.871e-06;

j_CO2m_CT = pars.j_CO2m * ct; % mol CO2 mol V-1 h-1, max volume specific CO2 uptake rate
k_I_CT = pars.k_I * ct; % mol gamma mol PSU-1 h-1, dissociationrate of photosynthetic products
j_ECAm_CT = pars.j_ECAm * ct; %mol C mol V-1 h-1, max volume specific carbon assimilation
j_ENAm_CT = pars.j_ENAm * ct; %mol N mol V-1 h-1, max volume specific nitrogen assimilation
k_EC_CT = pars.k_EC * ct; %h-1, C reserve turnover
k_EN_CT = pars.k_EN* ct; %h-1, N reserve turnover
j_ECM_CT = pars.j_ECM * ct; %mol C mol V-1 h-1, volume specific maitenance rate paid by C reserve
j_ENM_CT = pars.j_ENM * ct; %mol N mol V-1 h-1, volume specific maitenance rate paid by N reserve
% 