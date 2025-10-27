% ---------------------------------------------------------------------
%% Creation of experiments of nutrient uptake, NO3
% Maria José Lagunes
% First created : 2025/03/31
% ---------------------------------------------------------------------
function [simu, obs] = run_Nuptake(N_microConcentration, temp_C)
global pets

% temp_C = [22 26 28];
temp_K = C2K(temp_C); 
% N_microConcentration =  0:7:50 ; % micro mol NO3 L-1
N_molConcentration = N_microConcentration * 1e-6; 

pets = {'Sargassum_fluitans'};


pars_initnm = ['pars_init_', pets{1}];
resultsnm   = ['results_', pets{1}, '.mat'];
auxDatanm = ['auxData.', pets{1}];


[pars, ~, ~] = feval(pars_initnm, []);

pars_T = [pars.T_A; pars.T_L; pars.T_H; pars.T_AL; pars.T_AH]; % K

ct = tempcorr(temp_K, pars.T_ref, pars_T); %-


Ww0 = 2.5 ; %g Ww, initial wet biomass 
% Wd0 = Ww0 * pars.dwratio;  %g dW, initial dry weight biomass
% M_V_0 = Wd0 / (pars.w_V +   (pars. m_EN_0 * pars.w_EN) +   (pars.m_EC_0 * pars.w_EC)); % mol V, structural intial mass

Wd0_ash = Ww0 * pars.dwratio;  %g dW, initial dry weight biomass
Wd0 = Wd0_ash * (1 - pars.x_moist - pars.x_ash); 
M_V_0 = Wd0 / (pars.w_V +   (pars.m_EN_0 * pars.w_EN) +   (pars.m_EC_0 * pars.w_EC)); % mol V, structural intial mass
      


j_ENAm_unitTransformation = pars.j_ENAm * M_V_0 / Wd0 ; %mol N g dW-1 h-1,transformed from mol N mol V-1 h-1
j_ENAm_CT = j_ENAm_unitTransformation * ct;  % mol N g dW-1 h-1

j_EN_A = []; 

% Nitrogen into reserve
for i = 1:length(ct)
    simu(i).N_microConcentration = N_microConcentration;
    simu(i).temp = temp_C(i); 
    simu(i).Ww0 = Ww0; 
    simu(i).Wd0 = Wd0; 
    simu(i).M_V_0 = M_V_0; 
    for l = 1:length(N_molConcentration)
        j_EN_A(:,l) = j_ENAm_CT(i) * (N_molConcentration(l) / (N_molConcentration(l) + pars.K_N)); % mol N g dW-1 h-1
    end
    obs(i).nitrateUptake = j_EN_A;
    obs(i).j_EN_A_inmicroMol = j_EN_A / 1e-6; % micro mol N g dW-1 h-1
    plot(simu(i).N_microConcentration,obs(i).j_EN_A_inmicroMol, ...
        'LineWidth', 2); hold on; 

end

xlabel('N concentration (micro mol NO3 L-1)');
ylabel('N uptake rate (micro mol g dW-1 h-1');
% legend("22°C", "26°C", "28°C")

end