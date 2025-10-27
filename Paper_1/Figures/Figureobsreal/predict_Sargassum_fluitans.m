function [prdData, info] = predict_Sargassum_fluitans(pars, data, auxData)
global infosgr2



vars_pull(data);
vars_pull(auxData);

pars_T = [pars.T_A; pars.T_L; pars.T_H; pars.T_AL; pars.T_AH]; % K
 
TC_t_Ww_28 = tempcorr(temp.tWw_MG2023a_28, pars.T_ref, pars_T);  % -



%Transform parameters with z factor algal
dwratio_SF = pars.dwratio * pars.zalgal ; % gWd/gWw, Dry to wet (gdw/gww) weight ratio
smcoeff_SF = pars.smcoeff / pars.zalgal;  % g cm-2, Coefficient of linear regression from Chambon

% MG2023a
% time and total wet weight
Wd_0_t_Ww_28_ash= Ww0.tWw_MG2023a_28 * dwratio_SF;  %g dW, dry weight
Wd_0_t_Ww_28 = Wd_0_t_Ww_28_ash * (1 - pars.x_moist - pars.x_ash); 

M_V_0_t_Ww_28 = Wd_0_t_Ww_28  / (pars. w_V +   (pars.m_EN_0 * pars.w_EN) +   (pars.m_EC_0* pars.w_EC)); % mol V, structural initial mass 
CI_t_Ww_28 = [pars.m_EC_0, pars.m_EN_0, M_V_0_t_Ww_28];
t = [0 ; (tWw_MG2023a_28(2:end,1) * 24)]; %in hours


forcingVariables.CO2 = CO2.tWw_MG2023a_28; 
forcingVariables.irradiance = light.tWw_MG2023a_28; 
forcingVariables.N = nitrogen.tWw_MG2023a_28; 
forcingVariables.P =  phosphorous.tWw_MG2023a_28; 
forcingVariables.photoPeriod =  photoPeriod.tWw_MG2023a_28; 
forcingVariables.temp =  temp.tWw_MG2023a_28; 

forcingVariables.dt = linspace(t(1), t(end), length(temp.tWw_MG2023a_28));




[~, mECENV_t_Ww_28] = ode89(@(t,mECENV_t_Ww_28) get_mECENV(t,  mECENV_t_Ww_28, pars, TC_t_Ww_28, forcingVariables), t, CI_t_Ww_28);
% mECENV_t_Ww_28(1,:) = [] ;



m_EC = mECENV_t_Ww_28(:,1); % mol C mol V-1
m_EN = mECENV_t_Ww_28(:,2); % mol N mol V-1
M_V = mECENV_t_Ww_28(:,3);  % mol V

% Wd_MG2023a_28 = (pars.w_V +   (m_EN  * pars.w_EN) +   (m_EC * pars.w_EC)) .* M_V; %g dW, dry weight

Wd_MG2023a_28_free= (pars.w_V + pars.w_EN * m_EN + pars.w_EC * m_EC) .* M_V;
Wd_MG2023a_28 = Wd_MG2023a_28_free / (1 - pars.x_moist - pars.x_ash); 
Ww_MG2023a_28 = Wd_MG2023a_28 ./ dwratio_SF   ; %g wW, wet weight


%Saving flux as datatable
colName = string({"ct", "j_CO2", "j_I", "j_EC_A", "j_EN_A", "j_EC_C", "j_EN_C", "j_EC_MC", "j_EN_MN", "j_V_MC", ...
                    "j_V_MN", "j_V_M", "j_EC_G", "j_EN_G", "j_VG", "j_EC_R", "j_EN_R", "r", "I", "CO2", "N","P", "T"});

J = zeros(length(t),length(colName)); %initializing flux of organic vector

for i = 1:length(t)
    [~,J(i,:)] = get_mECENV (t(i), mECENV_t_Ww_28(i,:),pars, TC_t_Ww_28, forcingVariables); 
end


% Create a structure
J_struct= struct();
% 
% Assign each column of J to a field in the structure
for i = 1:length(colName)
    J_struct.(colName(i)) = J(:, i);
end



prdDatatest.tWw_MG2023a_28 = Ww_MG2023a_28; 
prdDatatest.Irradiance = J_struct.I; 
prdDatatest.Nitrogen = J_struct.N;
prdDatatest.Phosphorus = J_struct.P;
prdDatatest.Temperature = K2C(J_struct.T);

prdDatatest.m_EC = m_EC;
prdDatatest.m_EN = m_EN;
prdDatatest.M_V = M_V;



% prdData.tWw_MG2023a_28 = Ww_MG2023a_28;
prdData  = prdDatatest;

info = 1;  
end


