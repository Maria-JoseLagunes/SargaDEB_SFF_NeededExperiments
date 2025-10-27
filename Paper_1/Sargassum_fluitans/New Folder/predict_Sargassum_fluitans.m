function [prdData, info] = predict_Sargassum_fluitans(pars, data, auxData)

% %unpack parameters, data and auxData
% vars_pull(pars);
vars_pull(data);
vars_pull(auxData);

m_ECm = (pars.j_ECAm - pars.kap_EC * pars.j_ECM)/ (1 - pars.kap_EC)/ pars.k_EC; %Max C reserve density,  mol E_C mol V-1
m_ENm = (pars.j_ENAm - pars.kap_EN * pars.j_ENM)/ (1 - pars.kap_EN)/ pars.k_EN; %Max N reserve density, mol E_N mol V-1

RR_Nm = pars.n_NV + m_ENm; % (mol E_N  mol V-1  +  mol E_N mol V-1)


% Check temperature parameter values
if pars.T_L > pars.T_ref || pars.T_H < pars.T_ref
  info = 0; prdData = []; return;
end



% r = (0.1 / 24) * 0.95;  %in hours, rmax from studies
% 
% pars.m_EC_0 = (pars.j_ECAm - pars.kap_EC * pars.j_ECM - pars.kap_EC * pars.y_ECV * r)/ ...
%                 ((1 - pars.kap_EC) * pars.k_EC + pars.kap_EC  * r);
% 
% pars.m_EN_0 = (pars.j_ENAm - pars.kap_EN * pars.j_ENM - pars.kap_EN * pars.y_ENV * r)/ ...
%                 ((1 - pars.kap_EN) * pars.k_EN + pars.kap_EN * r);


% 
filterChecks =   pars.y_ENV < pars.n_NV   || ...     % mass conservation  n_HV > 263/106 || n_OV > 10/106 || n_NV > 16/106 || n_PV > 1/106 || ...     % assuming there is some reserve in Redfield ratio
                 pars.j_ECAm < pars.j_ECM || pars.j_ENAm < pars.j_ENM  || ...                 % species survival
                 pars.j_ECAm < 1.1 * pars.j_ECM || pars.j_ENAm < 1.1 * pars.j_ENM || ...              % species survival
                 RR_Nm > 1200|| ... % Ratios check up
                 pars.kap_EC < 0 || pars.kap_EC > 1 || pars.kap_EN < 0 || pars.kap_EN > 1 || ...  % energy conservation
                 pars.j_ECM < 1e-8||  pars.j_ENM < 1e-8 || ...                  
                 pars.m_EC_0 < 0 ||  pars.m_EN_0 < 0 ; 
                 
 
             
  if filterChecks  
    info = 0;
    prdData = [];
    return;
  end  

%compute temperature correction factors
pars_T = [pars.T_A; pars.T_L; pars.T_H; pars.T_AL; pars.T_AH]; % K
 
TC_t_Ww_28 = tempcorr(temp.tWw_MG2023a_28, pars.T_ref, pars_T);  % -
TC_t_Ww_31 = tempcorr(temp.tWw_MG2023a_31, pars.T_ref, pars_T);  %-

TC_t_Ww_b_22 = tempcorr(temp.tWw_MG2023b_22, pars.T_ref, pars_T);  %-
TC_t_Ww_b_25 = tempcorr(temp.tWw_MG2023b_25, pars.T_ref, pars_T);  %-
TC_t_Ww_b_28 = tempcorr(temp.tWw_MG2023b_28, pars.T_ref, pars_T);  %-
TC_t_Ww_b_31 = tempcorr(temp.tWw_MG2023b_31, pars.T_ref, pars_T);  %-



TC_HL = tempcorr(temp.I_J_O2_VE2023_HL, pars.T_ref,pars_T); % - 
% TC_LL = tempcorr(temp.I_J_O2_VE2023_LL, pars.T_ref,pars_T); % -

% TC_tWw_VE2023_LL = tempcorr(temp.tWw_VE2023_LL, pars.T_ref, pars_T);
% TC_tWw_VE2023_ML = tempcorr(temp.tWw_VE2023_ML, pars.T_ref, pars_T);
% TC_tWw_VE2023_HL = tempcorr(temp.tWw_VE2023_HL, pars.T_ref, pars_T);


%Correction of parameter values for Vazquez-Elizondo data
k_I_CT_HL = pars.k_I * TC_HL; % mol γ mol PSU–1 h–1
% k_I_CT_LL = pars.k_I * TC_LL; % mol γ mol PSU–1 h–



%Transform parameters with z factor algal
dwratio_SF = pars.dwratio * pars.zalgal ; % gWd/gWw, Dry to wet (gdw/gww) weight ratio
smcoeff_SF = pars.smcoeff / pars.zalgal;  % g cm-2, Coefficient of linear regression from Chambon

% zero-variate data
% HT2023
TC_HT2023 = tempcorr(temp.HT2023_CNratio, pars.T_ref, pars_T);
% 
% rmax = 0.1; 
% TC_pars = [pars.j_ECAm, pars.j_ECM, pars.k_EC, pars.j_ENAm, pars.j_ENM, pars.k_EN];
% pars_m0 = [TC_pars, pars.kap_EC, pars.y_ECV ]; 
% [m_EC_0, m_EN0] = get_mECEN0_CT(rmax, pars_m0), 


Wd_0 = Ww0.HT2023_CNratio * dwratio_SF;  %g dW, dry weight
M_V_0 = Wd_0  / (pars. w_V +   (pars.m_EN_0 * pars.w_EN) +   (pars.m_EC_0* pars.w_EC)); % mol V, structural initial mass 

CI_HT= [pars.m_EC_0, pars.m_EN_0, M_V_0];
t = [0 time.HT2023_CNratio*24]; % in hours


forcingVariables.CO2 = CO2.HT2023_CNratio; 
forcingVariables.lightIntensity = light.HT2023_CNratio; 
forcingVariables.N = nitrogen.HT2023_CNratio; 
forcingVariables.P =  phosphorous.HT2023_CNratio; 
forcingVariables.photoPeriod =  photoPeriod.HT2023_CNratio; 


[~, mECENV] = ode89(@(t,mECENV) get_mECENV(t,  mECENV, pars, TC_HT2023, forcingVariables), t, CI_HT);



m_EC = mECENV(:,1); % mol C mol V-1
m_EN = mECENV(:,2); % mol N mol V-1
M_V = mECENV(:,3);  % mol V

% Ratios
M_EC = m_EC .* M_V; % mol C
M_EN = m_EN .* M_V; % mol N



C_total = (pars.n_CV * M_V) + (pars.n_CEC * M_EC); % mol C  
N_total = (pars.n_NV * M_V) + (pars.n_NEN * M_EN); % mol N
P_total  = (pars.n_PV * M_V); % mol P


Wd_HT2023 = (pars.w_V + pars.w_EN * m_EN + pars.w_EC * m_EC) .* M_V;
 
% To verify order of magnitud
percC_gdW =  (C_total * pars.w_C  ./ Wd_HT2023) * 100; % %C (gdW)
percN_gdW =  (N_total  * pars.w_N ./ Wd_HT2023) * 100; % %N (gdW)
percP_gdW  = (P_total * pars.w_P  ./ Wd_HT2023) * 100; % %P (gdW)

percC_gdW_end = percC_gdW(end); %To verify order of magnitud
percN_gdW_end = percN_gdW(end); %To verify order of magnitud
percP_gdW_end = percP_gdW(end); %To verify order of magnitud

CNtotalRatio = C_total ./ N_total ;  % mol C mol N-1
CPtotalRatio = C_total ./ P_total ;  % mol C mol P-1
NPtotalRatio = N_total ./ P_total ;  % mol N mol P-1

CNendRatio = CNtotalRatio(end); % mol C mol N-1
CPendRatio = CPtotalRatio(end); % mol C mol P-1
NPendRatio = NPtotalRatio(end); % mol N mol P-1

prdData.HT2023_CNratio = CNendRatio;
prdData.HT2023_CPratio = CPendRatio;
prdData.HT2023_NPratio = NPendRatio;

% uni-variate data

% HS1987 !!!! SNATANS 

temps = [12, 18, 24, 30];
TC_HS1987_12 = tempcorr(temp.tWw_HS1987_12, pars.T_ref, pars_T);  %-
TC_HS1987_18 = tempcorr(temp.tWw_HS1987_18, pars.T_ref, pars_T);  %-
TC_HS1987_24 = tempcorr(temp.tWw_HS1987_24, pars.T_ref, pars_T);  %-
TC_HS1987_30 = tempcorr(temp.tWw_HS1987_30, pars.T_ref, pars_T);  %-

for i = 1:length(temps)
    T = temps(i);
    
    % Construct variable names dynamically
    tWw_var = eval(['tWw_HS1987_', num2str(T)]);
    CO2_var = eval(['CO2.tWw_HS1987_', num2str(T)]);
    light_var = eval(['light.tWw_HS1987_', num2str(T)]);
    N_var = eval(['nitrogen.tWw_HS1987_', num2str(T)]);
    P_var = eval(['phosphorous.tWw_HS1987_', num2str(T)]);
    photo_var = eval(['photoPeriod.tWw_HS1987_', num2str(T)]);
    TC_var = eval(['TC_HS1987_', num2str(T)]);
    Ww0_var = eval(['Ww0.tWw_HS1987_', num2str(T)]);
    
    % Initial conditions
    Wd_0 = Ww0_var * dwratio_SF;  % dry weight
    M_V_0 = Wd_0 / (pars.w_V + (pars.m_EN_0 * pars.w_EN) + (pars.m_EC_0 * pars.w_EC));
    CI = [pars.m_EC_0, pars.m_EN_0, M_V_0];
    
    % Time vector in hours
    t = [0; tWw_var(2:end,1) * 24]; % ensure t(1) = 0

    % Forcings
    forcingVariables.CO2 = CO2_var;
    forcingVariables.lightIntensity = light_var;
    forcingVariables.N = N_var;
    forcingVariables.P = P_var;
    forcingVariables.photoPeriod = photo_var;

    nt = length(t);
    
  
    [~, mECENV] = ode89(@(t, mECENV) get_mECENV(t, mECENV, pars, TC_var, forcingVariables), t, CI);

    if nt == 2
        mECENV = mECENV([1 end], :); 
    end

    m_EC = mECENV(:,1);
    m_EN = mECENV(:,2);
    M_V  = mECENV(:,3);

    Wd = (pars.w_V + pars.w_EN * m_EN + pars.w_EC * m_EC) .* M_V;
    Ww = Wd ./ dwratio_SF;

    % Store in prdData
    prdData.(['tWw_HS1987_', num2str(T)]) = Ww;
end


% MG2023a
% time and total wet weight
Wd_0_t_Ww_28 = Ww0.tWw_MG2023a_28 * dwratio_SF;  %g dW, dry weight
M_V_0_t_Ww_28 = Wd_0_t_Ww_28  / (pars. w_V +   (pars.m_EN_0 * pars.w_EN) +   (pars.m_EC_0* pars.w_EC)); % mol V, structural initial mass 
CI_t_Ww_28 = [pars.m_EC_0, pars.m_EN_0, M_V_0_t_Ww_28];
t = [0 ; (tWw_MG2023a_28(2:end,1) * 24)]; %in hours

forcingVariables.CO2 = CO2.tWw_MG2023a_28; 
forcingVariables.lightIntensity = light.tWw_MG2023a_28; 
forcingVariables.N = nitrogen.tWw_MG2023a_28; 
forcingVariables.P =  phosphorous.tWw_MG2023a_28; 
forcingVariables.photoPeriod =  photoPeriod.tWw_MG2023a_28; 


[~, mECENV_t_Ww_28] = ode89(@(t,mECENV_t_Ww_28) get_mECENV(t,  mECENV_t_Ww_28, pars, TC_t_Ww_28, forcingVariables), t, CI_t_Ww_28);
% mECENV_t_Ww_28(1,:) = [] ;


m_EC = mECENV_t_Ww_28(:,1); % mol C mol V-1
m_EN = mECENV_t_Ww_28(:,2); % mol N mol V-1
M_V = mECENV_t_Ww_28(:,3);  % mol V

Wd_MG2023a_28 = (pars.w_V +   (m_EN  * pars.w_EN) +   (m_EC * pars.w_EC)) .* M_V; %g dW, dry weight


Ww_MG2023a_28 = Wd_MG2023a_28 ./ dwratio_SF   ; %g wW, wet weight


prdData.tWw_MG2023a_28 = Ww_MG2023a_28;

%%

Wd_0_t_Ww_31 = Ww0.tWw_MG2023a_31 * dwratio_SF; %g dW, dry weight
M_V_0_t_Ww_31 = Wd_0_t_Ww_31  / (pars.w_V +   (pars.m_EN_0 * pars.w_EN) +   (pars.m_EC_0* pars.w_EC)); % mol V, structural initial mass 
CI_t_Ww_31 = [pars.m_EC_0, pars.m_EN_0, M_V_0_t_Ww_31];
t = [0 (tWw_MG2023a_31(2:end,1) * 24)']; %in hours

forcingVariables.CO2 = CO2.tWw_MG2023a_31; 
forcingVariables.lightIntensity = light.tWw_MG2023a_31; 
forcingVariables.N = nitrogen.tWw_MG2023a_31; 
forcingVariables.P =  phosphorous.tWw_MG2023a_31; 
forcingVariables.photoPeriod =  photoPeriod.tWw_MG2023a_31; 

[~, mECENV_t_Ww_31] = ode89(@(t,mECENV_t_Ww_31) get_mECENV(t, mECENV_t_Ww_31, pars, TC_t_Ww_31, forcingVariables), t, CI_t_Ww_31);
% mECENV_t_Ww_31(1,:) = [] ;


m_EC = mECENV_t_Ww_31(:,1); % mol C mol V-1
m_EN = mECENV_t_Ww_31(:,2); % mol N mol V-1
M_V = mECENV_t_Ww_31(:,3);  % mol V

Wd_MG2023a_31 = (pars.w_V +   (m_EN * pars.w_EN) +   (m_EC * pars.w_EC)) .* M_V; %g dW, dry weight
Ww_MG2023a_31 = Wd_MG2023a_31 ./ dwratio_SF ; %g wW, wet weight


prdData.tWw_MG2023a_31 = Ww_MG2023a_31;


% MG2023b
% time and total weight
% 22 
Wd_0_t_Ww = Ww0.tWw_MG2023b_22 * dwratio_SF;  %g dW, dry weight
M_V_0_t_Ww= Wd_0_t_Ww  / (pars.w_V +   (pars.m_EN_0 * pars.w_EN) +   (pars.m_EC_0* pars.w_EC)); % mol V, structural initial mass 
CI_t_Ww = [pars.m_EC_0, pars.m_EN_0, M_V_0_t_Ww];
t = [0 (tWw_MG2023b_22(2:end,1) * 24)']; %in hours

forcingVariables.CO2 = CO2.tWw_MG2023b_22; 
forcingVariables.lightIntensity = light.tWw_MG2023b_22; 
forcingVariables.N = nitrogen.tWw_MG2023b_22; 
forcingVariables.P =  phosphorous.tWw_MG2023b_22; 
forcingVariables.photoPeriod =  photoPeriod.tWw_MG2023b_22; 

nt = length(t);

[~, mECENV_t_Ww_22] = ode89(@(t,mECENV_t_Ww_22) get_mECENV(t,  mECENV_t_Ww_22, pars, TC_t_Ww_b_22, forcingVariables), t, CI_t_Ww);


if nt == 2
    mECENV_t_Ww_22 = mECENV_t_Ww_22([1 end],:);
end





m_EC = mECENV_t_Ww_22(:,1); % mol C mol V-1
m_EN = mECENV_t_Ww_22(:,2); % mol N mol V-1
M_V = mECENV_t_Ww_22(:,3);  % mol V

Wd_MG2023b_22 = (pars.w_V +   (m_EN * pars.w_EN) +   (m_EC * pars.w_EC)) .* M_V; %g dW, dry weight
Ww_MG2023b_22 = Wd_MG2023b_22 ./ dwratio_SF   ; %g wW, wet weight

% 25
Wd_0_t_Ww = Ww0.tWw_MG2023b_25 * dwratio_SF;  %g dW, dry weight
M_V_0_t_Ww= Wd_0_t_Ww  / (pars.w_V +   (pars.m_EN_0 * pars.w_EN) +   (pars.m_EC_0 * pars.w_EC)); % mol V, structural initial mass 
CI_t_Ww = [pars.m_EC_0, pars.m_EN_0, M_V_0_t_Ww];
t = [0 (tWw_MG2023b_25(2:end,1) * 24)']; %in hours

forcingVariables.CO2 = CO2.tWw_MG2023b_25; 
forcingVariables.lightIntensity = light.tWw_MG2023b_25; 
forcingVariables.N = nitrogen.tWw_MG2023b_25; 
forcingVariables.P =  phosphorous.tWw_MG2023b_25; 
forcingVariables.photoPeriod =  photoPeriod.tWw_MG2023b_25; 

nt = length(t);


[~, mECENV_t_Ww_25] = ode89(@(t,mECENV_t_Ww_25) get_mECENV(t,  mECENV_t_Ww_25, pars, TC_t_Ww_b_25, forcingVariables), t, CI_t_Ww);
% mECENV_t_Ww_25(1,:) = [] ;

if nt == 2
    mECENV_t_Ww_25 = mECENV_t_Ww_25([1 end],:);
end





m_EC = mECENV_t_Ww_25(:,1); % mol C mol V-1
m_EN = mECENV_t_Ww_25(:,2); % mol N mol V-1
M_V = mECENV_t_Ww_25(:,3);  % mol V

Wd_MG2023b_25 = (pars.w_V +   (m_EN * pars.w_EN) +   (m_EC * pars.w_EC)) .* M_V; %g dW, dry weight
Ww_MG2023b_25 = Wd_MG2023b_25 ./ dwratio_SF   ; %g wW, wet weight

% 28
Wd_0_t_Ww = Ww0.tWw_MG2023b_28* dwratio_SF;  %g dW, dry weight
M_V_0_t_Ww = Wd_0_t_Ww  / (pars.w_V +   (pars.m_EN_0 * pars.w_EN) +   (pars.m_EC_0* pars.w_EC)); % mol V, structural initial mass 
CI_t_Ww = [pars.m_EC_0, pars.m_EN_0, M_V_0_t_Ww];
t = [0 (tWw_MG2023b_28(2:end,1) * 24)']; %in hours

forcingVariables.CO2 = CO2.tWw_MG2023b_28; 
forcingVariables.lightIntensity = light.tWw_MG2023b_28; 
forcingVariables.N = nitrogen.tWw_MG2023b_28; 
forcingVariables.P =  phosphorous.tWw_MG2023b_28; 
forcingVariables.photoPeriod =  photoPeriod.tWw_MG2023b_28; 

nt = length(t);


[~, mECENV_t_Ww_28] = ode89(@(t,mECENV_t_Ww_28) get_mECENV(t,  mECENV_t_Ww_28, pars, TC_t_Ww_b_28, forcingVariables), t, CI_t_Ww);
% mECENV_t_Ww_28(1,:) = [] ;


if nt == 2
    mECENV_t_Ww_28 = mECENV_t_Ww_28([1 end],:);
end


m_EC = mECENV_t_Ww_28(:,1); % mol C mol V-1
m_EN = mECENV_t_Ww_28(:,2); % mol N mol V-1
M_V = mECENV_t_Ww_28(:,3);  % mol V

Wd_MG2023b_28 = (pars.w_V +   (m_EN * pars.w_EN) +   (m_EC * pars.w_EC)) .* M_V; %g dW, dry weight
Ww_MG2023b_28 = Wd_MG2023b_28 ./ dwratio_SF   ; %g wW, wet weight


% 31 
Wd_0_t_Ww = Ww0.tWw_MG2023b_31 * dwratio_SF;  %g dW, dry weight
M_V_0_t_Ww = Wd_0_t_Ww  / (pars.w_V +   (pars.m_EN_0 * pars.w_EN) +   (pars.m_EC_0* pars.w_EC)); % mol V, structural initial mass 
CI_t_Ww = [pars.m_EC_0, pars.m_EN_0, M_V_0_t_Ww];
t = [0 (tWw_MG2023b_31(2:end,1) * 24)']; %in hours

forcingVariables.CO2 = CO2.tWw_MG2023b_31; 
forcingVariables.lightIntensity = light.tWw_MG2023b_31; 
forcingVariables.N = nitrogen.tWw_MG2023b_31; 
forcingVariables.P =  phosphorous.tWw_MG2023b_31; 
forcingVariables.photoPeriod =  photoPeriod.tWw_MG2023b_31; 

nt = length(t);


[~, mECENV_t_Ww_31] = ode89(@(t,mECENV_t_Ww_31) get_mECENV(t,  mECENV_t_Ww_31, pars, TC_t_Ww_b_31, forcingVariables), t, CI_t_Ww);
% mECENV_t_Ww_31(1,:) = [] ;
if nt == 2
    mECENV_t_Ww_31 = mECENV_t_Ww_31([1 end],:);
end



m_EC = mECENV_t_Ww_31(:,1); % mol C mol V-1
m_EN = mECENV_t_Ww_31(:,2); % mol N mol V-1
M_V = mECENV_t_Ww_31(:,3);  % mol V

Wd_MG2023b_31 = (pars.w_V +   (m_EN * pars.w_EN) +   (m_EC * pars.w_EC)) .* M_V; %g dW, dry weight
Ww_MG2023b_31 = Wd_MG2023b_31 ./ dwratio_SF   ; %g wW, wet weight

%Pack output
prdData.tWw_MG2023b_22 = Ww_MG2023b_22;
prdData.tWw_MG2023b_25 = Ww_MG2023b_25;
prdData.tWw_MG2023b_28 = Ww_MG2023b_28;
prdData.tWw_MG2023b_31 = Ww_MG2023b_31;



% VE2023
% Photosynthesis - Specific relaxation rate

I_J_O2_VE2023_HL(:,1) =  I_J_O2_VE2023_HL(:,1) * 1e-6 * 3600; % convert light intensity from micro mol to mol and from seconds to hours
Wd = surface.I_J_O2_VE2023_HL * smcoeff_SF;  % g dW SM = g cm-2 1.5 in cm-2 from VE2023, same for HL and LL
M_V= Wd / (pars.w_V +   (pars.m_EN_0 * pars.w_EN) +   (pars.m_EC_0* pars.w_EC)); %mol V


J_I_HL = (pars.rho_PSU * I_J_O2_VE2023_HL(:,1) * pars.eps_I) ./ ...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
    (1 +  ((I_J_O2_VE2023_HL(:,1) * pars.eps_I) / k_I_CT_HL)) ; %% mol gamma mol V-1 h-1

J_O2_HL =  J_I_HL .* (M_V / Wd)  * pars.y_LO2 * smcoeff_SF * dwratio_SF; % mol 02 cm-2  h-1
% J_O2_HL =  J_I_HL .* (M_V / Wd)  * pars.y_LO2 ; % mol 02 gdW-1  h-1

J_O2_HL_micromol= J_O2_HL / 1e-6; % micro mol O2 cm-2 h-1
prdData.I_J_O2_VE2023_HL= J_O2_HL_micromol; 


% rho_PSU_LL = pars.rho_PSU * 0.5; 
% I_J_O2_VE2023_LL(:,1) =  I_J_O2_VE2023_LL(:,1) * 1e-6 * 3600; % convert light intensity from micro mol to mol and from seconds to hours
% 
% J_I_LL = (rho_PSU_LL* I_J_O2_VE2023_LL(:,1) * pars.eps_I) ./ ...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
%     (1 +  ((I_J_O2_VE2023_LL(:,1) * pars.eps_I) / k_I_CT_LL)) ; %% mol gamma mol V-1 h-1
% 
% % J_I_LL = (pars.rho_PSU * I_J_O2_VE2023_LL(:,1) * pars.eps_I) ./ ...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
%     % (1 +  ((I_J_O2_VE2023_LL(:,1) * pars.eps_I) / k_I_CT_LL)) ; %% mol gamma mol V-1 h-1
% 
% J_O2_LL =  J_I_LL .* (M_V / Wd)  * pars.y_LO2 * smcoeff_SF * dwratio_SF; % mol 02 cm-2  h-1
% J_O2_LL =  J_I_LL .* (M_V / Wd)  * pars.y_LO2; % mol 02 cm-2  h-1
% 
% J_O2_LL_micromol= J_O2_LL / 1e-6; %micro mol O2 cm-2 h-1
% 
% prdData.I_J_O2_VE2023_LL = J_O2_LL_micromol; 


% Irradiance condition - GR - CN content
% LL
% Wd_0_tWw_VE2023_LL = Ww0.tWw_VE2023_LL * dwratio_SF;  %g dW, dry weight
% M_V_0_tWw_VE2023_LL = Wd_0_tWw_VE2023_LL  / (pars.w_V+   (pars.m_EN_0* pars.w_EN) +   (pars.m_EC_0* pars.w_EC)); % mol V, structural initial mass 
% CI_tWw_VE2023_LL = [pars.m_EC_0, pars.m_EN_0, M_V_0_tWw_VE2023_LL];
% t = [1 (tWw_VE2023_LL(2:end,1) * 24)']; %in hours
% 
% forcingVariables.CO2 = CO2.tWw_VE2023_LL; 
% forcingVariables.lightIntensity = light.tWw_VE2023_LL; 
% forcingVariables.N = nitrogen.tWw_VE2023_LL; 
% forcingVariables.P =  phosphorous.tWw_VE2023_LL; 
% forcingVariables.photoPeriod =  photoPeriod.tWw_VE2023_LL; 
% 
% [~, mECENV_tWw_VE2023_LL] = ode89(@(t,mECENV_tWw_VE2023_LL) get_mECENV(t,  mECENV_tWw_VE2023_LL, pars, TC_tWw_VE2023_LL, forcingVariables), t, CI_tWw_VE2023_LL);
% % mECENV_Ww_28(1,:) = [] ;
% 
% 
% m_EC = mECENV_tWw_VE2023_LL(:,1); % mol C mol V-1
% m_EN = mECENV_tWw_VE2023_LL(:,2); % mol N mol V-1
% M_V = mECENV_tWw_VE2023_LL(:,3);  % mol V
% 
% Wd_tWw_VE2023_LL = (pars.w_V + (m_EN .* pars.w_EN) + (m_EC .* pars.w_EC)) .* M_V; %g dW, dry weight
% Ww_tWw_VE2023_LL = Wd_tWw_VE2023_LL ./ dwratio_SF   ; %g wW, wet weight
% 
% 
% % N:C ratio during simulation and total percentages of C and N
% M_EC = m_EC(end) .* M_V(end); % mol C
% M_EN = m_EN(end) .* M_V(end); % mol N
% 
% C_total = (pars.n_CV * M_V(end)) + (pars.n_CEC * M_EC); % mol C  
% N_total = (pars.n_NV * M_V(end)) + (pars.n_NEN * M_EN); % mol N
% % P_total  = (pars.n_PV * M_V(end)); % mol P
% 
% percC_gdW = (C_total * pars.w_C ./ Wd_tWw_VE2023_LL(end) ) * 100; % %C (gdW)
% percN_gdW = (N_total  * pars.w_N ./ Wd_tWw_VE2023_LL(end) ) * 100; % %N (gdW)
% % percP_gdW  = (P_total * pars.w_P./ Wd_tWw_VE2023_LL(end) ) * 100; % %P (gdW)
% 
% 
% 
% prdData.tWw_VE2023_LL = Ww_tWw_VE2023_LL;
% prdData.C_gdW_VE2023_LL = percC_gdW;
% prdData.N_gdW_VE2023_LL = percN_gdW;
% 
% % ML
% Wd_0_tWw_VE2023_ML = Ww0.tWw_VE2023_ML * dwratio_SF;  %g dW, dry weight
% M_V_0_tWw_VE2023_ML = Wd_0_tWw_VE2023_ML  / (pars.w_V+   (pars.m_EN_0* pars.w_EN) +   (pars.m_EC_0* pars.w_EC)); % mol V, structural initial mass 
% CI_tWw_VE2023_ML = [pars.m_EC_0, pars.m_EN_0, M_V_0_tWw_VE2023_ML];
% t = [1 (tWw_VE2023_ML(2:end,1) * 24)']; %in hours
% 
% forcingVariables.CO2 = CO2.tWw_VE2023_ML; 
% forcingVariables.lightIntensity = light.tWw_VE2023_ML; 
% forcingVariables.N = nitrogen.tWw_VE2023_ML; 
% forcingVariables.P =  phosphorous.tWw_VE2023_ML; 
% forcingVariables.photoPeriod =  photoPeriod.tWw_VE2023_ML; 
% 
% 
% [~, mECENV_tWw_VE2023_ML] = ode89(@(t,mECENV_tWw_VE2023_ML) get_mECENV(t,  mECENV_tWw_VE2023_ML, pars, TC_tWw_VE2023_ML, forcingVariables), t, CI_tWw_VE2023_ML);
% % mECENV_Ww_28(1,:) = [] ;
% 
% 
% m_EC = mECENV_tWw_VE2023_ML(:,1); % mol C mol V-1
% m_EN = mECENV_tWw_VE2023_ML(:,2); % mol N mol V-1
% M_V = mECENV_tWw_VE2023_ML(:,3);  % mol V
% 
% Wd_tWw_VE2023_ML = (pars.w_V + (m_EN .* pars.w_EN) + (m_EC .* pars.w_EC)) .* M_V; %g dW, dry weight
% Ww_tWw_VE2023_ML = Wd_tWw_VE2023_ML ./ dwratio_SF   ; %g wW, wet weight
% 
% 
% % N:C ratio during simulation and total percentages of C and N
% M_EC = m_EC(end) .* M_V(end); % mol C
% M_EN = m_EN(end) .* M_V(end); % mol N
% 
% C_total = (pars.n_CV * M_V(end)) + (pars.n_CEC * M_EC); % mol C  
% N_total = (pars.n_NV * M_V(end)) + (pars.n_NEN * M_EN); % mol N
% % P_total  = (pars.n_PV * M_V(end)); % mol P
% 
% percC_gdW = (C_total * pars.w_C ./ Wd_tWw_VE2023_ML(end) ) * 100; % %C (gdW)
% percN_gdW = (N_total  * pars.w_N ./ Wd_tWw_VE2023_ML(end) ) * 100; % %N (gdW)
% % percP_gdW  = (P_total * pars.w_P./ Wd_tWw_VE2023_ML(end) ) * 100; % %P (gdW)
% 
% 
% 
% prdData.tWw_VE2023_ML = Ww_tWw_VE2023_ML;
% prdData.C_gdW_VE2023_ML = percC_gdW;
% prdData.N_gdW_VE2023_ML = percN_gdW;
% 
% % HL
% Wd_0_tWw_VE2023_HL = Ww0.tWw_VE2023_HL * dwratio_SF;  %g dW, dry weight
% M_V_0_tWw_VE2023_HL = Wd_0_tWw_VE2023_HL  / (pars.w_V+   (pars.m_EN_0* pars.w_EN) +   (pars.m_EC_0* pars.w_EC)); % mol V, structural initial mass 
% CI_tWw_VE2023_HL = [pars.m_EC_0, pars.m_EN_0, M_V_0_tWw_VE2023_HL];
% t = [1 (tWw_VE2023_HL(2:end,1) * 24)']; %in hours
% 
% forcingVariables.CO2 = CO2.tWw_VE2023_HL; 
% forcingVariables.lightIntensity = light.tWw_VE2023_HL; 
% forcingVariables.N = nitrogen.tWw_VE2023_HL; 
% forcingVariables.P =  phosphorous.tWw_VE2023_HL; 
% forcingVariables.photoPeriod =  photoPeriod.tWw_VE2023_HL; 
% 
% 
% [~, mECENV_tWw_VE2023_HL] = ode89(@(t,mECENV_tWw_VE2023_HL) get_mECENV(t,  mECENV_tWw_VE2023_HL, pars, TC_tWw_VE2023_HL, forcingVariables), t, CI_tWw_VE2023_HL);
% % mECENV_Ww_28(1,:) = [] ;
% 
% 
% m_EC = mECENV_tWw_VE2023_HL(:,1); % mol C mol V-1
% m_EN = mECENV_tWw_VE2023_HL(:,2); % mol N mol V-1
% M_V = mECENV_tWw_VE2023_HL(:,3);  % mol V
% 
% Wd_tWw_VE2023_HL = (pars.w_V + (m_EN .* pars.w_EN) + (m_EC .* pars.w_EC)) .* M_V; %g dW, dry weight
% Ww_tWw_VE2023_HL = Wd_tWw_VE2023_HL ./ dwratio_SF   ; %g wW, wet weight
% 
% 
% % N:C ratio during simulation and total percentages of C and N
% M_EC = m_EC(end) .* M_V(end); % mol C
% M_EN = m_EN(end) .* M_V(end); % mol N
% 
% C_total = (pars.n_CV * M_V(end)) + (pars.n_CEC * M_EC); % mol C  
% N_total = (pars.n_NV * M_V(end)) + (pars.n_NEN * M_EN); % mol N
% % P_total  = (pars.n_PV * M_V(end)); % mol P
% 
% percC_gdW = (C_total * pars.w_C ./ Wd_tWw_VE2023_HL(end) ) * 100; % %C (gdW)
% percN_gdW = (N_total  * pars.w_N ./ Wd_tWw_VE2023_HL(end) ) * 100; % %N (gdW)
% % percP_gdW  = (P_total * pars.w_P./ Wd_tWw_VE2023_HL(end) ) * 100; % %P (gdW)
% 
% 
% 
% prdData.tWw_VE2023_HL = Ww_tWw_VE2023_HL;
% prdData.C_gdW_VE2023_HL = percC_gdW;
% prdData.N_gdW_VE2023_HL = percN_gdW;

%% Created datasets
 options = odeset('RelTol',1e-10,'AbsTol',1e-12);
if contains(pwd,'Identifiability')
    %% TC_WW 
    datasets=[]; 

    for i = 1:11
        datasetName = sprintf("TC_Ww_%d", i);
        datasets = [datasets datasetName]; % List of datasets
    end
    
    
    
    for i = 1:length(datasets)
        dataset = char(datasets(i));
        TC = tempcorr(temp.(dataset), pars.T_ref, pars_T);  %-
        
        Wd_0 = Ww0.(dataset) * dwratio_SF;  
        M_V_0 = Wd_0 / (pars.w_V + (pars.m_EN_0 * pars.w_EN) + (pars.m_EC_0 * pars.w_EC));  
        CI = [pars.m_EC_0, pars.m_EN_0, M_V_0];  
       
        t = [0 (eval([dataset,'(2:end,1)'])*24)'];  % Convert time to hours

        forcingVariables.CO2 = CO2.(dataset);
        forcingVariables.lightIntensity = light.(dataset);
        forcingVariables.N = nitrogen.(dataset);
        forcingVariables.P = phosphorous.(dataset);
        forcingVariables.photoPeriod = photoPeriod.(dataset);
      
        [~, mECENV] = ode89(@(t, mECENV) get_mECENV(t, mECENV, pars, TC, forcingVariables), t, CI, options);

      
        m_EC = mECENV(:,1);
        m_EN = mECENV(:,2);
        M_V = mECENV(:,3);
    
        Wd = (pars.w_V + (m_EN * pars.w_EN) + (m_EC * pars.w_EC)) .* M_V;%g dW, dry weight
        Ww = Wd ./ dwratio_SF;
     

       
        prdData.(dataset) = Ww;
    end

    %% TC_WW_twoDays
    datasets=[]; 

    for i = 1:11
        datasetName = sprintf("TC_Ww_twoDays_%d", i);
        datasets = [datasets datasetName]; % List of datasets
    end
    
    
    
    for i = 1:length(datasets)
        dataset = char(datasets(i));
        TC = tempcorr(temp.(dataset), pars.T_ref, pars_T);  %-
        Wd_0 = Ww0.(dataset) * dwratio_SF;  
        M_V_0 = Wd_0 / (pars.w_V + (pars.m_EN_0 * pars.w_EN) + (pars.m_EC_0 * pars.w_EC));  
        CI = [pars.m_EC_0, pars.m_EN_0, M_V_0];  

        t = [0 (eval([dataset,'(2:end,1)'])*24)'];  % Convert time to hours
        nt = length(t);
        
        forcingVariables.CO2 = CO2.(dataset);
        forcingVariables.lightIntensity = light.(dataset);
        forcingVariables.N = nitrogen.(dataset);
        forcingVariables.P = phosphorous.(dataset);
        forcingVariables.photoPeriod = photoPeriod.(dataset);
        
        [~, mECENV_twodays] = ode89(@(t, mECENV) get_mECENV(t, mECENV, pars, TC, forcingVariables), t, CI, options);        
        
        if nt == 2
            mECENV_twodays = mECENV_twodays([1 end],:);
        end

    
    
        m_EC = mECENV_twodays(:,1);
        m_EN = mECENV_twodays(:,2);
        M_V = mECENV_twodays(:,3);
    
        Wd = (pars.w_V + (m_EN * pars.w_EN) + (m_EC * pars.w_EC)) .* M_V;
        Ww = Wd ./ dwratio_SF;
    
        prdData.(dataset) = Ww;
    end

    %% Starvation
    datasets=[]; 

    for i = 1
        datasetName = sprintf("Starvation_%d", i);
        datasets = [datasets datasetName]; % List of datasets
    end
    
    
    
    for i = 1:length(datasets)
        dataset = char(datasets(i));
        TC = tempcorr(temp.(dataset), pars.T_ref, pars_T);  %-
  
        Wd_0 = Ww0.(sprintf('tWw_%s', dataset)) * dwratio_SF;  
        M_V_0 = Wd_0 / (pars.w_V + (pars.m_EN_0 * pars.w_EN) + (pars.m_EC_0 * pars.w_EC));  
        CI = [pars.m_EC_0, pars.m_EN_0, M_V_0];  
        t = [0 (eval([(sprintf('tWw_%s', dataset)),'(2:end,1)'])*24)'];  % Convert time to hours
    
        forcingVariables.CO2 = CO2.(dataset);
        forcingVariables.lightIntensity = light.(dataset);
        forcingVariables.N = nitrogen.(dataset);
        forcingVariables.P = phosphorous.(dataset);
        forcingVariables.photoPeriod = photoPeriod.(dataset);
      
        [~, mECENV] = ode89(@(t, mECENV) get_mECENV(t, mECENV, pars, TC, forcingVariables), t, CI, options);

        if infosgr2 == 0
          info = 0; 
          prdData = []; 
        end

    
        m_EC = mECENV(:,1);
        m_EN = mECENV(:,2);
        M_V = mECENV(:,3);
    
        Wd = (pars.w_V + (m_EN * pars.w_EN) + (m_EC * pars.w_EC)) .* M_V;
        Ww = Wd ./ dwratio_SF;

        %% N:C ratio during simulation and total percentages of C and N
        M_EC = m_EC .* M_V; % mol C
        M_EN = m_EN .* M_V; % mol N
        
        C_total = (pars.n_CV * M_V) + (pars.n_CEC * M_EC); % mol C  
        N_total = (pars.n_NV * M_V) + (pars.n_NEN * M_EN); % mol N
        P_total  = (pars.n_PV * M_V); % mol P

        CNtotalRatio = C_total ./ N_total ;  % mol C mol N-1
        CPtotalRatio = C_total ./ P_total ;  % mol C mol N-1
        NPtotalRatio = N_total ./ P_total ;  % mol C mol N-1



  
        prdData.(sprintf('tWw_%s', dataset)) = Ww;
        prdData.(sprintf('CNratio_%s', dataset)) = CNtotalRatio; 

     
    end
    

    
    %% Nutrient uptake
    
    % N
    datasets = []; 
    for i = 1:3
        datasetName = sprintf("N_uptake_%d", i);
        datasets = [datasets datasetName]; % List of datasets
    end
    
    for i = 1:length(datasets)
        dataset = char(datasets(i));
        TC_N = tempcorr(temp.(dataset), pars.T_ref, pars_T);  %-
        Ww0_N = Ww0.(dataset); %g Ww, initial wet biomass 
        Wd0_N = Ww0_N * pars.dwratio;  %g dW, initial dry weight biomass
        M_V_0_N = Wd0_N/ (pars.w_V +   (pars. m_EN_0 * pars.w_EN) +   (pars.m_EC_0 * pars.w_EC)); % mol V, structural intial mass
        j_ENAm_CT = pars.j_ENAm * TC_N; 
        j_EN_A =  j_ENAm_CT * (eval([dataset '(:,1)']) ./ (eval([dataset '(:,1)']) + pars.K_N)); % mol N mol V-1 h-1
        j_EN_A_dryweight =  j_EN_A * M_V_0_N / Wd0_N;  %mol N g dW-1 h-1
        prdData.(dataset) = j_EN_A_dryweight; 
    
    end 
end 



info = 1;  % Laure check info = 1 
%(Maria, sometimes there's conditions for parameters or results inside the predict, so it doesn't continue if the condition isn't filled)
% in that case,info == 0 to stop estimation 
end


