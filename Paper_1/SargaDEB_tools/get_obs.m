%% Main file to do simulations
%
% Maria José Lagunes
% First created : 2024/03/14
% ---------------------------------------------------------------------


function obs = get_obs(mECENV,J_struct, simulationConditions)


vars_pull(simulationConditions)
vars_pull(par);

m_EC = mECENV(:,1); % mol EC mol V-1
m_EN = mECENV(:,2); % mol EN mol V-1
M_V = mECENV(:,3);  % mol V

% obs.Wd = (w_V + (m_EC .*w_EC) + (m_EN .* w_EN)) .* M_V; %g dW, dry weight
% obs.Ww = obs.Wd / dwratio; %g wW, wet weight

obs.Wd_free = (w_V + (m_EC .*w_EC) + (m_EN .* w_EN)) .* M_V; %g dW, dry weight
obs.Wd = obs.Wd_free / (1 - par.x_moist - par.x_ash); 
obs.Ww = obs.Wd / dwratio; %g wW, wet weight



obs.Ww_t = obs.Ww(end) ; % g dW,  final wet weight  
obs.Ww_0 = obs.Ww(1) ;   % g  dW, initial wet weight 

obs.timeDays = round(t(end)/24); %time in days

%% Growth rate measure
obs.r_D = log2(obs.Ww_t /obs.Ww_0) / obs.timeDays; % doublings per d-1, relative growth rate
obs.t_D = 1 / obs.r_D; % time to double biomass


%% N:C ratio during simulation and total percentages of C and N
M_EC = m_EC .* M_V; % mol C
M_EN = m_EN .* M_V; % mol N

obs.C_total = (n_CV * M_V) + (n_CEC * M_EC); % mol C  
obs.N_total = (n_NV * M_V) + (n_NEN * M_EN); % mol N
obs.P_total  = (n_PV * M_V); % mol P


perc_gEC = (M_EC * w_EC ./ obs.Wd_free) * 100; %  % EC gdW
perc_gEN = (M_EN * w_EN./ obs.Wd_free) * 100; %   % EN gdW
perc_gV = (M_V * w_V./ obs.Wd_free) * 100;  %     % MV gdW

% obs.C_total = 1 + m_EC; % mol C mol V-1
% obs.N_total = n_NV + m_EN; % mol N mol V-1

obs.percC_gdW = (obs.C_total * w_C ./ obs.Wd_free) * 100; % %C (gdW)
obs.percN_gdW = (obs.N_total  * w_N ./ obs.Wd_free) * 100; % %N (gdW)
obs.percP_gdW  = (obs.P_total * w_P./ obs.Wd_free) * 100; % %P (gdW)




obs.CNtotalRatio = obs.C_total ./ obs.N_total ;  % mol C mol N-1
obs.CPtotalRatio = obs.C_total ./ obs.P_total ;  % mol C mol N-1
obs.NPtotalRatio = obs.N_total ./ obs.P_total ;  % mol C mol N-1

% obs.NCtotalRatio = obs.N_totalRatio ./ obs.C_totalRatio;  % mol N mol C-1


%% N:C ratio at the end of simulation
obs.C_endRatio = 1 + m_EC(end); % mol C mol V-1
obs.N_endRatio= n_NV + m_EN(end); % mol N mol V-1

obs.CNendRatio = obs.C_endRatio ./ obs.N_endRatio; %mol C mol N -1
% obs.NCendRatio = obs.N_endRatio ./ obs.C_endRatio;  % mol N mol C-1
%% Oxygen production rate
obs.J_O2 =  J_struct.j_I .* M_V   * y_LO2; % absolut O2 production rate, mol 02 h-1
obs.j_O2 =  J_struct.j_I .* (M_V / obs.Wd_free) * y_LO2 ; % mol 02 gdW-1  h-1

obs.j_O2 = obs.j_O2 / 1e-6; % micro mol O2 gdW-1 h-1


%% % Composition in V
obs.Wd_V = M_V * w_V; % g dW_V

 
obs.M_CV = (n_CV * M_V); % mol C in MV
obs.M_NV = (n_NV * M_V); % mol N in MV
obs.M_PV = (n_PV * M_V); % mol P in MV


obs.M_CVVV = (n_CV * M_V) ./ M_V; % mol C in MV / mol M_V --> n_CV
obs.M_NVVV = (n_NV * M_V) ./ M_V; % mol N in MV / mol M_V --> n_CEC
obs.M_PVVV = (n_PV * M_V) ./ M_V; % mol P in MV / mol M_V --> n_CEN


obs.Wd_CV = (obs.M_CV * w_C); % gdW CV
obs.Wd_NV = (obs.M_NV  * w_N); %  gdW NV
obs.Wd_PV = (obs.M_PV * w_P); %gdW PV
% 
% obs.Wd_V_with_ash = obs.Wd_V * 1 / x_orga; 

obs.percCV_gdW = (obs.Wd_CV  ./ obs.Wd_V) * 100;
obs.percNV_gdW = (obs.Wd_NV  ./ obs.Wd_V) * 100;
obs.percPV_gdW = (obs.Wd_PV  ./ obs.Wd_V) * 100;

% 
% obs.percCV_gdW_with_ash = (obs.Wd_CV  ./ obs.Wd_V_with_ash) * 100;
% obs.percNV_gdW_with_ash = (obs.Wd_NV  ./ obs.Wd_V_with_ash) * 100;
% obs.percPV_gdW_with_ash = (obs.Wd_PV  ./ obs.Wd_V_with_ash) * 100;
% 

%% Ratio structure

obs.CN_V = obs.M_CV ./ obs.M_NV ; 
obs.CP_V = obs.M_CV ./ obs.M_PV ; 
obs.NP_V = obs.M_NV ./ obs.M_PV ; 


%% % Composition in EC
obs.Wd_EC = M_EC * w_EC; % E_C in g dW 

 
obs.M_CEC = (n_CEC * M_EC); % mol C in MV
obs.M_NEC = (n_NEC * M_EC); % mol N in MV
obs.M_PEC = (n_PEC * M_EC); % mol P in MV

obs.Wd_CEC = (obs.M_CEC * w_C); % gdW CV
obs.Wd_NEC= (obs.M_NEC * w_N); %  gdW NV
obs.Wd_PEC = (obs.M_PEC * w_P); %gdW PV


obs.percCEC_gdW = (obs.Wd_CEC  ./ obs.Wd_EC) * 100;
obs.percNEC_gdW = (obs.Wd_NEC  ./ obs.Wd_EC) * 100;
obs.percPEC_gdW = (obs.Wd_PEC  ./ obs.Wd_EC) * 100;


%% %Composition in EN
obs.Wd_EN = M_EN * w_EN; % E_C in g dW 

 
obs.M_CEN = (n_CEN * M_EN); % mol C in MV
obs.M_NEN = (n_NEN * M_EN); % mol N in MV
obs.M_PEN = (n_PEN * M_EN); % mol P in MV

obs.Wd_CEN = (obs.M_CEN * w_C); % gdW CV
obs.Wd_NEN = (obs.M_NEN * w_N); %  gdW NV
obs.Wd_PEN = (obs.M_PEN * w_P); %gdW PV


obs.percCEN_gdW = (obs.Wd_CEN  ./ obs.Wd_EN) * 100;
obs.percNEN_gdW = (obs.Wd_NEN  ./ obs.Wd_EN) * 100;
obs.percPEN_gdW = (obs.Wd_PEN  ./ obs.Wd_EN) * 100;

%%
g_N_dW_min = 0.04; 
g_N_dW_max  = 0.2; 
g_C_dW =  1; 

mol_N_min_Lorena = g_N_dW_min / w_N; 
mol_N_max_Lorena = g_N_dW_max  / w_N;

mol_C_Lorena = g_C_dW / w_C; 

mol_NC_min_Lorena = mol_N_min_Lorena / mol_C_Lorena;
mol_NC_max_Lorena = mol_N_max_Lorena / mol_C_Lorena;


% g_N_dW_min_range_min = 0.016; 
% g_N_dW_min_range_max = 0.029; 
g_N_dW_min_range_av = 0.027; 

% g_N_dW_max_range_min = 0.033; 
% g_N_dW_max_range_max = 0.058; 
g_N_dW_max_range_av = 0.034; 



g_P_dW_min_range_av = 0.003; 

g_P_dW_max_range_av = 0.008; 

g_C_dW =  1; 


mol_N_min_Jouanno_average = g_N_dW_min_range_av  / w_N; 
mol_N_max_Jouanno_average = g_N_dW_max_range_av  / w_N;

mol_C_Jouanno_average= g_C_dW / w_C; 

mol_NC_min_Jouanno_average = mol_N_min_Jouanno_average / mol_C_Jouanno_average;
mol_NC_max_Jouanno_average = mol_N_max_Jouanno_average / mol_C_Jouanno_average;

mol_CN_min_Jouanno_average = 1 / mol_NC_min_Jouanno_average; 
mol_CN_max_Jouanno_average = 1 / mol_NC_max_Jouanno_average; 

%
mol_P_min_Jouanno_average = g_P_dW_min_range_av  / w_P; 
mol_P_max_Jouanno_average = g_P_dW_max_range_av  / w_P;

mol_C_Jouanno_average= g_C_dW / w_C; 

mol_PC_min_Jouanno_average = mol_P_min_Jouanno_average / mol_C_Jouanno_average;
mol_PC_max_Jouanno_average = mol_P_max_Jouanno_average / mol_C_Jouanno_average;

mol_CP_min_Jouanno_average = 1 / mol_PC_min_Jouanno_average; 
mol_CP_max_Jouanno_average = 1 / mol_PC_max_Jouanno_average; 



end


