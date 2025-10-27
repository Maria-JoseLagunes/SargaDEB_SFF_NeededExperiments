function [par, metaPar, txtPar] = pars_init_Sargassum_fluitans_Arrhenius(metaData)

metaPar.model = 'algal'    ; 

par.j_ENAm = 0.00015;     free.j_ENAm = 0;   units.j_ENAm = 'mol N mol V–1 h–1';    label.j_ENAm = 'Maximum volume specific nitrogen assimilation';
par.K_N = 2.5e-06;        free.K_N = 0;      units.K_N = 'mol NO3– and NO2– L–1';    label.K_N = 'Half-saturation concentration for NO3– and NO2– uptake';
par.j_CO2m = 0.0075;      free.j_CO2m = 0;   units.j_CO2m = 'mol CO2 mol V–1 h–1';   label.j_CO2m = 'Maximum volume specific CO2 uptake rate';
par.K_C = 4e-07;          free.K_C = 0;      units.K_C = 'mol CO2 L–1';              label.K_C = 'Half-saturation concentration for CO2 uptake';
par.rho_PSU = 0.5;        free.rho_PSU = 0;  units.rho_PSU = 'mol PSU mol V–1';      label.rho_PSU = 'Photosynthetic unit (PSU) density';
% par.beta_I = 0.5;     free.beta_I = 0;   units.beta_I = '-';                     label.beta_I = 'Binding probability of a photon to a free light SU';
% par.alpha_I = 1;    free.alpha_I = 0;  units.alpha_I = 'm2 mol PSU–1';         label.alpha_I = 'Specific photon arrival cross section';
par.eps_I = 0.5 ;     free.eps_I = 0;  units.eps_I = 'm2 mol PSU–1';         label.eps_I = 'Effective photon binding efficiency, compound parameter';
par.k_I = 0.075;       free.k_I = 0;      units.k_I = 'mol γ mol PSU–1 h–1';      label.k_I = 'Dissociation rate of photosynthetic products';
par.y_IC = 10;        free.y_IC = 0;     units.y_IC = 'mol γ mol C–1';           label.y_IC = 'Yield of C reserve per photon';
par.y_CO2C = 1;           free.y_CO2C = 0;   units.y_CO2C = 'mol CO2 mol C–1';       label.y_CO2C = 'Yield of C reserve per CO2';
par.y_LO2 = 0.125;        free.y_LO2 = 0;    units.y_LO2 = 'mol O2 mol γ–1';         label.y_LO2 = 'Yield factor of photon to O2';
par.j_ECAm = 0.282;       free.j_ECAm = 0;   units.j_ECAm = 'mol C mol V–1 h–1';     label.j_ECAm = 'Maximum volume specific carbon assimilation';
par.k_EN = 0.04;          free.k_EN = 0;     units.k_EN = 'h–1';                     label.k_EN = 'N reserve turnover rate';
par.k_EC = 0.02;          free.k_EC = 0;     units.k_EC = 'h–1';                     label.k_EC = 'C reserve turnover rate';
par.j_ENM = 4e-06;        free.j_ENM = 0;    units.j_ENM = 'mol N mol V–1 h–1';      label.j_ENM = 'Volume specific maintenance cost paid by N reserve';
par.j_ECM = 1e-06;        free.j_ECM = 0;    units.j_ECM = 'mol C mol V–1 h–1';      label.j_ECM = 'Volume specific maintenance cost paid by C reserve';
par.y_ENV = 0.04;         free.y_ENV = 0;    units.y_ENV = 'mol N mol V–1';          label.y_ENV = 'Yield factor of N reserve to structure';
par.y_ECV = 1;            free.y_ECV = 0;    units.y_ECV = 'mol C mol V-1';          label.y_ECV = 'Yield factor of C reserve to structure';
par.kap_EN = 0.9;         free.kap_EN = 0;   units.kap_EN = '-';                     label.kap_EN = 'Fraction of rejection flux incorporated back in i-reserve';
par.kap_EC = 0.9;         free.kap_EC = 0;   units.kap_EC = '-';                     label.kap_EC = 'Fraction of rejection flux incorporated back in i-reserve';

load("simulationParameters.mat")

par.T_A = simulationParameters.T_A;           free.T_A = 1;      units.T_A = 'K';                        label.T_A = 'Arrhenius temperature';
par.T_ref = C2K(simulationParameters.T_ref);       free.T_ref = 0;    units.T_ref = 'K';                      label.T_ref = 'Reference temperature';
par.T_H = C2K(simulationParameters.T_H);         free.T_H = 1;      units.T_H = 'K';                        label.T_H = 'Upper boundary of temp;rature tolerance';
par.T_L = C2K(simulationParameters.T_L);         free.T_L = 1;      units.T_L = 'K';                        label.T_L = 'Lower boundary of temperature tolerance';
par.T_AH = simulationParameters.T_AH;         free.T_AH = 1;     units.T_AH = 'K';                       label.T_AH = 'Arrhenius temperature outside of T_H';
par.T_AL = simulationParameters.T_AL;         free.T_AL = 1;     units.T_AL = 'K';                       label.T_AL = 'Arrhenius temperature outside of T_L';



par.dwratio = 0.178;      free.dwratio = 0;   units.dwratio = '-';        label.dwratio = 'Dry to wet (gdw/gww) weight ratio';
 
par.smcoeff = 0.1727;     free.smcoeff = 0;   units.smcoeff = 'g ww cm-2';        label.smcoeff = 'coefficient of linear regression from Chambon'; 
    
par.zalgal = 1;           free.zalgal = 0;     units.zalgal = '-';        label.zalgal = 'Zoom factor for algal species'; 

par.x_orga = 0.5;         free.x_orga = 0;   units.x_orga = '-';   label.x_orga = 'Organic proportion in total algae dry weight'; 

% %Initial conditions of algae
par.m_EC_0 = 0.002;   free.m_EC_0 = 0; units.m_EC_0 = 'mol C / mol V'; label.m_EC_0 = 'initial carbon reserve density'; 
par.m_EN_0 = 0.01; free.m_EN_0 = 0; units.m_EN_0 = 'mol N / mol V'; label.m_EN_0 = 'initial nitrogen reserve density'; 



%% set chemical parameters 
[par, units, label, free] = add_chem_Sargasssum_Arrhenius(par,units,label,free) ; 

%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
end
