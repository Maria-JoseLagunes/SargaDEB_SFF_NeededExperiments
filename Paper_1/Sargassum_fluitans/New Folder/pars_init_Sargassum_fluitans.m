function [par, metaPar, txtPar] = pars_init_Sargassum_fluitans(metaData)

metaPar.model = 'algal'; 

par.j_ENAm = 2.152616415995544e-04;     free.j_ENAm = 1;   units.j_ENAm = 'mol N mol V–1 h–1';    label.j_ENAm = 'Maximum volume specific nitrogen assimilation';
%  par.j_ENAm = 0.00015 ;     free.j_ENAm = 1;   units.j_ENAm = 'mol N mol V–1 h–1';    label.j_ENAm = 'Maximum volume specific nitrogen assimilation';
par.K_N = 1e-06;        free.K_N = 1;      units.K_N = 'mol NO3– and NO2– L–1';    label.K_N = 'Half-saturation concentration for NO3– and NO2– uptake';
% par.K_N = 2.5e-06;        free.K_N = 1;      units.K_N = 'mol NO3– and NO2– L–1';    label.K_N = 'Half-saturation concentration for NO3– and NO2– uptake';
par.j_CO2m = 0.01076;      free.j_CO2m = 0;   units.j_CO2m = 'mol CO2 mol V–1 h–1';   label.j_CO2m = 'Maximum volume specific CO2 uptake rate';
% par.j_CO2m = 0.0075;      free.j_CO2m = 1;   units.j_CO2m = 'mol CO2 mol V–1 h–1';   label.j_CO2m = 'Maximum volume specific CO2 uptake rate';
par.K_C = 4e-07;          free.K_C = 0;      units.K_C = 'mol CO2 L–1';              label.K_C = 'Half-saturation concentration for CO2 uptake';

% Before parameters (just to keep record of the change for future init
% simulations)
% par.rho_PSU = 0.5;        free.rho_PSU = 0;  units.rho_PSU = 'mol PSU mol V–1';      label.rho_PSU = 'Photosynthetic unit (PSU) density';
% par.beta_I = 0.60688;     free.beta_I = 0;   units.beta_I = '-';                     label.beta_I = 'Binding probability of a photon to a free light SU';
% par.alpha_I = 0.31039;    free.alpha_I = 0;  units.alpha_I = 'm2 mol PSU–1';         label.alpha_I = 'Specific photon arrival cross section';
% par.k_I = 0.024972;       free.k_I = 0;      units.k_I = 'mol γ mol PSU–1 h–1';      label.k_I = 'Dissociation rate of photosynthetic products';
% par.y_IC = 2.8446;        free.y_IC = 0;     units.y_IC = 'mol γ mol C–1';           label.y_IC = 'Yield of C reserve per photon';

% New estimation parameters
par.rho_PSU = 0.5;        free.rho_PSU = 0;  units.rho_PSU = 'mol PSU mol V–1';      label.rho_PSU = 'Photosynthetic unit (PSU) density';

% par.beta_I = 1;           free.beta_I = 0;   units.beta_I = '-';                     label.beta_I = 'Binding probability of a photon to a free light SU';
% par.alpha_I = 0.109 ;     free.alpha_I = 0;  units.alpha_I = 'm2 mol PSU–1';         label.alpha_I = 'Specific photon arrival cross section';

% par.alpha_I = 0.5;     free.alpha_I = 0;  units.alpha_I = 'm2 mol PSU–1';         label.alpha_I = 'Specific photon arrival cross section';

par.eps_I = 0.5 ;     free.eps_I = 1;  units.eps_I = 'm2 mol PSU–1';         label.eps_I = 'Effective photon binding efficiency, compound parameter';

par.k_I = 0.1076 * 2;        free.k_I = 0;      units.k_I = 'mol γ mol PSU–1 h–1';      label.k_I = 'Dissociation rate of photosynthetic products';
% par.k_I = 0.075;        free.k_I = 1;      units.k_I = 'mol γ mol PSU–1 h–1';      label.k_I = 'Dissociation rate of photosynthetic products';


par.y_IC = 10;            free.y_IC = 0;     units.y_IC = 'mol γ mol C–1';           label.y_IC = 'Yield of C reserve per photon';
par.y_CO2C = 1;           free.y_CO2C = 0;   units.y_CO2C = 'mol CO2 mol C–1';       label.y_CO2C = 'Yield of C reserve per CO2';
par.y_LO2 = 0.125;        free.y_LO2 = 0;    units.y_LO2 = 'mol O2 mol γ–1';         label.y_LO2 = 'Yield factor of photon to O2';

par.j_ECAm = 0.4047;       free.j_ECAm = 0;   units.j_ECAm = 'mol C mol V–1 h–1';     label.j_ECAm = 'Maximum volume specific carbon assimilation';
% par.j_ECAm = 0.282;       free.j_ECAm = 1;   units.j_ECAm = 'mol C mol V–1 h–1';     label.j_ECAm = 'Maximum volume specific carbon assimilation';

% par.k_EN = 0.04;          free.k_EN = 1;     units.k_EN = 'h–1';                     label.k_EN = 'N reserve turnover rate';
% par.k_EC = 0.02;          free.k_EC = 1;     units.k_EC = 'h–1';                     label.k_EC = 'C reserve turnover rate';

par.k_EN = 0.0574 ;          free.k_EN = 1;     units.k_EN = 'h–1';                     label.k_EN = 'N reserve turnover rate';
par.k_EC = 0.0287 ;          free.k_EC = 1;     units.k_EC = 'h–1';                     label.k_EC = 'C reserve turnover rate';

% par.j_ENM = 4e-06;        free.j_ENM = 1;    units.j_ENM = 'mol N mol V–1 h–1';      label.j_ENM = 'Volume specific maintenance cost paid by N reserve';
par.j_ENM = 5.7403e-06;        free.j_ENM = 1;    units.j_ENM = 'mol N mol V–1 h–1';      label.j_ENM = 'Volume specific maintenance cost paid by N reserve';

% par.j_ECM = 1e-06;        free.j_ECM = 1;    units.j_ECM = 'mol C mol V–1 h–1';      label.j_ECM = 'Volume specific maintenance cost paid by C reserve';
par.j_ECM = 1.4351e-06;        free.j_ECM = 1;    units.j_ECM = 'mol C mol V–1 h–1';      label.j_ECM = 'Volume specific maintenance cost paid by C reserve';

par.y_ENV = 0.04;           free.y_ENV = 0;    units.y_ENV = 'mol N mol V–1';          label.y_ENV = 'Yield factor of N reserve to structure';
par.y_ECV = 1.25;            free.y_ECV = 0;    units.y_ECV = 'mol C mol V-1';          label.y_ECV = 'Yield factor of C reserve to structure';

par.kap_EN = 0.9;         free.kap_EN = 0;   units.kap_EN = '-';                     label.kap_EN = 'Fraction of rejection flux incorporated back in i-reserve';
par.kap_EC = 0.9;         free.kap_EC = 0;   units.kap_EC = '-';                     label.kap_EC = 'Fraction of rejection flux incorporated back in i-reserve';


% Without modification of Arrhenius vector at 28°C
% par.T_A = 3821.1716;      free.T_A = 0;      units.T_A = 'K';                        label.T_A = 'Arrhenius temperature';
% par.T_ref = 298.15;       free.T_ref = 0;    units.T_ref = 'K';                      label.T_ref = 'Reference temperature';
% par.T_H = 304.9042;       free.T_H = 0;      units.T_H = 'K';                        label.T_H = 'Upper boundary of temp;rature tolerance';
% par.T_L = 290.1967;       free.T_L = 0;      units.T_L = 'K';                        label.T_L = 'Lower boundary of temperature tolerance';
% par.T_AH = 57585.4945;    free.T_AH = 0;     units.T_AH = 'K';                       label.T_AH = 'Arrhenius temperature outside of T_H';
% par.T_AL = 39692.2025;    free.T_AL = 0;     units.T_AL = 'K';                       label.T_AL = 'Arrhenius temperature outside of T_L';

% Venolia with Tref TH and TL modified by hand 
% par.T_A = 6314.3;      free.T_A = 1;      units.T_A = 'K';                        label.T_A = 'Arrhenius temperature';
% par.T_ref = 298.15;       free.T_ref = 0;    units.T_ref = 'K';                      label.T_ref = 'Reference temperature';
% par.T_H = C2K(31);       free.T_H = 0;      units.T_H = 'K';                        label.T_H = 'Upper boundary of temp;rature tolerance';
% par.T_L = C2K(18);       free.T_L = 0;      units.T_L = 'K';                        label.T_L = 'Lower boundary of temperature tolerance';
% par.T_AH = 18702;    free.T_AH = 1;     units.T_AH = 'K';                       label.T_AH = 'Arrhenius temperature outside of T_H';
% par.T_AL = 4391.9;    free.T_AL = 1;     units.T_AL = 'K';                       label.T_AL = 'Arrhenius temperature outside of T_L';

par.T_A = 6603;          free.T_A = 0;      units.T_A = 'K';                        label.T_A = 'Arrhenius temperature';
par.T_ref = 298.15;      free.T_ref = 0;    units.T_ref = 'K';                      label.T_ref = 'Reference temperature';
par.T_H = 304.3;         free.T_H = 0;      units.T_H = 'K';                        label.T_H = 'Upper boundary of temp;rature tolerance';
par.T_L = 288.9;         free.T_L = 0;      units.T_L = 'K';                        label.T_L = 'Lower boundary of temperature tolerance';
par.T_AH = 5.801e+04;    free.T_AH = 0;     units.T_AH = 'K';                       label.T_AH = 'Arrhenius temperature outside of T_H';
par.T_AL = 5.462e+04;    free.T_AL = 0;     units.T_AL = 'K';                       label.T_AL = 'Arrhenius temperature outside of T_L';


par.dwratio = 0.178;      free.dwratio = 0;   units.dwratio = '-';        label.dwratio = 'Dry to wet (gdw/gww) weight ratio';

par.smcoeff = 0.1727;     free.smcoeff = 0;   units.smcoeff = 'g ww cm-2';        label.smcoeff = 'coefficient of linear regression from Chambon'; 
par.zalgal = 1;           free.zalgal = 0;     units.zalgal = '-';        label.zalgal = 'Zoom factor for algal species'; 

par.x_orga = 0.5;         free.x_orga = 0;   units.x_orga = '-';   label.x_orga = 'Organic proportion in total algae dry weight'; 

% %Initial conditions of algae
% par.m_EC_0 = 0.002  ;   free.m_EC_0 = 1; units.m_EC_0 = 'mol C / mol V'; label.m_EC_0 = 'initial carbon reserve density'; 
% par.m_EN_0 = 0.01 ; free.m_EN_0 = 1; units.m_EN_0 = 'mol N / mol V'; label.m_EN_0 = 'initial nitrogen reserve density'; 

par.m_EC_0 = 0.002  ;   free.m_EC_0 = 0; units.m_EC_0 = 'mol C / mol V'; label.m_EC_0 = 'initial carbon reserve density'; 
par.m_EN_0 = 0.01 ; free.m_EN_0 = 0; units.m_EN_0 = 'mol N / mol V'; label.m_EN_0 = 'initial nitrogen reserve density'; 

par.m_EC_0_bis = 0.002;   free.m_EC_0_bis = 0; units.m_EC_0_bis = 'mol C / mol V'; label.m_EC_0_bis = 'initial carbon reserve density'; 
par.m_EN_0_bis = 0.01; free.m_EN_0_bis = 0; units.m_EN_0_bis = 'mol N / mol V'; label.m_EN_0_bis = 'initial nitrogen reserve density'; 



%% set chemical parameters 
[par, units, label, free] = add_chem_Sargasssum(par,units,label,free) ; 

%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
end
