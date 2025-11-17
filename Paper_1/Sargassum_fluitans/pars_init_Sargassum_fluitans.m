function [par, metaPar, txtPar] = pars_init_Sargassum_fluitans(metaData)

%% Parameters at T_ref = 25°C
metaPar.model = 'algal';
par.j_ENAm     = 1.013610402210901e-04;    free.j_ENAm     = 1;   units.j_ENAm     = 'mol M_E_N mol M_V^{-1} h^{-1}';     label.j_ENAm     = 'Maximum volume specific nitrogen assimilation';
par.K_N        = 3.816063417553861e-07;                    free.K_N        = 1;   units.K_N        = 'mol N L^{-1}';                     label.K_N        = 'Half-saturation concentration for NO3- and NH4+ uptake';
par.j_Cm       = 0.01076;                  free.j_Cm       = 0;   units.j_Cm       = 'mol CO2 mol M_V^{-1} h^{-1}';       label.j_Cm       = 'Maximum volume specific CO2 uptake rate';
par.K_C        = 4e-07;                    free.K_C        = 0;   units.K_C        = 'mol CO2 L^{-1}';                    label.K_C        = 'Half-saturation concentration for CO2 uptake';
par.rho_PSU    = 0.5;                      free.rho_PSU    = 0;   units.rho_PSU    = 'mol PSU mol M_V^{-1}';              label.rho_PSU    = 'Photosynthetic unit (PSU) density';
par.eps_I      = 1.390630083197465;                      free.eps_I      = 1;   units.eps_I      = 'm^2 mol PSU^{-1}';                  label.eps_I      = 'Effective photon binding efficiency, compound parameter';
par.k_I        = 0.068313024589564;                   free.k_I        = 1;   units.k_I        = 'mol photon mol PSU^{-1} h^{-1}';    label.k_I        = 'Dissociation rate of photosynthetic products';
par.y_IEC      = 10;                       free.y_IEC      = 0;   units.y_IEC      = 'mol photon mol M_E_C^{-1}';         label.y_IEC      = 'Yield of C reserve per photon';
par.y_CEC      = 1;                        free.y_CEC      = 0;   units.y_CEC      = 'mol CO2 mol M_E_C^{-1}';            label.y_CEC      = 'Yield of C reserve per CO2';
par.y_OI       = 0.125;                    free.y_OI       = 0;   units.y_OI       = 'mol O2 mol photon^{-1}';             label.y_OI       = 'Yield factor of photon to O2';
par.j_ECAm     = 0.4047;                   free.j_ECAm     = 0;   units.j_ECAm     = 'mol M_E_C mol M_V^{-1} h^{-1}';     label.j_ECAm     = 'Maximum volume specific carbon assimilation';
par.k_EN       = 0.078921429491994;                   free.k_EN       = 1;   units.k_EN       = 'h^{-1}';                            label.k_EN       = 'N reserve turnover rate';
par.k_EC       = 0.029650354797646;                   free.k_EC       = 1;   units.k_EC       = 'h^{-1}';                            label.k_EC       = 'C reserve turnover rate';
par.j_ENM      = 3.707313010922020e-06;               free.j_ENM      = 1;   units.j_ENM      = 'mol M_E_N mol M_V^{-1} h^{-1}';     label.j_ENM      = 'Volume specific maintenance cost paid by N reserve';
par.j_ECM      = 3.844120163223560e-07;               free.j_ECM      = 1;   units.j_ECM      = 'mol M_E_C mol M_V^{-1} h^{-1}';     label.j_ECM      = 'Volume specific maintenance cost paid by C reserve';
par.y_ENV      = 0.04;                     free.y_ENV      = 0;   units.y_ENV      = 'mol M_E_N mol M_V^{-1}';            label.y_ENV      = 'Yield factor of N reserve to structure';
par.y_ECV      = 1.25;                     free.y_ECV      = 0;   units.y_ECV      = 'mol M_E_C mol M_V^{-1}';            label.y_ECV      = 'Yield factor of C reserve to structure';
par.kap_EN     = 0.9;                      free.kap_EN     = 0;   units.kap_EN     = '-';                                 label.kap_EN     = 'Fraction of rejection flux incorporated back in N reserve';
par.kap_EC     = 0.9;                      free.kap_EC     = 0;   units.kap_EC     = '-';                                 label.kap_EC     = 'Fraction of rejection flux incorporated back in C reserve';
par.T_A        = 6603;                     free.T_A        = 0;   units.T_A        = 'K';                                 label.T_A        = 'Arrhenius temperature';
par.T_ref      = 298.15;                   free.T_ref      = 0;   units.T_ref      = 'K';                                 label.T_ref      = 'Reference temperature';
par.T_H        = 304.3;                    free.T_H        = 0;   units.T_H        = 'K';                                 label.T_H        = 'Upper boundary of temperature tolerance';
par.T_L        = 288.9;                    free.T_L        = 0;   units.T_L        = 'K';                                 label.T_L        = 'Lower boundary of temperature tolerance';
par.T_AH       = 5.801e+04;                free.T_AH       = 0;   units.T_AH       = 'K';                                 label.T_AH       = 'Arrhenius temperature beyond T_H';
par.T_AL       = 5.462e+04;                free.T_AL       = 0;   units.T_AL       = 'K';                                 label.T_AL       = 'Arrhenius temperature beyond T_L';
par.phi_dw     = 0.178;                    free.phi_dw     = 0;   units.phi_dw     = '-';                                 label.phi_dw     = 'Dry to wet weight ratio (gdw/gww)';
par.psi_s      = 0.1727;                   free.psi_s      = 0;   units.psi_s      = 'g ww cm^{-2}';                      label.psi_s      = 'Coefficient of linear regression from Chambon';
par.zalgal     = 1;                        free.zalgal     = 0;   units.zalgal     = '-';                                 label.zalgal     = 'Zoom factor for algal species';
par.x_ash      = 0.4;                      free.x_ash      = 0;   units.x_ash      = '-';                                 label.x_ash      = 'Organic proportion in total algae dry weight';
par.x_moist    = 0.1;                      free.x_moist    = 0;   units.x_moist    = '-';                                 label.x_moist    = 'Moisture proportion in total algae dry weight';
% Initial conditions of algae
par.m_EC_0     = 0.002;                    free.m_EC_0     = 0;   units.m_EC_0     = 'mol C mol V^{-1}';                  label.m_EC_0     = 'Initial carbon reserve density';
par.m_EN_0     = 0.01;                     free.m_EN_0     = 0;   units.m_EN_0     = 'mol N mol V^{-1}';                  label.m_EN_0     = 'Initial nitrogen reserve density';
par.m_EC_0_bis = 0.002;                    free.m_EC_0_bis = 0;   units.m_EC_0_bis = 'mol C mol V^{-1}';                  label.m_EC_0_bis = 'Initial carbon reserve density (bis)';
par.m_EN_0_bis = 0.01;                     free.m_EN_0_bis = 0;   units.m_EN_0_bis = 'mol N mol V^{-1}';                  label.m_EN_0_bis = 'Initial nitrogen reserve density (bis)';




%% set chemical parameters 
[par, units, label, free] = add_chem_Sargasssum(par,units,label,free) ; 

%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 