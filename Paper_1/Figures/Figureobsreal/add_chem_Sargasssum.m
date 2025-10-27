% ---------------------------------------------------------------------
%% Addchemical parameters
% Maria José Lagunes
% First created : 2024/03/11


% Parameters (values, unit and description)
% Output 
% chemicalParameters = structure with parameters associated to the studied
% species

% ---------------------------------------------------------------------
function [par, units, label, free] = add_chem_Sargassum(par, units, label, free)

par.n_CXN = 0;     free.n_CXN = 0;   units.n_CXN = '-';   label.n_CXN = 'Chemical index of carbon in nitrogen nutrient (NO3)';
par.n_HXN = 0;     free.n_HXN = 0;   units.n_HXN = '-';   label.n_HXN = 'Chemical index of hydrogen in nitrogen nutrient (NO3)';
par.n_OXN = 3;     free.n_OXN = 0;   units.n_OXN = '-';   label.n_OXN = 'Chemical index of oxygen in nitrogen nutrient (NO3)';
par.n_NXN = 1;     free.n_NXN = 0;   units.n_NXN = '-';   label.n_NXN = 'Chemical index of nitrogen in nitrogen nutrient (NO3)';

par.n_CXC = 1;     free.n_CXC = 0;   units.n_CXC = '-';   label.n_CXC = 'Chemical index of carbon in carbon nutrient (CO2+HCO3-)';
par.n_HXC = 0.5;   free.n_HXC = 0;   units.n_HXC = '-';   label.n_HXC = 'Chemical index of hydrogen in carbon nutrient (CO2+HCO3-)';
par.n_OXC = 2.5;   free.n_OXC = 0;   units.n_OXC = '-';   label.n_OXC = 'Chemical index of oxygen in carbon nutrient (CO2+HCO3-)';
par.n_NXC = 0;     free.n_NXC = 0;   units.n_NXC = '-';   label.n_NXC = 'Chemical index of nitrogen in carbon nutrient (CO2+HCO3-)';

% par.n_CV = 106/106;    free.n_CV = 0;    units.n_CV = '-';    label.n_CV = 'Chemical index of carbon in structure from Redfield ratio (106:263:110:16:1)';
% par.n_HV = 263/106;    free.n_HV = 0;    units.n_HV = '-';    label.n_HV = 'Chemical index of hydrogen in structure from Redfield ratio (106:263:110:16:1)';
% par.n_OV = 110/106;    free.n_OV = 0;    units.n_OV = '-';    label.n_OV = 'Chemical index of oxygen in structure from Redfield ratio (106:263:110:16:1)';
% par.n_NV = 16/106;     free.n_NV = 0;    units.n_NV = '-';    label.n_NV = 'Chemical index of nitrogen in structure from Redfield ratio (106:263:110:16:1)';
% par.n_PV = 1/106;     free.n_PV = 0;    units.n_PV = '-';    label.n_PV = 'Chemical index of phosphorous in structure from Redfield ratio (106:263:110:16:1)';


%Modifying elemental composition of structure
 par.n_CV = 1;    free.n_CV = 0;    units.n_CV = '-';    label.n_CV = 'Chemical index of carbon in structure from Redfield ratio (106:263:110:16:1)';
 par.n_HV = 1.33;    free.n_HV = 0;    units.n_HV = '-';    label.n_HV = 'Chemical index of hydrogen in structure from Redfield ratio (106:263:110:16:1)';
 par.n_OV = 1;    free.n_OV = 0;    units.n_OV = '-';    label.n_OV = 'Chemical index of oxygen in structure from Redfield ratio (106:263:110:16:1)';
 par.n_NV = 0.04;     free.n_NV = 0;    units.n_NV = '-';    label.n_NV = 'Chemical index of nitrogen in structure from Redfield ratio (106:263:110:16:1)';
 par.n_PV = 0.0012;     free.n_PV = 0;    units.n_PV = '-';    label.n_PV = 'Chemical index of phosphorous in structure from Redfield ratio (106:263:110:16:1)';
 

par.n_CEN = 0;     free.n_CEN = 0;   units.n_CEN = '-';   label.n_CEN = 'Chemical index of carbon in nitrogen reserve (NO3 + NH4)';
par.n_HEN = 2;     free.n_HEN = 0;   units.n_HEN = '-';   label.n_HEN = 'Chemical index of hydrogen in nitrogen reserve (NO3 + NH4)';
par.n_OEN = 1.5;   free.n_OEN = 0;   units.n_OEN = '-';   label.n_OEN = 'Chemical index of oxygen in nitrogen reserve (NO3 + NH4)';
par.n_NEN = 1;     free.n_NEN = 0;   units.n_NEN = '-';   label.n_NEN = 'Chemical index of nitrogen in nitrogen reserve (NO3 + NH4)';
par.n_PEN = 0;     free.n_PEN = 0;   units.n_PEN = '-';   label.n_PEN = 'Chemical index of phosphorous in nitrogen reserve (NO3 + NH4)';

par.n_CEC = 1;     free.n_CEC = 0;   units.n_CEC = '-';   label.n_CEC = 'Chemical index of carbon in carbon reserve (glucose C6H12O6 --> CH2O)';
par.n_HEC = 2;     free.n_HEC = 0;   units.n_HEC = '-';   label.n_HEC = 'Chemical index of hydrogen in carbon reserve (glucose C6H12O6 --> CH2O)';
par.n_OEC = 1;     free.n_OEC = 0;   units.n_OEC = '-';   label.n_OEC = 'Chemical index of oxygen in carbon reserve (glucose C6H12O6 --> CH2O)';
par.n_NEC = 0;     free.n_NEC = 0;   units.n_NEC = '-';   label.n_NEC = 'Chemical index of nitrogen in carbon reserve (glucose C6H12O6 --> CH2O)';
par.n_PEC = 0;     free.n_PEC = 0;   units.n_PEC = '-';   label.n_PEC = 'Chemical index of phosphorous in carbon reserve (glucose C6H12O6 --> CH2O)';


par.n_CEP = 0;     free.n_CEP = 0;   units.n_CEP = '-';   label.n_CEP = 'Chemical index of carbon in phosphorus reserve (PO4)';
par.n_HEP = 0;     free.n_HEP = 0;   units.n_HEP = '-';   label.n_HEP = 'Chemical index of hydrogen in phosphorus reserve (PO4)';
par.n_OEP = 4;     free.n_OEP = 0;   units.n_OEP = '-';   label.n_OEP = 'Chemical index of oxygen in phosphorus reserve (PO4)';
par.n_NEP = 0;     free.n_NEP = 0;   units.n_NEP = '-';   label.n_NEP = 'Chemical index of nitrogen in phosphorus reserve (PO4)';
par.n_PEP = 1;     free.n_PEP = 0;   units.n_PEP = '-';   label.n_PEP = 'Chemical index of phosphorous in phosphorus reserve (PO4)';


%% atomic masses and molar weights
par.w_C = 12;   free.w_C = 0;   units.w_C = 'g/mol';   label.w_C = 'Atomic mass of carbon';
par.w_H = 1;    free.w_H = 0;   units.w_H = 'g/mol';   label.w_H = 'Atomic mass of hydrogen';
par.w_O = 16;   free.w_O = 0;   units.w_O = 'g/mol';   label.w_O = 'Atomic mass of oxygen';
par.w_N = 14;   free.w_N = 0;   units.w_N = 'g/mol';   label.w_N = 'Atomic mass of nitrogen';
par.w_P = 31;   free.w_P = 0;   units.w_P = 'g/mol';   label.w_P = 'Atomic mass of phosphorous';

par.w_V = (par.n_CV * par.w_C) + (par.n_HV * par.w_H) + (par.n_OV * par.w_O) ...
    + (par.n_NV * par.w_N) + (par.n_PV * par.w_P) ; % Can I include P in structure when it isn't in any of the reserves ? 
free.w_V = 0;    units.w_V = 'g C-mol V^-1';   label.w_V = 'Molar weight of structure';

par.w_EC = (par.n_CEC * par.w_C) + (par.n_HEC * par.w_H) + (par.n_OEC * par.w_O)...
    + (par.n_NEC * par.w_N) + (par.n_PEC * par.w_P) ; 
free.w_EC = 0;   units.w_EC = 'g C mol C^-1';  label.w_EC = 'Molar weight of C reserve';


par.w_EN = (par.n_CEN * par.w_C) + (par.n_HEN * par.w_H) + (par.n_OEN * par.w_O) ...
    + (par.n_NEN * par.w_N) + (par.n_PEN * par.w_P) ; 
free.w_EN = 0;   units.w_EN = 'g N mol N^-1';  label.w_EN = 'Molar weight of N reserve';

par.w_EP = (par.n_CEP * par.w_C) + (par.n_HEP * par.w_H) + (par.n_OEP * par.w_O) ...
    + (par.n_NEP * par.w_N) + (par.n_PEP * par.w_P) ; 
free.w_EP = 0;   units.w_EP= 'g P mol N^-1';  label.w_EP = 'Molar weight of P reserve';

par.w_O2 = par.w_O * 2; 
free.w_O2 = 0;   units.w_O2 = 'g O2 mol O2^-1';   label.w_O2 = 'Molar weight of O2';

 

end