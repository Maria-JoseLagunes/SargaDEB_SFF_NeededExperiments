% Initialize time, parameters, initial conditions, ...
pets = {'Sargassum_fluitans'};
pars_init_method = 2; %1 = from initial parameters %2 from .mat estimation of parameters



pars_initnm = ['pars_init_', pets{1}];
resultsnm   = ['results_', pets{1}, '.mat'];
auxDatanm = ['auxData.', pets{1}];


if pars_init_method == 1
   [par, metaPar, txtPar] = feval(pars_initnm, []);
elseif pars_init_method == 2
   load(resultsnm, 'par');
end


vars_pull(par)


%compute temperature correction factors
pars_T = [T_A; T_L; T_H; T_AL; T_AH]; % K
TC_HL = tempcorr(C2K(23), T_ref,pars_T); % - 
TC_LL = tempcorr(C2K(23), T_ref,pars_T); % -

%Correction of parameter values for Vazquez-Elizondo data
k_I_CT_HL = k_I * TC_HL; % mol γ mol PSU–1 h–1
k_I_CT_LL = k_I * TC_LL; % mol γ mol PSU–1 h–1

%Transform parameters with z factor algal
dwratio_SF = dwratio * zalgal ; % gWd/gWw, Dry to wet (gdw/gww) weight ratio
smcoeff_SF = smcoeff / zalgal;  % g cm-2, Coefficient of linear regression from Chambon


I_light = 0:1500; % micro mol m-2 s-1
I = I_light * 3600 * 1e-6; % mol m-2 h-1


% VE2023
% Photosynthesis - Specific relaxation rate
Wd = 1.5 * smcoeff_SF * dwratio_SF;  % g dW SM = g cm-2 1.5 in cm-2 from VE2023, same for HL and LL
M_V= Wd / (w_V+   (m_EN_0*w_EN) +   (m_EC_0*w_EC)); %mol V

rho_PSU_HL =  0.5;

J_I_HL = (rho_PSU_HL * I * eps_I) ./ ...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
    (1 +  (I * eps_I) / k_I_CT_HL); %% mol gamma mol V-1 h-1

J_O2_HL =  J_I_HL .* (M_V / Wd)  * y_LO2 * smcoeff_SF * dwratio_SF / 1e-6; %micro mol 02 cm-2  h-1


rho_PSU_LL =  0.17;

J_I_LL = (rho_PSU_LL* I* eps_I) ./ ...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
    (1 +  (I* eps_I) / k_I_CT_LL) ; %% mol gamma mol V-1 h-1

J_O2_LL =  J_I_LL .* (M_V / Wd)  * y_LO2 * smcoeff_SF * dwratio_SF / 1e-6; % micro mol 02 cm-2  h-1





%% Create figure and save
figure()
plot(I_light, J_O2_HL, 'b', 'LineWidth', 4.5); hold on; 
% plot(I_light, J_O2_LL, 'b', 'LineWidth', 1.5); hold on; 
xlabel('Irradiance (micro mol m^{-2}  s^{-1})'); ylabel('Photosynthetic rate (micro mol 0_2 cm^{-2}  h^{-1})');
% leg_tex = legend('HL acclimated ($\rho_{psu}$ = 0.5)', 'LL acclimated ($\rho_{psu}$ = 0.17)'); 
% set(leg_tex,'Interpreter','latex');
% set(leg_tex,'FontSize',12);
 ax = gca;
    ax.XAxis.FontSize = 24;
    ax.YAxis.FontSize = 24;
%% Save figure
figName = '../../Figures/light_O2production';  % relative path to the Figures folder
set(gcf,'PaperOrientation','landscape');
% plotting in A4 format and positioning the plot in the PDF
set(gcf,'PaperUnits','centimeters'); 
set(gcf,'PaperType','A4');
set(gcf,'PaperPosition',[1. 1. 22. 14.]);
% print the plot as a vector in a PDF (if '-vector' is not working, try '-painters')
print(figName,'-vector','-bestfit','-dpdf')
