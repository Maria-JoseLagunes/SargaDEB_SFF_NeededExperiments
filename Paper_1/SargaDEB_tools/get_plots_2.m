% ---------------------------------------------------------------------
%% Make plots
% Maria José Lagunes
% First created : 2024/03/15
% Input
% t : vector of time
% mECENV : 3 - vector of state variables
% obs : structure with observables
% Output
% plots
% ---------------------------------------------------------------------
function get_plots_2(t,mECENV,J_struct,obs, simulationConditions)
global pets
%% Plot of results
vars_pull(simulationConditions);

m_EC = mECENV(:,1); % mol C mol V-1
m_EN = mECENV(:,2); % mol N mol V-1
M_V = mECENV(:,3);  % mol V



% Convert forcing variables in known units
I = J_struct.I / 1e-6 / 3600 ; % micro mol E m-2 s-1
N = J_struct.N / 1e-6; % micro mol N L-1
P = J_struct.P / 1e-6; % micro mol P L-1
%% Plot forcing variables

%
t = t/24 ; %t in days 

%% Forcing variables

T_set = simulationConditions.temp;
figure(1);

col =  [0.6275    0.5333    0.8196]; 
% time (t), Temperature (T°C)
subplot(3,5,1);hold on;
plot(t, repelem(T_set-273.15,numel(t)), 'Color', col, 'LineWidth', 1.2); xlabel("Time (days)"); ylabel("Temperature (°C)" );
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ax.XLim = [t(1,:) t(end)];


% time (t), Irradiance (micro mol E m-2 s-1)
subplot(3,5,2);hold on;
plot(t, I, 'Color', col,'LineWidth', 1.2); xlabel("Time (days)"); ylabel("micro mol E m^{-2} s^{-1}");
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ax.XLim = [t(1,:) t(end)];


% time (t), CO2 concentration (mol CO2 L-1)
subplot(3,5,3);hold on;
plot(t, J_struct.CO2 ,'Color',col,'LineWidth', 1.2); xlabel("Time (days)"); ylabel("mol DIC L^{-1}");
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ax.XLim = [t(1,:) t(end)];

% time (t), N concentration (mol NO3 + NH4 L-1)
subplot(3,5,4);hold on;
plot(t, N ,'Color',col,'LineWidth', 1.2); xlabel("Time (days)"); ylabel("micro mol NO_3^- and NH_4^+ L-1");
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;

% time (t), P concentration (mol PO43- L-1)
subplot(3,5,5);hold on;
plot(t, P ,'Color',col,'LineWidth', 1.2); xlabel("Time (days)"); ylabel("micro mol PO_4^{3-} L-1");
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ax.XLim = [t(1,:) t(end)];
%% Plot state variables
% time (t), structure (mol V)
subplot(3,5,6);hold on;
plot(t, M_V,'Color',col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("M_V (mol V)");
ax = gca;
% ax.YLim = [min(M_V) max(M_V) ];
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ax.XLim = [t(1,:) t(end)];

% time (t), C reserve density (mol E_C mol V-1)
subplot(3,5,7);hold on;
plot(t, m_EC,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("m_{E_C} (mol C molV-1)" );
ax = gca;
% ax.YLim = [min(m_EC) max(m_EC) ];
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ax.XLim = [t(1,:) t(end)];

% time (t), N reserve density (mol E_N mol V-1)
subplot(3,5,8);hold on;
plot(t, m_EN,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("m_{E_N} (mol N molV-1)");
ax = gca;
% ax.YLim = [min(m_EN) max(m_EN) ];
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ax.XLim = [t(1,:) t(end)];

% time (t), C reserve  (mol E_C )
subplot(3,5,9);hold on;
plot(t, m_EC.*M_V,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("M_{E_C} (mol C)" );
ax = gca;
% ax.YLim = [min(m_EC.*M_V) max(m_EC.*M_V) ];
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ax.XLim = [t(1,:) t(end)];

% time (t), N reserve  (mol E_N )
subplot(3,5,10);hold on;
plot(t, m_EN.*M_V,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("M_{E_N}(mol N )");
ax = gca;
% ax.YLim = [min(m_EN.*M_V) max(m_EN.*M_V) ];
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ax.XLim = [t(1,:) t(end)];

%% Plot fluxes
% time, j_AEN and j_AEC (assimilation)
subplot(3,5,11); hold on;
plot(t, J_struct.j_EC_A,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("Carbon assimilation flux (mol C mol V^{-1} h^{-1})");
ax = gca;
% ax.YLim = [min(J_struct.j_EC_A) max(J_struct.j_EC_A) ];
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 8;
ax.XLim = [t(1,:) t(end)];

subplot(3,5,12); hold on;
plot(t, J_struct.j_EN_A,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("Nitrogen assimilation flux (mol N mol V^{-1} h^{-1})");
ax = gca;
% ax.YLim = [min(J_struct.j_EN_A) max(J_struct.j_EN_A) ];
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 8;
ax.XLim = [t(1,:) t(end)];

subplot(3,5,13); hold on;
plot(t, J_struct.j_VG,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("Specific growth flux  (h^{-1})");
ax = gca;
% ax.YLim = [min(J_struct.j_VG) max(J_struct.j_VG) ];
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ax.XLim = [t(1,:) t(end)];

subplot(3,5,14); hold on;
plot(t, J_struct.j_EC_R,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("Rejection rate C  (mol C mol V^{-1} h^{-1})");
ax = gca;
% ax.YLim = [min(J_struct.j_EC_R) max(J_struct.j_EC_R) ];
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 8;
ax.XLim = [t(1,:) t(end)];

subplot(3,5,15); hold on;
plot(t, J_struct.j_EN_R,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("Rejection rate N  (mol N mol V^{-1} h^{-1})");
ax = gca;
% ax.YLim = [min(J_struct.j_EN_R) max(J_struct.j_EN_R)];
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 8;
ax.XLim = [t(1,:) t(end)];


set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','centimeters'); 
set(gcf,'PaperType','A4');
set(gcf,'PaperPosition',[1. 1. 22. 14.]);
[~,D,X] = fileparts(pwd);
D = strcat(D,X); 
print(['Figure1', D, '.pdf'],'-vector','-bestfit','-dpdf')
print(['Figure1', D, 'png'],'-dpng', '-r600')



%% plot observables in time
figure(2)
% unpacking observable vector
vars_pull(obs)
% time (t), wet weight (g Ww)
subplot(2,2,1);hold on;
plot(t, Ww, 'Color',col, 'LineWidth', 1.2);
xlabel("Time (days)");
ylabel('Weight (g_{Ww})');
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;

% time (t),  O_2 production rate
subplot(2,2,2);hold on;
plot(t, j_O2, 'Color',col, 'LineWidth', 1.2);
xlabel("Time (days)");
ylabel('O_2 production rate (micro mol {O_2} g_{Wd}^{-1} h^{-1})');
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;

% % time (t), O_2 production rate
% subplot(2,2,3);hold on;
% plot(t, j_O2, 'Color',col,'LineWidth', 1.2);
% 
% xlabel("Time (days)");
% ylabel('O_2 production rate (mol {O_2} g_{Ww}^{-1} h^{-1})');
% ax = gca;
% ax.XAxis.FontSize = 12;
% ax.YAxis.FontSize = 12;


% time (t), O_2 production rate
subplot(2,2,3);hold on;
plot(t, CNtotalRatio, 'Color',col,'LineWidth', 1.2);

xlabel("Time (days)");
ylabel('C:N ratio (mol C / mol N) ');
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;




% temperature, Relative growth rate (doublings day-1)
subplot(2,2,4);hold on;
temperatureTested = categorical(K2C(temp)); 
bar(temperatureTested, r_D, 'FaceColor',col);

xlabel("Temperature (°C)");
ylabel('Growth rate (doublings day^{-1})');
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;


set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','centimeters'); 
set(gcf,'PaperType','A4');
set(gcf,'PaperPosition',[1. 1. 22. 14.]);
[~,D,X] = fileparts(pwd);
D = strcat(D,X); 
print(['Figure2', D, '.pdf'],'-vector','-bestfit','-dpdf')
print(['Figure2', D, 'png'],'-dpng', '-r600')


%% Nutrient composition 
str = split(string(pets));
str = strrep(str, '_', ' ');
sgtitle(['\it ' str]);

figure(3)
% time C% dW
subplot(3,2,1);hold on;
plot(t, percC_gdW,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("%C (gdW)");
% plot(t, percCV_gdW,'Color', 'r','LineWidth', 1.5); xlabel("Time (days)"); ylabel("%C_V (gdW)");
% plot(t, percCV_gdW_with_ash,'Color', 'r', 'Linestyle', ':', 'LineWidth', 1.5); xlabel("Time (days)"); ylabel("%P_V (gdW)");
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;

% time N% dW
subplot(3,2,3);hold on;
plot(t, percN_gdW,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("%N (gdW)");
% plot(t, percNV_gdW, 'Color', 'r','LineWidth', 1.5); xlabel("Time (days)"); ylabel("%N_V (gdW)"); 
% plot(t, percNV_gdW_with_ash,'Color', 'r', 'Linestyle', ':','LineWidth', 1.5); xlabel("Time (days)"); ylabel("%P_V (gdW)");
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;

% time P% dW
subplot(3,2,5);hold on;
plot(t, percP_gdW,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("%P (gdW)");
% plot(t, percPV_gdW,'Color', 'r','LineWidth', 1.5); xlabel("Time (days)"); ylabel("%P_V (gdW)");
% plot(t, percPV_gdW_with_ash,'Color', 'r', 'Linestyle', ':','LineWidth', 1.5); xlabel("Time (days)"); ylabel("%P_V (gdW)");
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;

% time CN ratio
subplot(3,2,2);hold on;
plot(t, CNtotalRatio,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("CN ratio");
% plot(t, CN_V,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("CN ratio");
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;

% time CP ratio
subplot(3,2,4);hold on;
plot(t, CPtotalRatio,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("CP ratio");
%plot(t, CP_V,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("CP ratio");
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;


% time NP ratio
subplot(3,2,6);hold on;
plot(t, NPtotalRatio,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("NP ratio");
%plot(t, NP_V,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("NP ratio");
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;

set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','centimeters'); 
set(gcf,'PaperType','A4');
set(gcf,'PaperPosition',[1. 1. 22. 14.]);
[~,D,X] = fileparts(pwd);
D = strcat(D,X); 
print(['Figure3', D, '.pdf'],'-vector','-bestfit','-dpdf')
print(['Figure3', D, 'png'],'-dpng', '-r600')

%% Analyse propriétés du modèle
figure(4)
subplot(3,2,1);hold on;
plot(t, percCV_gdW,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("%C (gdW)");
plot(t, percCEC_gdW,'Color', col, 'Linestyle', '--','LineWidth', 1.5); xlabel("Time (days)"); ylabel("%C (gdW)");
plot(t, percCEN_gdW,'Color', col,'Linestyle', ':','LineWidth', 1.5); xlabel("Time (days)"); ylabel("%C (gdW)");
% plot(t, percC_gdW,'Color', 'r','LineWidth', 1.5); xlabel("Time (days)"); ylabel("%C (gdW)");
% legend("%C_V", "%C_{EC}", "%C_{EN}","%C_total", 'Location', 'best');
legend("%C_V", "%C_{EC}", "%C_{EN}",'Location', 'best');
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;

% time N% dW
subplot(3,2,3);hold on;
plot(t, percNV_gdW,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("%N (gdW)");
plot(t, percNEC_gdW, 'Color', col ,'Linestyle', '--','LineWidth', 1.5); xlabel("Time (days)"); ylabel("%N (gdW)"); 
plot(t, percNEN_gdW,'Color', col, 'Linestyle', ':','LineWidth', 1.5); xlabel("Time (days)"); ylabel("%N (gdW)");
% plot(t, percN_gdW,'Color', 'r','LineWidth', 1.5); xlabel("Time (days)"); ylabel("%N (gdW)");
legend("%N_V", "%N_{EC}", "%N_{EN}", 'Location', 'best');
% legend("%N_V", "%N_{EC}", "%N_{EN}","%N_total", 'Location', 'best');
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
% time P% dW
subplot(3,2,5);hold on;
plot(t, percPV_gdW,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("%P (gdW)");
plot(t, percPEC_gdW, 'Color', col ,'Linestyle', '--','LineWidth', 1.5); xlabel("Time (days)"); ylabel("%P (gdW)"); 
plot(t, percPEN_gdW,'Color', col, 'Linestyle', ':','LineWidth', 1.5); xlabel("Time (days)"); ylabel("%P (gdW)");
% plot(t, percP_gdW,'Color', 'r','LineWidth', 1.5); xlabel("Time (days)"); ylabel("%P (gdW)");
% legend("%P_V", "%P_{EC}", "%P_{EN}","%P_total", 'Location', 'best');
legend("%P_V", "%P_{EC}", "%P_{EN}", 'Location', 'best');
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
legend()

% time CN ratio
subplot(3,2,2);hold on;
% plot(t, CNtotalRatio,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("CN ratio");
plot(t, CN_V,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("CN ratio");
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;

% time CP ratio
subplot(3,2,4);hold on;
% plot(t, CPtotalRatio,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("CP ratio");
plot(t, CP_V,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("CP ratio");
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;

% time NP ratio
subplot(3,2,6);hold on;
% plot(t, NPtotalRatio,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("NP ratio");
plot(t, NP_V,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("NP ratio");
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;

set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','centimeters'); 
set(gcf,'PaperType','A4');
set(gcf,'PaperPosition',[1. 1. 22. 14.]);
[~,D,X] = fileparts(pwd);
D = strcat(D,X); 
print(['Figure4', D, '.pdf'],'-vector','-bestfit','-dpdf')
print(['Figure4', D, 'png'],'-dpng', '-r600')



%%
figure(5)
%% Plot state variables
% time (t), structure (mol V)
subplot(2,2,1);hold on;
plot(t, M_V,'Color',col,'LineWidth', 1.); xlabel("Time (days)"); ylabel("M_V (mol V)");
ax = gca;
% ax.YLim = [min(M_V) max(M_V) ];
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ax.XLim = [t(1,:) t(end)];

subplot(2,2,2);hold on;

plot(t, J_struct.j_VG,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("Specific growth flux  (h^{-1})");
ax = gca;
% ax.YLim = [min(J_struct.j_VG) max(J_struct.j_VG) ];
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ax.XLim = [t(1,:) t(end)];

% time (t), C reserve density (mol E_C mol V-1)
subplot(2,2,3);hold on;
plot(t, m_EC,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("m_{E_C} (mol C molV-1)" );
ax = gca;
% ax.YLim = [min(m_EC) max(m_EC) ];
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ax.XLim = [t(1,:) t(end)];

% time (t), N reserve density (mol E_N mol V-1)
subplot(2,2,4);hold on;
plot(t, m_EN,'Color', col,'LineWidth', 1.5); xlabel("Time (days)"); ylabel("m_{E_N} (mol N molV-1)");
ax = gca;
% ax.YLim = [min(m_EN) max(m_EN) ];
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ax.XLim = [t(1,:) t(end)];


set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','centimeters'); 
set(gcf,'PaperType','A4');
set(gcf,'PaperPosition',[1. 1. 22. 14.]);
[~,D,X] = fileparts(pwd);
D = strcat(D,X); 
print(['Figure5', D, '.pdf'],'-vector','-bestfit','-dpdf')
print(['Figure5', D, 'png'],'-dpng', '-r600')



%% 










