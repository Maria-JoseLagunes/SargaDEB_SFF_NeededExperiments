% ---------------------------------------------------------------------
%% Temperature function analyses
% Maria José Lagunes
% First created : 2024/05/31
% Last modified : 2025/04
% ---------------------------------------------------------------------

%% global pets
cd 
global pets
pets = {'Sargassum_fluitans'};

[data, auxData, metaData, txtData, weights] = mydata_pets;

% % Temperature low and high limit according to literaturec
T_L = 18; %°C, lower temperature limit below which growth ceases/decreases
T_H = 32; %°C, higher temperature limit below which growth ceases/decreases



T_min = 10; T_max = 40;  
T = T_min:2.6:T_max; %°C, temperatures tested, 2.6 to have maximum 12 points for estimation procedure of DEBtoolM
T = T_min:0.5:T_max; %°C, temperatures tested, 2.6 to have maximum 12 points for estimation procedure of DEBtoolM

% Arrhenius values to create simulated cT points
T_A =  	3000; 
T_AH = 60000;
T_AL = 40000;

pars_T_calibration_hand = [T_A; C2K(T_L); C2K(T_H); T_AL; T_AH]; % K

%% DEB estimation using the datasets 
% This are the resulting parameters, to replicate this go to file 
% Miscelaneous : ArrheniusCrudeEstimation

T_A = 5643;
T_L = 288.4;
T_H = 302.1;
T_AH = 3.456e+04;
T_AL = 7.982e+04;
pars_T_classical_estimation = [T_A; T_L; T_H; T_AL; T_AH]; % K

% Reference temperature at which we encounter Sargassum (modified from
% typical DEB 20, as 20°C is the lower range of temperatures experienced by
% Sargassum)

T_ref = 25;
T_opt = [27.5]; %°C, optimum temperature and T_ref
steepness = [2]; 

T_L_jouanno = 20; %°C, lower temperature limit below which growth ceases/decreases
T_H_jouanno = 31; %°C, higher temperature limit below which growth ceases/decreases

% colors_pal = [          0.2980    0.0314    0.8314;
%                       0    0.5020    0.1333;
%                     0.6784    0.0549    0.0549] ;
% 
% 
%% To obtain datasets from mydata_pet
dataSetsFields = fieldnames(auxData.(pets{1}).temp);
temperature = [];
time = [];
RGR = [];  
datasets = [];

for j = 1:length(dataSetsFields)
    fieldName = dataSetsFields{j};
    if contains(fieldName, 'tWw')
    datasets = [datasets; {txtData.(pets{1}).bibkey.(fieldName)}]; 
    temperature = [temperature; auxData.(pets{1}).temp.(fieldName)]; 
    time = [time; auxData.(pets{1}).time.(fieldName)]; 
    RGR = [RGR; log2(data.(pets{1}).(fieldName)(end,2) / data.(pets{1}).(fieldName)(1,2)) / (auxData.(pets{1}).time.(fieldName)) ]; 
    end
end

dataTable = table(datasets, temperature, RGR);

% Create plot of datasets with GR
uniqueDatasets = unique(datasets); 

colors_palette = flipud(get_idf_palette(length(uniqueDatasets)+5));

refTemperature_finder = find(dataTable.temperature == C2K(T_ref));
normalizationRGR = dataTable.RGR(refTemperature_finder);
dataForm = {'o', 's', '^'}; 

fig = figure(); hold on;
for r=1:length(uniqueDatasets)
    idx = strcmp(dataTable.datasets, uniqueDatasets{r});
     scatter(K2C(dataTable.temperature(idx)), ...
            dataTable.RGR(idx) / normalizationRGR, ...
            80, ...                               
            'filled', ...
            'MarkerFaceColor', colors_palette(r+1, :), ...
            'Marker', dataForm{r});              
end
set(gca, 'FontSize', 20)

xlabel("Temperature (°C)", 'FontSize', 22); ylabel("Normalized growth rate (-)", 'FontSize',22)
% legend(["HS1987", "MG2023a" , "MG2023b"], "Location",'best')

legend(["Hanisak and Samuel 1987", "Magaña-Gallegos et al 2023a" , "Magaña-Gallegos et al 2023b"], "Location",'best')
xlim([ 10    40]);     % hand choosen to create figure
ylim([  0    1.4000]); % hand choosen to create figure


%% Arrhenius temperature function 

% % pars_T = [T_A; C2K(T_L); C2K(T_H); T_AL; T_AH]; % K
% pars_T = [T_A; T_L; T_H; T_AL; T_AH]; % K

temperatureEffect_classical_estimation  =  zeros(length(T),1);
temperatureEffect_calibration_hand  =  zeros(length(T),1);
temperatureEffect_Jouanno =  zeros(length(T),1);


   
for l = 1:length(T)
        temp = T(l); %°C
        ct_classical_estimation = tempcorr(C2K(temp), C2K(T_ref), pars_T_classical_estimation);
        ct_calibration_hand = tempcorr(C2K(temp), C2K(T_ref), pars_T_calibration_hand);

        % ct = tempcorr(C2K(temp), T_ref, pars_T);
        temperatureEffect_classical_estimation(l,1) = ct_classical_estimation; % -
        temperatureEffect_calibration_hand(l,1) = ct_calibration_hand; % -

        fT = Jouanno(temp,T_L_jouanno,T_H_jouanno,T_opt,steepness);
        temperatureEffect_Jouanno(l,1) = fT;

end



% for r=1:length(uniqueDatasets)
%     idx = strcmp(dataTable.datasets, uniqueDatasets{r});
%     scatter(K2C(dataTable.temperature(idx)), dataTable.RGR(idx)  , ...
%         80, 'filled', 'MarkerFaceColor', colors(r, :)); hold on;
% end
% xlabel("Temperature (°C)"); ylabel("Growth rate (doubling day-1)")
% legend(["HS1987", "MG2023a" , "MG2023b"], "Location",'best')





%% Create plot of function cT vector before estimation
% colors = (get_idf_palette(8));
% plot(T,temperatureEffect, 'r-', 'MarkerSize', 12); hold on;
% xlabel("Temperature (°C)"); ylabel("cT (-)")
% legend(["HS1987", "MG2023a" , "MG2023b","Arrhenius function"], "Location",'best')

% plot(T,temperatureEffect_calibration_hand, '-', 'Color', colors_palette(5,:), ...
%     'MarkerSize', 12, 'LineWidth', 1.7); hold on;
% xlabel("Temperature (°C)"); ylabel("cT (-)")
% legend(["HS1987", "MG2023a" , "MG2023b","Arrhenius function (hand calibration)"], "Location",'best')
% 
% % plot(T,temperatureEffect_Jouanno, 'b--', 'MarkerSize', 12); hold on;
% % xlabel("Temperature (°C)"); ylabel("cT (-)")
% 
% plot(T,temperatureEffect_classical_estimation, '--', 'Color', colors_palette(6,:), ...
%     'MarkerSize', 12, 'LineWidth', 1.7); hold on;
% xlabel("Temperature (°C)"); ylabel("cT (-)")
% % legend(["HS1987", "MG2023a" , "MG2023b","Arrhenius function (classical estimation)"], "Location",'best')
% legend(["HS1987", "MG2023a" , "MG2023b","Arrhenius function (hand calibration)","Arrhenius function (classical estimation)"], "Location",'best')


% legend(["HS1987", "MG2023a" , "MG2023b","Arrhenius function", "Jouanno function"], "Location",'best')

%% Load vector used for estimation with modified data point
pets = {'Sargassum_fluitans_Arrhenius'};

[data_Arr, auxData_Arr, metaData_Arr, txtData_Arr, weights_Arr] = mydata_Sargassum_fluitans_Arrhenius;
% plot(data_Arr.TC_RG(:,1), data_Arr.TC_RG(:,2), 'b*', 'MarkerSize', 12);

scatter(data_Arr.TC_RG(:,1), data_Arr.TC_RG(:,2)  , ...
        80, 'MarkerEdgeColor', colors_palette(6, :), 'LineWidth', 2.2); hold on;
% 
% legend(["HS1987", "MG2023a" , "MG2023b","Arrhenius function", ...
%     "Corrected Arrhenius function"], "Location",'best') 


% legend(["HS1987", "MG2023a" , "MG2023b","Created Arrhenius vector"], "Location",'best')
legend(["Hanisak and Samuel 1987", "Magaña-Gallegos et al 2023a" , "Magaña-Gallegos et al 2023b","Created Arrhenius vector"], "Location",'best')


% legend(["HS1987", "MG2023a" , "MG2023b","Arrhenius function (hand calibration)",...
    % "Arrhenius function (classical estimation)", "Corrected Arrhenius function"], "Location",'best')



%% Plot estimation with vector created
resultsnm   = ['results_', pets{1}, '.mat'];
load(resultsnm, 'par');
load(resultsnm, 'metaPar');
load(resultsnm, 'txtPar');

vars_pull(par)
pars_T_estimated_vector =  [T_A; T_L; T_H; T_AL; T_AH]; % K
temperatureEffect_estimated_vector = zeros(length(T),1); 

for l = 1:length(T)
        temp = T(l); %°C
        ct_estimated_vector = tempcorr(C2K(temp), T_ref, pars_T_estimated_vector);
        % ct = tempcorr(C2K(temp), T_ref, pars_T);
        temperatureEffect_estimated_vector(l,1) = ct_estimated_vector; % -     
end

plot(T,temperatureEffect_estimated_vector, '-', 'Color', colors_palette(7,:), ...
     'LineWidth', 1.7); hold on;
xlabel("Temperature (°C)", 'FontSize', 22); ylabel("Normalized growth rate(-)", 'FontSize', 22)

% legend(["HS1987", "MG2023a" , "MG2023b","Created Arrhenius vector", "Arrhenius parameter estimation"], "Location",'best')
% legend(["Hanisak and Samuel 1987", "Magaña-Gallegos et al 2023a" , "Magaña-Gallegos et al 2023b","Created Arrhenius vector", "Arrhenius parameter estimation"], "Location",'best')
legend(["Hanisak and Samuel 1987", "Magaña-Gallegos et al 2023a" , "Magaña-Gallegos et al 2023b", "Arrhenius parameter estimation"], "Location",'best')

% legend(["HS1987", "MG2023a" , "MG2023b","Arrhenius function (hand calibration)",...
    % "Arrhenius function (classical estimation)", "Corrected Arrhenius function", ...
    % "Estimation corrected function"], "Location",'best')




r = figure(1);
ax = findall(r, 'Type', 'axes');  % Find the axes inside the figure
x_limits = xlim(ax);
y_limits = ylim(ax);
set(ax, 'FontSize', 22); 
xlim(x_limits);
ylim(y_limits);
%%



% %% Create vector of points to use as input in mydata file
% vectorData= ([T' temperatureEffect]); 
% save(fullfile(pwd,'1_Arrhenius_parameters/data_TC_GR'), "vectorData")
% %% Export graphics
% u = 1; 
% %Saving parameters used
% simulationParameters. T = T; 
% simulationParameters. T_A = T_A; 
% simulationParameters. T_ref = T_ref; 
% simulationParameters. T_AH = T_AH; 
% simulationParameters. T_AL = T_AL; 
% simulationParameters. T_H = T_H; 
% simulationParameters. T_L = T_L;
% 
% %%
saveDir   = ['Figure1/'];

%Get today's date and time in MATLAB's serial date number format
timeStamp = char(datetime('now','Format','dd-MMM-uuuu'));
timeStamp = replace(timeStamp, ":", "_");


figName = strcat("ArrheniusEstimation", timeStamp,".pdf");
figName2 = strcat("ArrheniusEstimation", timeStamp, ".png");

%% Save figure

set(gcf,'PaperOrientation','landscape');
% plotting in A4 format and positioning the plot in the PDF
set(gcf,'PaperUnits','centimeters'); set(gcf,'PaperType','A4'); set(gcf,'PaperPosition',[1. 1. 22. 14.]);
% print the plot as a vector in a PDF (if '-vector' is not working, try '-painters')
print(fullfile(saveDir,figName),'-vector','-bestfit','-dpdf')

figName = "ArrheniusEstimation.pdf";
print(fullfile(saveDir,figName),'-vector','-bestfit','-dpdf')

exportgraphics(gcf, fullfile(saveDir, figName2) )

%%
% save(fullfile(saveDir, 'simulationParameters.mat'), 'simulationParameters');
