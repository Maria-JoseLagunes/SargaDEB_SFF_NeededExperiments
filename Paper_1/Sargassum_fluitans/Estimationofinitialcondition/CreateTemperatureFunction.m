% ---------------------------------------------------------------------
%% Temperature function analyses
% Maria José Lagunes
% First created : 2024/05/31
% Last modified : 2025/04
% ---------------------------------------------------------------------

%% global pets
global pets
pets = {'Sargassum_fluitans'};

[data, auxData, metaData, txtData, weights] = mydata_pets;

% % Temperature low and high limit according to literature
T_L = 18; %°C, lower temperature limit below which growth ceases/decreases
T_H = 32; %°C, higher temperature limit below which growth ceases/decreases



T_min = 10; T_max = 40;  
T = T_min:2.6:T_max; %°C, temperatures tested, 2.6 to have maximum 12 points for estimation procedure of DEBtoolM

% Arrhenius values to create simulated cT points
T_A =  	3000; 
T_AH = 60000;
T_AL = 40000;


% Reference temperature at which we encounter Sargassum (modified from
% typical DEB 20, as 20°C is the lower range of temperatures experienced by
% Sargassum)
T_ref = 25;
T_opt = [27.5]; %°C, optimum temperature and T_ref
steepness = [2]; 

T_L_jouanno = 20; %°C, lower temperature limit below which growth ceases/decreases
T_H_jouanno = 31; %°C, higher temperature limit below which growth ceases/decreases

colors_pal = [          0.2980    0.0314    0.8314;
                      0    0.5020    0.1333;
                    0.6784    0.0549    0.0549] ;


% To obtain datasets from mydata_pet
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
temperatureEffect =  zeros(length(T),1);
temperatureEffect_Jouanno =  zeros(length(T),1);
%% Arrhenius temperature function 

pars_T = [T_A; C2K(T_L); C2K(T_H); T_AL; T_AH]; % K
% pars_T = [T_A; T_L; T_H; T_AL; T_AH]; % K

fig = figure(); hold on;
   
for l = 1:length(T)
        temp = T(l); %°C
        ct = tempcorr(C2K(temp), C2K(T_ref), pars_T);
        % ct = tempcorr(C2K(temp), T_ref, pars_T);
        temperatureEffect(l,1) = ct; % -
        fT = Jouanno(temp,T_L_jouanno,T_H_jouanno,T_opt,steepness);
        temperatureEffect_Jouanno(l,1) = fT;

end


% Create plot of datasets with GR
uniqueDatasets = unique(datasets); 
colors = flipud(get_idf_palette(length(uniqueDatasets)+2));
refTemperature_finder = find(dataTable.temperature == C2K(25));
normalizationRGR = dataTable.RGR(refTemperature_finder);


for r=1:length(uniqueDatasets)
    idx = strcmp(dataTable.datasets, uniqueDatasets{r});
    scatter(K2C(dataTable.temperature(idx)), dataTable.RGR(idx) / normalizationRGR , ...
        80, 'filled', 'MarkerFaceColor', colors(r, :)); hold on;
end

% for r=1:length(uniqueDatasets)
%     idx = strcmp(dataTable.datasets, uniqueDatasets{r});
%     scatter(K2C(dataTable.temperature(idx)), dataTable.RGR(idx)  , ...
%         80, 'filled', 'MarkerFaceColor', colors(r, :)); hold on;
% end
% xlabel("Temperature (°C)"); ylabel("Growth rate (doubling day-1)")
% legend(["HS1987", "MG2023a" , "MG2023b"], "Location",'best')

xlabel("Temperature (°C)"); ylabel("Normalized growth rate (-)")
legend(["HS1987", "MG2023a" , "MG2023b"], "Location",'best')

xlim([ 10    40]);
ylim([  0    1.4000]);








%% Create plot of function cT vector before estimation

% plot(T,temperatureEffect, 'r-', 'MarkerSize', 12); hold on;
% xlabel("Temperature (°C)"); ylabel("cT (-)")
% legend(["HS1987", "MG2023a" , "MG2023b","Arrhenius function"], "Location",'best')

plot(T,temperatureEffect, 'r-', 'MarkerSize', 12); hold on;
xlabel("Temperature (°C)"); ylabel("cT (-)")

% plot(T,temperatureEffect_Jouanno, 'b--', 'MarkerSize', 12); hold on;
% xlabel("Temperature (°C)"); ylabel("cT (-)")

% legend(["HS1987", "MG2023a" , "MG2023b","Arrhenius function", "Jouanno function"], "Location",'best')

%% Load vector used for estimation with modified data point
pets = {'Sargassum_fluitans_Arrhenius'};

[data_Arr, auxData_Arr, metaData_Arr, txtData_Arr, weights_Arr] = mydata_Sargassum_fluitans_Arrhenius;
plot(data_Arr.TC_RG(:,1), data_Arr.TC_RG(:,2), 'r*', 'MarkerSize', 12);

legend(["HS1987", "MG2023a" , "MG2023b","Arrhenius function", ...
    "Corrected Arrhenius function"], "Location",'best') 

%%
resultsnm   = ['results_', pets{1}, '.mat'];
load(resultsnm, 'par');
load(resultsnm, 'metaPar');
load(resultsnm, 'txtPar');

vars_pull(par)
pars_T   
for l = 1:length(T)
        temp = T(l); %°C
        ct = tempcorr(C2K(temp), C2K(T_ref), pars_T);
        % ct = tempcorr(C2K(temp), T_ref, pars_T);
        temperatureEffect(l,1) = ct; % -
        fT = Jouanno(temp,T_L_jouanno,T_H_jouanno,T_opt,steepness);
        temperatureEffect_Jouanno(l,1) = fT;

end


r = figure(1);
ax = findall(r, 'Type', 'axes');  % Find the axes inside the figure
x_limits = xlim(ax);
y_limits = ylim(ax);

xlim(x_limits);
ylim(y_limits);

%% Create vector of points to use as input in mydata file
vectorData= ([T' temperatureEffect]); 
save(fullfile(pwd,'1_Arrhenius_parameters/data_TC_GR'), "vectorData")
%% Export graphics
u = 1; 
%Saving parameters used
simulationParameters. T = T; 
simulationParameters. T_A = T_A; 
simulationParameters. T_ref = T_ref; 
simulationParameters. T_AH = T_AH; 
simulationParameters. T_AL = T_AL; 
simulationParameters. T_H = T_H; 
simulationParameters. T_L = T_L;

%Get today's date and time in MATLAB's serial date number format
timeStamp = char(datetime('now','Format','dd-MMM-uuuu'));
timeStamp = replace(timeStamp, ":", "_");


figName = strcat("ArrheniusEstimation", timeStamp,".pdf");
figName2 = strcat("ArrheniusEstimation", timeStamp, ".png");

saveDir   = ['1_Arrhenius_parameters/'];



set(gcf,'PaperOrientation','landscape');
% plotting in A4 format and positioning the plot in the PDF
set(gcf,'PaperUnits','centimeters'); set(gcf,'PaperType','A4'); set(gcf,'PaperPosition',[1. 1. 22. 14.]);
% print the plot as a vector in a PDF (if '-vector' is not working, try '-painters')
print(fullfile(saveDir,figName),'-vector','-bestfit','-dpdf')
exportgraphics(gcf, fullfile(saveDir, figName2) )

%%
save(fullfile(saveDir, 'simulationParameters.mat'), 'simulationParameters');
