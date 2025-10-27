%Work
addpath(genpath('/home/LAGUNES/Documents/GitHub/SargaDEB_working/DEBtool_M-master'))
cd /home/LAGUNES/Documents/GitHub/SargaDEB_working/matlab/multi_DEB/Paper_1
addpath(genpath("/home/LAGUNES/Documents/GitHub/SargaDEB_working/matlab/multi_DEB/Paper_1/"))

%Laure
  % addpath(genpath('/Users/laurepecquerie/Documents/GitHub/SargaDEB_working/DEBtool_M-master'))
  % cd /Users/laurepecquerie/Documents/GitHub/SargaDEB_working/matlab/multi_DEB/Paper_1/
  % addpath(genpath("/Users/laurepecquerie/Documents/GitHub/SargaDEB_working/matlab/multi_DEB/Paper_1"))
 % 
%Home
% addpath(genpath('C:\Users\Majo\Documents\PhD\SargaDEB_home\DEBtool_M-master'))
% cd C:\Users\Majo\Documents\PhD\SargaDEB_home\matlab\multi_DEB\Cleaningcodes
% addpath("C:\Users\Majo\Documents\PhD\SargaDEB_home\matlab\multi_DEB\Cleaningcodes")


cd("Identifiability/ExperimentalCreation")
%% Selection of datasets not to estimate
clc; clear; close all;
%% Datasets available for estimation
% The PE refers to : pseudo numerical experiment
%
% 1 : HT2023
% 2 : HS1987
% 3 : MG2023a
% 4 : MG2023b
% 5 : VE2023
% 6 : PE : TC_Ww temperatures (10:3:40) with datapoints at [0 2 4 6 8 10] days 
% 7 : PE : N uptake at three temperatures
% 8 : PE : TC_Ww temperatures (10:3:40) with datapoints at [0 10] days 
% 9 : PE : Starvation with datapoints at [0 10 20 30] days Ww and CNtotalratio at 25°C

%% Estimation
% Selection of which data sets to use for the estimation procedure. Those
% not included will be set with weight = 0 in the EstimationProcedure_NE
% function
global pets

list_datasets = ["HT2023", "HS1987", "MG2023a","MG2023b","VE2023", "TC_Ww",...
    "N_uptake", "TC_WW_Twodays", 'Starvation']; %List of existent datasets on mydata_pets



pets={'Sargassum_fluitans'};
%% Analysis if another date to retrieve
clf


    all_which_data = {

    [1 3 4]
    [1 3 4 5]

    [1 3 4 5 7]
    [1 3 4 5 9]
    [1 3 4 5 7 9]

    [1 3 4 7 9]
    };

  

f = figure(2); % Same figure to accumulate plots
hold on; % Keep all the previous plots

% timeStamp = char(datetime('now','Format','dd_MM_uuuu'));
timeStamp = char('22_07_2025');
% timeStamp = char('18_09_2025');
% timeStamp = char('23_09_2025');
numParams = 8; 

varFactor_list = [0.05 0.1 0.2] ; % Variation factor, 10%
varFactor_list = [ 0.1 ] ; % Variation factor, 10%

for j =1:length(varFactor_list)
varFactor = varFactor_list(j); % Variation factor, 10%

    for i = 1:length(all_which_data)
        which_data   = all_which_data{i};
        config_label = strjoin(string(which_data), '_');
        config_label = string(numParams) + '_' + config_label ;
    
        filename = strjoin(['SSE_Analysis/SSE_results', config_label, 'day', timeStamp, ...
        replace(string(varFactor), '.', '_'), '.mat'], '_'); 
        if isfile(filename)
            load(filename); % Should load SSE_results
            plot_SSE_variation(SSE_results.SSE_structure, which_data, list_datasets, i, varFactor_list, j);
        else
            disp('no file with this name')
        end
    end 
end

%% Save figure
timeStamp = "24102025"; 
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','centimeters'); 
set(gcf,'PaperType','A4');
set(gcf,'PaperPosition',[3 5 15 20]);
set(gcf,'PaperSize',[21,20]);
figName = strcat('idf_analsysis_diffscale', timeStamp , '.pdf'); 
figName_png = strcat('idf_analsysis_diffscale',timeStamp, '.png'); 

print(figName,'-vector','-bestfit','-dpdf')
print(figName_png,'-dpng', '-r600')
%%
fig1 = gcf;
fig2 = copyobj(fig1, 0); 
figure(fig2);            

% Get all axes in the current figure
ax = findall(fig2, 'Type', 'axes');


ylims = zeros(length(ax), 2);

for i = 1:length(ax)
    ylims(i,:) = ylim(ax(i)); 
end

ymin = min(ylims(:,1));
ymax = max(ylims(:,2));

% Compute range and margin (5%)
yrange = ymax - ymin;
margin = 0.05 * yrange;
%%
for i = 1:length(ax)
    ylim(ax(i), [ymin-margin ymax+margin]);
    ylim(ax(i), [-0.15 ymax+margin]);
end
%%
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','centimeters'); 
set(gcf,'PaperType','A4');
set(gcf,'PaperPosition',[3 5 15 20]);
set(gcf,'PaperSize',[21,20]);
figName = strcat('idf_analsysis_samescale', timeStamp , '.pdf'); 
figName_png = strcat('idf_analsysis_samescale',timeStamp, '.png'); 

print(figName,'-vector','-bestfit','-dpdf')
print(figName_png,'-dpng', '-r600')