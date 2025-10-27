% ---------------------------------------------------------------------
%% Temperature function analyses
% Maria José Lagunes
% First created : 2024/05/31
% ---------------------------------------------------------------------

%% Jouanno temperature dependent function
% parameters 
T_L_list = [24 10 20]; %°C, lower temperature limit below which growth ceases
T_H_list = [28 40 32]; %°C, higher temperature limit below which growth ceases
T_opt_list = [27.5 26 25]; %°C, optimum temperature at which growth ceases and T_ref
steepness_list = [2 1/2 1]; 
T_min = 10; T_max = 40;  
T = T_min:0.05:T_max; %°C, temperatures tested

T_A = 4000; 
T_AH = 50000;
T_AL = 40000;
T_ref = 25;


colors_pal = [          0.2980    0.0314    0.8314;
                      0    0.5020    0.1333;
                    0.6784    0.0549    0.0549

                       ] ;




% temperature functions 

for i = 1:length(T_L_list)
    T_L = T_L_list(i); %°C, lower temperature limit below which growth ceases
    T_H = T_H_list(i); %°C, higher temperature limit below which growth ceases
    T_opt = T_opt_list(i); %°C, optimum temperature at which growth ceases
    steepness = steepness_list(i);
    
    
    for l = 1:length(T)
        temp = T(l);
        if i ~= 3 
            fT = Jouanno(temp,T_L,T_H,T_opt,steepness);
            temperatureEffect(l,i) = fT;
            
        elseif i==3
            ct = arrhenius(temp, T_L, T_H, T_opt, T_A,T_AL,T_AH);
            temperatureEffect(l,i) = ct;
            
        end
    end

   fig = figure(1); hold on;
   plot(T,temperatureEffect(:,i), 'LineWidth', 1.5,'Color', colors_pal(i,:));


   
    % 
end
%

legend(["Witold", "Jouanno", "Arrhenius", "MG2023a", "", "MG2023b", "", "SLL2024", ""])
%% Export graphics
u = 1; 
%Saving parameters used
simulationParameters. T = T; 
simulationParameters. T_A = T_A; 
simulationParameters. T_ref = T_ref; 
simulationParameters. T_opt = T_opt;
simulationParameters. T_AH = T_AH; 
simulationParameters. T_AL = T_AL; 
simulationParameters. T_H = T_H; 
simulationParameters. T_L = T_L;

%Get today's date and time in MATLAB's serial date number format
% currentDate = now;
% 
% % Convert the serial date number to a formatted date string
% formattedDate = datestr(currentDate, 'ddmmyyyy');
% 
% figName = strcat(string(u),'temperatureEffect_',formattedDate,".png");
% directoryPath = '/home/wpodlejski/Documents/GitHub/MODEL_std/matlab/multi_DEB/Algae_test_forotherforcing/Temperature_analyses'; 
% parametersName=strcat(string(u),'simulationParameters_', formattedDate,'.mat');
% 
% % Save each figure separately
% exportgraphics(figure(1), fullfile(directoryPath,figName), 'Resolution',750); 
% save(fullfile(directoryPath, parametersName), 'simulationParameters');

currentDate = now;

% Convert the serial date number to a formatted date string
formattedDate = datestr(currentDate, 'ddmmyyyy');

figName = strcat("Temperatureequations","pdf");
directoryPath =  'C:\Users\Majo\Documents\PhD\SargaDEB_home\matlab\multi_DEB\Algae_test_forotherforcing'
folderName = strcat("1rstCSI",'_',formattedDate); 
mkdir(fullfile(directoryPath,folderName))


set(gcf,'PaperOrientation','landscape');
% plotting in A4 format and positioning the plot in the PDF
set(gcf,'PaperUnits','centimeters'); set(gcf,'PaperType','A4'); set(gcf,'PaperPosition',[1. 1. 22. 14.]);
% print the plot as a vector in a PDF (if '-vector' is not working, try '-painters')
print('temperatureEffect.pdf','-vector','-bestfit','-dpdf')


