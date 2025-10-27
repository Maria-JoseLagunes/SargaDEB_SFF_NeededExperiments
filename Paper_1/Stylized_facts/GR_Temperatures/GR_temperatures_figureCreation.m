% Temperature, Relative growth rate (doublings day-1)
clc; clear; close all; 
load("simu.mat")
load("obs.mat")

figure(1)
for i = 1:length(obs)
    sz = 100; 
    temperatureTested = categorical(K2C(simu(i).temp)); 

    % scatter(temperatureTested, obs(i).r_D, sz, 'MarkerFaceColor', simu(4).col, ...
    %     'MarkerEdgeColor', simu(4).col); hold on; 
    
    p1 = plot(temperatureTested, obs(i).r_D); hold on; 
        p1.Color = simu(2).col;
        p1.LineStyle = '--';
        p1.Marker = 'o';
        p1.MarkerFaceColor = simu(2).col;
        p1.MarkerSize = 15; 
    xlabel("Temperature (°C)");
    ylabel('Growth rate (doublings day^{-1})');
    ax = gca;
    ax.XAxis.FontSize = 24;
    ax.YAxis.FontSize = 24;
    
end



%% Export figure

figName = '../../Figures/GR_temperature';  % relative path to the Figures folder
set(gcf,'PaperOrientation','landscape');
% plotting in A4 format and positioning the plot in the PDF
set(gcf,'PaperUnits','centimeters'); 
set(gcf,'PaperType','A4');
set(gcf,'PaperPosition',[1. 1. 22. 14.]);
% print the plot as a vector in a PDF (if '-vector' is not working, try '-painters')
print(figName,'-vector','-bestfit','-dpdf')

