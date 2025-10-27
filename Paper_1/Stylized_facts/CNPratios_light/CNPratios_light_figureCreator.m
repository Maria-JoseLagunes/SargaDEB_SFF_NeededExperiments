% Light intensity, N % C % CN ratio (doublings day-1)
clc; clear; close all; 
load("simu.mat")
load("obs.mat")

figure(1)
for i = 1:length(obs)
    lightTested = categorical((simu(i).lightIntensity)/1e-6/3600); 

    % yyaxis left
    p1 = plot(lightTested, obs(i).CNendRatio); hold on; 
        p1.Color = simu(3).col;
        p1.LineStyle = '--';
        p1.Marker = 'o';
        p1.MarkerFaceColor = simu(3).col;
        p1.MarkerSize = 7.5; 
    




        xlabel("Light intensity (micro mol m^{-2} s^{-1})");
        ylabel('CN ratio (mol C mol N^{-1})');

     ax = gca;
     ax.YAxis.FontSize = 12;

    % 
    %      p2 = plot(lightTested, obs(i).percC_gdW(end)); hold on; 
    %     p2.Color = simu(i).col;
    %     p2.LineStyle = '--';
    %     p2.Marker = 'square';
    % 
    % p3 = plot(lightTested, obs(i).percN_gdW(end)); hold on; 
    %     p3.Color = simu(i).col;
    %     p3.LineStyle = '--';
    %     p3.Marker = 'square';
    % 
    % ax = gca;
    % ax.XAxis.FontSize = 12;
    % ax.YAxis(2).FontSize = 12;
    
end

%% Export figure

figName = '../../Figures/CNratio_lightIntensity';  % relative path to the Figures folder
set(gcf,'PaperOrientation','landscape');
% plotting in A4 format and positioning the plot in the PDF
set(gcf,'PaperUnits','centimeters'); 
set(gcf,'PaperType','A4');
set(gcf,'PaperPosition',[1. 1. 22. 14.]);
% print the plot as a vector in a PDF (if '-vector' is not working, try '-painters')
print(figName,'-vector','-bestfit','-dpdf')




%% Tash opttion
% 
% 
% figure(1)
% for i = 1:length(obs)
%     lightTested = categorical((simu(i).lightIntensity)/1e-6/3600); 
%     CNendRatio_vector = obs(i).CNendRatio; 
%     colorPlot_vector = simu(i).col; 
% 
% 
%     yyaxis left
%     p1 = plot(lightTested, obs(i).CNendRatio); hold on; 
%         p1.Color = simu(i).col;
%         p1.LineStyle = '--';
%         p1.Marker = 'o';
% 
%     xlabel("Temperature (°C)");
%     ylabel('Growth rate (doublings day^{-1})');
% 
% 
%     ax = gca;
%     ax.XAxis.FontSize = 12;
%     ax.YAxis(1).FontSize = 12;
% 
% 
%     yyaxis right
% 
%     p2 = plot(lightTested, obs(i).percC_gdW(end)); hold on; 
%         p2.Color = simu(i).col;
%         p2.LineStyle = '--';
%         p2.Marker = 'square';
% 
% 
%     p3 = plot(lightTested, obs(i).percN_gdW(end)); hold on; 
%         p3.Color = simu(i).col;
%         p3.LineStyle = '--';
%         p3.Marker = '^';
% 
% 
%     ylabel('%Nutrient gdW^{-1})');
%      ax.YAxis(2).FontSize = 12;
% 
% 
% end