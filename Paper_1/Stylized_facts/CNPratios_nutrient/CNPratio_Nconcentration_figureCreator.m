%  N  CN ratio 
clc; clear; close all; 
load("simu.mat")
load("obs.mat")

figure(1)
for i = 1:length(obs)
    lightTested = categorical((simu(i).N)/1e-6); 

    % yyaxis left
     p1 = plot(lightTested, obs(i).CNendRatio); hold on; 
        p1.Color = simu(2).col;
        p1.LineStyle = '--';
        p1.Marker = 'o';
        p1.MarkerFaceColor = simu(2).col;
        p1.MarkerSize = 15;
    




        xlabel("N concentration (micro mol N  L^{-1})");
        ylabel('CN ratio (mol C mol N^{-1})');

     ax = gca;
      ax.XAxis.FontSize = 24;
    ax.YAxis.FontSize = 24;
     ax.YAxis.FontSize = 24;

   
    
end

%% Export figure

figName = '../../Figures/CNratio_Nconcentration';  % relative path to the Figures folder
set(gcf,'PaperOrientation','landscape');
% plotting in A4 format and positioning the plot in the PDF
set(gcf,'PaperUnits','centimeters'); 
set(gcf,'PaperType','A4');
set(gcf,'PaperPosition',[1. 1. 22. 14.]);
% print the plot as a vector in a PDF (if '-vector' is not working, try '-painters')
print(figName,'-vector','-bestfit','-dpdf')



