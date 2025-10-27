%Simulate nutrient uptake
[simu, obs] = run_Nuptake(0:0.2:50, [22 26 28]);

%Plot nutrient uptake
figure();
colorMatrix = colormap(abyss(5));
for i = 1:length(simu)
    plot(simu(i).N_microConcentration, obs(i).j_EN_A_inmicroMol, 'LineWidth', 4, ...
        'Color',  colorMatrix(i+2,:)  ); hold on;
    legends{i} = (sprintf('T = %d °C',simu(i).temp)); hold on; 

end

ax = gca; 
set(ax,'FontSize',24)

xlabel('N concentration (micro mol N L^{-1})', 'FontSize', 24);
ylabel('N uptake rate (micro mol N g dW^{-1} h^{-1})', 'FontSize', 24);
legend(legends, 'Location', 'best');
%% Saving figure

figName = '../../Figures/Nconcentration_Nuptake';  % relative path to the Figures folder
set(gcf,'PaperOrientation','landscape');
% plotting in A4 format and positioning the plot in the PDF
set(gcf,'PaperUnits','centimeters'); 
set(gcf,'PaperType','A4');
set(gcf,'PaperPosition',[1. 1. 22. 14.]);
% print the plot as a vector in a PDF (if '-vector' is not working, try '-painters')
print(figName,'-vector','-bestfit','-dpdf')
