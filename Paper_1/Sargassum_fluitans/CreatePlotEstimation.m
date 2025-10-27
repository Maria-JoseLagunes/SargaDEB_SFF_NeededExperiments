% Create subplot with figures
figure
t = tiledlayout(2,2,'Padding','none','TileSpacing','none');
labels = {'A)', 'B)', 'C)', 'D)'};

for i = 1:4
    nexttile
    im = imread(sprintf('results_Sargassum_fluitans_0%d.png', i));
    imshow(im)
    axis off
    % Add panel label (top-left corner)
    text(500, 95, labels{i}, 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k')
end



set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','centimeters'); 
set(gcf,'PaperType','A4');
set(gcf,'PaperPosition',[1. 1. 22. 14.]);
print('total_results.pdf','-vector','-bestfit','-dpdf')
print('total_results.png','-dpng', '-r600')