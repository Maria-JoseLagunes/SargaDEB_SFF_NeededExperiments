
n = 5;% try 5, 7, 9, 12...
cmap = get_idf_palette(n);

% Visualize
figure;
imagesc(1:n);
colormap(cmap);
colorbar;
title(sprintf('%d-color palette with interpolation between base colors', n));

