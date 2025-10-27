%%
global pets pars_init_method n_pets results_output

% Select pars_init_method 1 to obtain initial parameters
% Select pars_init_method 2 to obtain parameters from .mat from estimated
% parameters
% pets == in this case only equal to Sargassum_fluitans
%
% 
pars_init_method = 1;


pets = {'Sargassum_fluitans'};
n_pets = length(pets);
pars_initnm = ['pars_init_', pets{1}];
resultsnm   = ['results_', pets{1}, '.mat'];
auxDatanm = ['auxData.', pets{1}];


if pars_init_method == 1
   [par, metaPar, txtPar] = feval(pars_initnm, []);
elseif pars_init_method == 2
   load(resultsnm, 'par');
   load(resultsnm, 'metaPar');
   load(resultsnm, 'txtPar');
end

[data, auxData, metaData, txtData, weights] = mydata_pets;
results_output = 3;

[prdData, dataVec, data2plot, prdData_x] = obtainingPredictData(par, metaPar, txtPar, data, auxData, metaData, txtData, weights); %function to obatin data from predictions
%With 100 points
x_data = data2plot.(pets{1}); %100 points of xdata
y_data = prdData_x.(pets{1}); %Obtaining structure with each data set
fldData = fieldnames(y_data);


% Original data
x_data_original  =  data.(pets{1}); %x and y original data
y_data_original = prdData.(pets{1}); %ydata of the predicted data with param values
colors = get_idf_palette(5); 

for i=1:length(fldData)
    figure(i)
    plot(x_data.(fldData{1}), y_data.(fldData{i}), '-', ...
        'Color', colors(2,:), 'LineWidth', 5)
    hold on; 
    scatter(x_data_original.(fldData{1})(:,1), y_data_original.(fldData{i}), ...
             150, ...                               
            'filled', ...
            'MarkerFaceColor', colors(5, :));              
    ylabel(fldData{i}, 'FontSize',22)
end


figure(1)
plot(x_data.(fldData{1}), y_data.(fldData{1}), '-', ...
    'Color', colors(2,:), 'LineWidth', 5)
hold on; 
scatter(x_data_original.(fldData{1})(:,1), x_data_original.(fldData{1})(:,2), ...
         150, ...                               
        'filled', ...
        'MarkerFaceColor', colors(5, :));              
ylabel(fldData{1}, 'FontSize',22)
%%
saveDir = pwd; 

for i = 1:length(fldData)
    figure(i);  
    

    timeStamp = char(datetime('now','Format','dd-MMM-uuuu'));
    
  
    pdfName = sprintf('%s_%s.pdf', fldData{i}, timeStamp);
    pngName = sprintf('%s_%s.png', fldData{i}, timeStamp);
    
    % Set paper and figure properties
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperType','A4');
    set(gcf,'PaperPosition',[1 1 22 14]);
    
    % Save as PDF (vector)
    print(fullfile(saveDir, pdfName), '-vector', '-bestfit', '-dpdf');
    
    % Save as PNG (raster)
    exportgraphics(gcf, fullfile(saveDir, pngName));
    
    fprintf('Saved figure %d: %s\n', i, fldData{i});
end