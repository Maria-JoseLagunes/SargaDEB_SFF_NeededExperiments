%% plot_SSE_variation
% Create plot of identifiability analysis based on datasets evaluated for
% each test and for each parameter evaluated


function plot_SSE_variation(SSE_structure, which_data, list_datasets, plotIndex, varFactor_list, varFactorIndex)
%% Description
%  Create plot of identifiability analysis based on datasets evaluated for
% each test and for each parameter evaluated
%
% Input
%
% * SSE_structure: structure with the SSE calculation for the original
% parameters, positive parameters and negative parameters
% * which_data: n-vector with numbers according to data selected 
% to estimate parameters based on this data 
% * plotIndex : vector to create a new x-axis entry for each set of
% datasets evaluated
% Output
% plot with n subplots for each parameter estimated and for each set of
% datasets evaluated
% 
%% Calls
%
% Called by : 
% * run_ID_analysis
%

global pets
% Determine the number of datasets
numDatasets = length(which_data);

% Load data and initial parameters
% [data, auxData, metaData, txtData, weights] = feval(['mydata_', pets{1}]);
% fieldsData = fieldnames(data);
% [par, metaPar, txtPar] = feval(['pars_init_',  pets{1}], metaData);
% [prdData0, info] = feval(['predict_', pets{1}], par, data, auxData);
% filternm = ['filter_', metaPar.model];

% Identify free = 1  parameters
% paramNames = fieldnames(par.free);
% freeParams = paramNames(cellfun(@(p) par.free.(p) == 1, paramNames));

freeParams = fieldnames(SSE_structure); 
numParams = length(freeParams);





% Define x-axis categories

selected_datasets = list_datasets(which_data);


%% Plot + setting up plot labels

numCols = ceil(sqrt(numParams));
numRows = ceil(numParams / numCols);
% label_str = selected_datasets(end); % use last dataset name as label
label_str = "Simu"; % use last dataset name as label


maxLim = max(struct2array(SSE_structure));
minLim = min(struct2array(SSE_structure));

varFactor = varFactor_list(varFactorIndex); 


numCols = 4;
numRows = 2;

for x = 1:numParams
    subplot(numCols,numRows,x);
    hold on;
    nameParam = freeParams{x};


    xVal = categorical(label_str + string(plotIndex));
switch varFactorIndex
    case 1
        colorPointGreen =  [179, 234, 159]; 
        colorPointBlue =  [128, 148, 227]; 
        colorPointRed =  [239, 134, 129]; 
    case 2
        colorPointGreen = [53, 198, 0 ]; 
        colorPointBlue =  [0, 40, 198]; 
        colorPointRed =  [223, 12, 2]; 
    case 3
        colorPointGreen = [20, 74, 0]; 
        colorPointBlue =  [36, 56, 134]; 
        colorPointRed =  [111, 6, 1]; 
end
        if numel(SSE_structure) == 4
            % Plot neutral variations (red)
            % scatter(xVal, SSE_structure(2).(nameParam), 75,  colorPointGreen./255, 'filled');
        
            % Plot positive variations
            colorPointBlue =  [    52   157    115]; %green
            scatter(xVal, SSE_structure(3).(nameParam), 75,  colorPointBlue./255, 'filled');
        
            % Plot negative variations
            colorPointRed =  [  247  96   54]; %orange
            scatter(xVal, SSE_structure(4).(nameParam), 75,  colorPointRed./255, 'filled');
            % set(gca, 'XTickLabel', label_str);
        else
    
            % Plot positive variations (red)
            scatter(xVal, SSE_structure(2).(nameParam), 100, 'r', 'filled');
        
            % Plot negative variations (blue)
            scatter(xVal, SSE_structure(3).(nameParam), 100, 'b', 'filled');
        

    end

    % set(gca, 'XTickLabel', label_str);
    % Formatting
    % ylim([minLim-0.1 maxLim+0.1]);
    yline(0, 'k--'); % Horizontal line at y=0
    % xlabel('Dataset');
    ylabel('SSE variation');
    
    nameParam_format = split(nameParam, '_');
    if ~contains(nameParam_format,'nContNum')
        latexnameParam = sprintf([nameParam_format{1}, '_{', nameParam_format{2}, '}']);
    else
        latexnameParam = nameParam_format;
    end
    title(latexnameParam) % 
   

end
