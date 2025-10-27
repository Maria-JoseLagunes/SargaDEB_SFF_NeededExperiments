% Create plots of SSE error calculation
function plot_SSE(which_data)
%%
global pets
%% Plot SSE

load('SSE_structure.mat')

list_datasets = ["HT2023","MG2023a","MG2023b","VE2023", "TC_Ww", "Nutrient_uptake"]; %List of existent datasets on mydata_pets
data_setweight_Estimate = list_datasets(which_data); %Select dataset to set weight = 0, not consider during estimation

%%
% Load data and initial parameters
[data, auxData, metaData, txtData, weights] = feval(['mydata_', pets{1}]); 
fieldsData = fieldnames(data); 
[par, metaPar, txtPar] = feval(['pars_init_',  pets{1}], metaData);
[prdData0, info] = feval(['predict_', pets{1}], par, data, auxData); 
filternm = ['filter_', metaPar.model];
% Identify free = 1  parameters
paramNames = fieldnames(par.free);
freeParams = paramNames(cellfun(@(p) par.free.(p) == 1, paramNames));
%%

for j = 1:length(freeParams)
    pName = freeParams{j};
    
    % Extract the values from each structure
    SSE_ori_vector.(pName) = SSE_structure(1).(pName);
    SSE_pos_vector.(pName) = SSE_structure(2).(pName);
    SSE_neg_vector.(pName) = SSE_structure(3).(pName);
end
%%

%% Convert into matrices for plotting
numDatasets = length(list_datasets);
numParams = length(freeParams);

SSE_storage.SSE_ori = NaN(numDatasets, numParams);
SSE_storage.SSE_pos = NaN(numDatasets, numParams);
SSE_storage.SSE_neg = NaN(numDatasets, numParams);

for j = 1:numParams
    pName = freeParams{j};
    SSE_storage.SSE_ori(:, j) = SSE_ori_vector.(pName)(:)'; % force column vector
    SSE_storage.SSE_pos(:, j) = SSE_pos_vector.(pName)(:)';
    SSE_storage.SSE_neg(:, j) = SSE_neg_vector.(pName)(:)';
end

%% Now plot

figure(1);
clf;

numCols = ceil(sqrt(numParams)); 
numRows = ceil(numParams / numCols);

for x = 1:numParams
    subplot(numCols,numRows,x);
    hold on; 
    categoricalDatasetNames = categorical(list_datasets);
    categoricalDatasetNames = reordercats(categoricalDatasetNames, list_datasets);

    % Plot positive and negative variations
    scatter(categoricalDatasetNames, SSE_storage.SSE_pos(:, x), 'r', 'filled');
    scatter(categoricalDatasetNames, SSE_storage.SSE_neg(:, x), 'b', 'filled');

   
    ylim([-1 1]);
    yline(0, 'k--'); 
    title(freeParams{x}, 'Interpreter', 'none');
    xlabel('Dataset');
    ylabel('SSE variation');
end
% 
% sgtitle('SSE Variations per Parameter');
% % legend('Original', 'Positive Variation', 'Negative Variation', 'Location', 'Best');
end
