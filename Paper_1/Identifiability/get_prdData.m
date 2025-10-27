% cd /home/LAGUNES/Documents/GitHub/SargaDEB_working/matlab/multi_DEB/Paper_1

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




%%
[nm nst] = fieldnmnst_st(x_data_original);

for i = 1:nst
    st_name = nm{i}; %Name of the field of data
    varData = x_data_original.(st_name);
    k = size(varData, 2);
        if k == 1
            x_data_original.(st_name) = y_data_original.(st_name);
        elseif k == 2
            x_data_original.(st_name)(:,2) = y_data_original.(st_name);
        elseif k == 3
            x_data_original.(st_name)(:,2:3) = y_data_original.(st_name);
        end
end

dataAvailable = [x_data_original];
%% Save the available pseudo observations
save('Identifiability/ExperimentalCreation/Experiments/dataAvailable', 'dataAvailable')
%% Plots the available datasets
results_pets_newColors(par, metaPar, txtPar, data, auxData, metaData, txtData, weights)


%% Plots the available pseudo observations 
% If ran in the Identifaibility binder it will also plot the new pseudo observations
data2.Sargassum_fluitans = dataAvailable;

results_pets_newColors(par, metaPar, txtPar,data2, auxData, metaData, txtData, weights)

% results_pets_nonlydata(par, metaPar, txtPar,data2, auxData, metaData, txtData, weights)
