function st = obtainingDatasetPrediction(selectDataset, nm_points)
global pets

n_pets = length(pets);
pars_initnm = ['pars_init_', pets{1}];
resultsnm   = ['results_', pets{1}, '.mat'];
auxDatanm = ['auxData.', pets{1}];



[par, metaPar, txtPar] = feval(pars_initnm, []);

[data, auxData, metaData, txtData, weights] = mydata_pets;


[data2plot, prdData_x] = obtainingPredictData(par, metaPar, txtPar, data, auxData, metaData, txtData, weights); %function to obatin data from predictions
x_data = data2plot.(pets{1}); 
y_data = prdData_x.(pets{1}); %Obtaining structure with each data set
fldData = fieldnames(y_data);

st = struct(); % To save values of predicted data

for i = 1:length(fldData)
    fldnm = fldData{i};
    if contains(fldnm, selectDataset)
    x_dataset = x_data.(fldnm);
    y_dataset = y_data.(fldnm);
    k = size(y_dataset, 1);
    if k==1   
        continue;
    else
        % Select points 
        idx_selectPoints = round(linspace(1,length(y_dataset), nm_points)); %Selecting points based on a linspace
        x_select = x_dataset(idx_selectPoints);
        y_select = y_dataset(idx_selectPoints); %Selecting points from dataset
    end
    end

    st.(fldnm) = y_select;

end

