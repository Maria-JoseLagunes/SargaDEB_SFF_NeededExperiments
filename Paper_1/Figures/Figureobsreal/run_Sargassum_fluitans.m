close all; 
global pets 

pets = {'Sargassum_fluitans'}; 

estim_options('default'); 
estim_options('max_step_number', 5e2); 
estim_options('max_fun_evals', 5e3); 


estim_options('pars_init_method',2); 
estim_options('results_output',3)
estim_options('method', 'no'); 





% Calculate the time it takes running the calibration
tStart = tic; 
estim_pars_algae; 
tEnd = toc(tStart);
fprintf('THIS CALIBRATION TOOK: %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));


%% Save all the results to working folder
figures = dir('*.png');
html    = dir('*.html');
results = dir('*.mat');
mfiles = dir('*.m');

results_names = {results.name}; 
mfiles_names = {mfiles.name}; 

% Use regexp to filter names containing "results"
filtered_results = results_names(~cellfun('isempty', regexp(results_names, 'results', 'once')));
filtered_results = filtered_results{1};

filtered_pars_init = mfiles_names(~cellfun('isempty', regexp(mfiles_names, 'pars_init', 'once')));
filtered_pars_init= filtered_pars_init{1};

filtered_mydata = mfiles_names(~cellfun('isempty', regexp(mfiles_names, 'mydata', 'once')));
filtered_mydata = filtered_mydata{1};

filtered_predict = mfiles_names(~cellfun('isempty', regexp(mfiles_names, 'predict', 'once')));
filtered_predict= filtered_predict{1};


filtered_addchem = mfiles_names(~cellfun('isempty', regexp(mfiles_names, 'add_chem', 'once')));
filtered_addchem = filtered_addchem{1};



results_file = {figures.name,...non
    html.name,...
   filtered_results,...
   %filtered_predict, ...
  % filtered_mydata,...
   %filtered_pars_init ...
   %filtered_addchem ...
    };   

timeStamp = char(datetime('now','Format','dd-MMM-uuuu'));
timeStamp = replace(timeStamp, ":", "_");
saveDir   = ['aumentedphotoweight_estimation_procedure_saving', timeStamp, '/'];
subDir    = [saveDir];
mkdir(subDir);

    for i = 1:length(results_file)
        copyfile(results_file{i},[subDir, results_file{i}])
    end





%%
 %    'pars_init_method':
  %      0: get initial estimates from automatized computation 
  %      1: read initial estimates from .mat file (for continuation)
  %      2: read initial estimates from pars_init file (default)
  %
  %    'results_output':
  %      0     - only saves data to .mat (no printing to html or screen and no figures) - use this for (automatic) continuations 
  %      1, -1 - no saving to .mat file, prints results to html (1) or screen (-1), shows figures but does not save them
  %      2, -2 - saves to .mat file, prints results to html (2) or screen (-2), shows figures but does not save them
  %      3, -3 - like 2 (or -2), but also prints graphs to .png files (default is 3)
  %      4, -4 - like 3 (or -3), but also prints html with implied traits
  %      5, -5 - like 4 (or -4), but includes related species in the implied traits
  %      6     - like 5, but also prints html with population traits
  %
  % 'method':
  % 'no': do not estimate
  % 'nm': Nelder-Mead simplex method (default)
  % 'mmea': multimodal evolutionary algorithm
