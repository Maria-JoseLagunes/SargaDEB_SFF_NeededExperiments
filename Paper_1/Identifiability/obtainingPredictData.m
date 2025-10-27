function  [prdData, dataVec, data2plot, prdData_x] = obtainingPredictData(par, metaPar, txtPar, data, auxData, metaData, txtData, weights)
global results_output pets n_pets

if ~results_output == 0 % plot figures
    txt0 = '0'; % text to prepend to figure counter if counter < 10 (for intuitive listing sequence in the directory)
    data2plot = data; % with pets as first field names
    close all % to avoid saving figures generated prior the current run
    for i = 1:n_pets % scan pets
      st = data2plot.(pets{i}); if isfield(st,'psd'); st = rmfield(st,'psd'); end
      [nm, nst] = fieldnmnst_st(st);
      for j = 1:nst  % replace univariate data by plot data 
        fieldsInCells = textscan(nm{j},'%s','Delimiter','.');
        varData = getfield(st, fieldsInCells{1}{:});   % scalar, vector or matrix with data in field nm{j}
        k = size(varData, 2);
        if k == 1
          st = rmfield(st,nm{j}); % remove zero-variate data from structure
        else % k > 1: uni- or bi-variate data set
          auxDataFields = fields(auxData.(pets{i}));
          dataCode = fieldsInCells{1}{:};
          univarAuxData = {};
          for ii = 1:length(auxDataFields) % add to univarAuxData all auxData for the data set that has length > 1
            if isfield(auxData.(pets{i}), 'treat') && isfield(auxData.(pets{i}).treat, dataCode) && auxData.(pets{i}).treat.(dataCode){1} == 0
              univarAuxData{end + 1} = auxDataFields{ii};
            end
          end
          aux = getfield(st, fieldsInCells{1}{:});
          dataVec = aux(:,1); 
          if isfield(auxData.(pets{i}), 'treat') && isfield(auxData.(pets{i}).treat, dataCode) && auxData.(pets{i}).treat.(dataCode){1} == 0
            % if auxData.(pets{i}).treat is specified for this dataCode and equals 1, interpolation is allowed
            xAxis = dataVec;
            
          else
            n_x = 100; xAxis = linspace(min(dataVec), max(dataVec), n_x)';
          end
          st = setfield(st, fieldsInCells{1}{:}, xAxis);
        end
      end % end of scan data sets of pets{i}
      data2plot.(pets{i}) = st;
      
    end % end of scan pets
    
    %
    
    prdData = predict_pets(par, data, auxData); % data prediction
    prdData_x = predict_pets(par, data2plot, auxData); % data prediction with xAxis
    counter_fig = 0;
    
end
