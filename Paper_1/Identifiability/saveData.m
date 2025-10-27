function saveData(paramStruct, outputDir, label)
    % Create directory and save parameters

    timeStamp = char(datetime('now','Format','dd_MMM_uuuu_HH_ss'));
    folderName = [label, '_ParVariation_', timeStamp];

    saveDir = fullfile(outputDir, folderName);

    if ~exist(saveDir, 'dir'), mkdir(saveDir); end
    save(fullfile(saveDir, 'resultsParameters.mat'), 'paramStruct');
end
