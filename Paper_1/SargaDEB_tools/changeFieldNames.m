function newStruct = changeFieldNames(oldStruct)
    oldFields = fieldnames(oldStruct);
    newStruct = struct();
    for i = 1:length(oldFields)
        oldName = oldFields{i};
        % Extract part before the first underscore
        newName = regexp(oldName, '^[^_]+', 'match', 'once');
        newStruct.(newName) = oldStruct.(oldName);
    end
end


