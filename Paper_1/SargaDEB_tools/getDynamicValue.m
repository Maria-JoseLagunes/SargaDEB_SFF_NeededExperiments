function val = getDynamicValue(variable, index)
    if numel(variable) > 1
        val = variable(index);
    else
        val = variable;
    end
end