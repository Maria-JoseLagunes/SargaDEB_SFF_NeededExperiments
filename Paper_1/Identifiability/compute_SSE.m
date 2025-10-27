function SSE = compute_SSE(originalPar, varPar, freeParams)
    SSE = struct();
    for i = 1:length(freeParams)
        pName = freeParams{i};
        SSE.(pName) = get_SSE(originalPar.(pName), varPar.(pName));
        SSE.(pName)(isnan(SSE.(pName))) = 0; % Handle NaN values
    end
end
