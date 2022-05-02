function extParams = addBiasNules(thenParams, nRules, n, m)
% Set bias parameters of tsModel as null.
    extParams = zeros(nRules * (n+m+1), n);
    firstLine = 1;
    for iRule=1:nRules
        extLines = firstLine : firstLine+(n+m-1);
        thenLines = 1+(iRule-1)*(n+m) : iRule*(n+m);
        extParams(extLines, :) = thenParams(thenLines, :);
        firstLine = firstLine + (n+m+1);
    end
end
