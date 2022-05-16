function extParams = addBiasNules(thenParams, nRules, n, r)
% Set bias parameters of tsModel as null.
    extParams = zeros(nRules * (n+r+1), n);
    firstLine = 1;
    for iRule=1:nRules
        extLines = firstLine : firstLine+(n+r-1);
        thenLines = 1+(iRule-1)*(n+r) : iRule*(n+r);
        extParams(extLines, :) = thenParams(thenLines, :);
        firstLine = firstLine + (n+r+1);
    end
end
