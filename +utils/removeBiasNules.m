function res = removeBiasNules(extParams, nRules, n, m)
% Remove bias parameters of tsModel
    linesToRemove = (n+m+1) : (n+m+1) : nRules*(n+m+1);
    extParams(linesToRemove, :) = [];
    res = extParams;
end