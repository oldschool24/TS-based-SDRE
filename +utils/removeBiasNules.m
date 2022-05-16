function res = removeBiasNules(extParams, nRules, n, r)
% Remove bias parameters of tsModel
    linesToRemove = (n+r+1) : (n+r+1) : nRules*(n+r+1);
    extParams(linesToRemove, :) = [];
    res = extParams;
end