function [output,fuzzifiedIn,ruleOut,aggregatedOut,ruleFiring] = evalProjection( ...
    tsModel, input, modelRange)
% extended evalfis: works even outside the modelRange
    [~, n] = size(modelRange);
    projection = zeros(size(input));
    for k=1:n
        if input(k) < modelRange(1, k)
            projection(k) = modelRange(1, k);
        elseif input(k) > modelRange(2, k)
            projection(k) = modelRange(2, k);
        else
            projection(k) = input(k);
        end
    end
    [output,fuzzifiedIn,ruleOut,aggregatedOut,ruleFiring] = evalfis( ...
        tsModel, projection);
end
