function [output,fuzzifiedIn,ruleOut,aggregatedOut,ruleFiring] = evalProjection( ...
    tsModel, input, modelRange, isWrap, sysName)
% extended evalfis: works even outside the modelRange
    [~, n] = size(modelRange);
    if isWrap
        if strcmp(sysName, 'invPend') 
            input = sys.invPendWrapper(input);
        elseif strcmp(sysName, 'flex2link')
            input = sys.flex2linkWrapper(input);
        end
    end
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
    projection(n+1:end) = input(n+1:end);
    [output,fuzzifiedIn,ruleOut,aggregatedOut,ruleFiring] = evalfis( ...
        tsModel, projection);
end
