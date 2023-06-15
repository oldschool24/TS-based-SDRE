function [output,fuzzifiedIn,ruleOut,aggregatedOut,ruleFiring, pure, processed] = evalProjection( ...
    tsModel, input, modelRange, isWrap, sysName)
% extended evalfis: works even outside the modelRange
    
    if isWrap
        if strcmp(sysName, 'invPend') 
            [wrapped, periods] = sys.invPendWrapper(input);
        elseif strcmp(sysName, 'flex2link')
            [wrapped, periods] = sys.flex2linkWrapper(input);
        end
    end
    
    projection = zeros(size(input));
    [~, n] = size(modelRange);
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

    if isWrap
        if strcmp(sysName, 'invPend') 
            pure = [input(1); projection(2:end)];
            processed = [wrapped(1); projection(2:end)];
        elseif strcmp(sysName, 'flex2link')
            pure = [input(1:2); projection(3:end)];
            processed = [wrapped(1:2); projection(3:end)];
        end
    else
        pure = projection;
        processed = projection;
    end

    [output,fuzzifiedIn,ruleOut,aggregatedOut,ruleFiring] = evalfis( ...
        tsModel, processed);
    output = output + periods';
end
