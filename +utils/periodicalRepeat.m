function data = periodicalRepeat(data, nPeriod, sysName)
    if strcmp('flex2link', sysName)
        period = 2*pi;
        pIdxs = [1, 2];
        n = 8; 
        r = 2;
    end
    nP = length(pIdxs);  % number of periodic components

    [nSamples, ~] = size(data);
    shift = period * randi([-nPeriod nPeriod], nSamples, nP);
    data(:, pIdxs) = data(:, pIdxs) + shift;
    data(:, pIdxs+n+r) = data(:, pIdxs+n+r) + shift;     
end
