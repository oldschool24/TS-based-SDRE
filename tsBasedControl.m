function [u, fHat, hatB] = tsBasedControl(x, sysName, tsModel, learnStep)
% this function finds SDRE-control based on tsModel, x - state vector
    if strcmp(sysName, 'invPend') 
        x = sys.invPendWrapper(x);
    end

    nRules = length(tsModel.Rules);
    n = length(x);
    nInputs = length(tsModel.Inputs);
    r = nInputs - n;

    % 1. Calculate waveA, waveB: 
    %    tsModel(x(k), u(k)) = waveA * x(k) + waveB * u(k)
    waveA = zeros(n, n);
    waveB = zeros(n, r);
    [~, ~, ~, ~, ruleFiring] = evalfis(tsModel, [x; zeros(r, 1)]);
    
    if sum(ruleFiring) < 1e-3
        disp(['The rules do not work. x:' num2str(x')])
        u = zeros(r, 1);
        fHat = zeros(n, 1);
        hatB = zeros(n, 1);
        return
    end
%     assert(sum(ruleFiring)>1e-3, ['The rules do not work. x:' num2str(x')])
    
    ruleFiring = ruleFiring / sum(ruleFiring);
    [~, out] = getTunableSettings(tsModel);
    thenParams = getTunableValues(tsModel, out);
    % column <-> component of state vector:
    thenParams = reshape(thenParams, [], n);
    thenParams = utils.removeBiasNules(thenParams, nRules, n, r);
%     theta = zeros(nRules, n, n+r);
%     % theta = reshaped thenParams: rule <-> matrix of parameters
    for iRule=1:nRules
        temp = thenParams(1+(iRule-1)*(n+r) : iRule*(n+r), :)';
        A = temp(:, 1:n);
        B = temp(:, n+1:end);
        waveA = waveA + ruleFiring(iRule) * A;
        waveB = waveB + ruleFiring(iRule) * B;
%         theta(iRule, :, :) = thenParams(1+(iRule-1)*(n+r):iRule*(n+r), :)';
    end

    % 2. Calculate hatA, hatB: estimates of A(x), B(x)
    hatA = 1/learnStep * (waveA - eye(n));
    hatB = 1/learnStep * waveB;
    fHat = hatA * x;

    % 3. Calculate u = SDRE(hatA, hatB)
    Q = 10;
    R = 5;
    P = icare(hatA, hatB, Q, R);
    u = -inv(R) * hatB' * P * x;
end
