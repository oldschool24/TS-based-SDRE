function [u, fHat, hatB, errorFlag] = tsBasedControl(x, sysName, tsModel, ...
                                                     dt, known, Q, R)
% this function finds SDRE-control based on tsModel, x - state vector
    errorFlag = false;
    if strcmp(sysName, 'motorLink')
        n = 2;
        r = 1;
    elseif strcmp(sysName, 'invPend') 
        x = sys.invPendWrapper(x);
        n = 4;
        r = 1;
    elseif strcmp(sysName, 'flex2link')
        x = sys.flex2linkWrapper(x);
        n = 8;
        r = 2;
    end
    nRules = length(tsModel.Rules);

    if nargin < 5
        known = [];
        Q = 10 * eye(n);
        R = 5 * eye(r);
    end

    % 1. Calculate waveA, waveB: 
    %    tsModel(x(k), u(k)) = waveA * x(k) + waveB * u(k)
    waveA = zeros(n, n);
    waveB = zeros(n, r);
    [~, ~, ~, ~, ruleFiring] = evalfis(tsModel, [x; zeros(r, 1)]);
    
    if sum(ruleFiring) < 1e-20
        disp(['The rules do not work. x:' num2str(x')])
        u = zeros(r, 1);
        fHat = zeros(n, 1);
        hatB = zeros(n, r);
        errorFlag = true;
        return
    end
    
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
    hatA = 1/dt * (waveA - eye(n));
    hatB = 1/dt * waveB;
    if ~isempty(known)
        A = sys.get_A(x, sysName);
        B = sys.get_B(x, sysName);
        hatA(known, :) = A(known, :);
        hatB(known, :) = B(known, :);
    end
    fHat = hatA * x;

    % 3. Calculate u = SDRE(hatA, hatB)
    P = icare(hatA, hatB, Q, R);
    if isempty(P)
        disp('Solution of SDRE has not been found')
        u = zeros(r, 1);
        errorFlag = true;
    else
        u = -inv(R) * hatB' * P * x;
    end
end

    % analysis of tsModel coefs 
%     [~, out] = getTunableSettings(tsModel);
%     theta = getTunableValues(tsModel, out);
% uWeights = theta(5:6:360);
% uWeights = reshape(uWeights, 4, []);
% sum(abs(uWeights), 2), max(abs(uWeights), [], 2), mean(abs(uWeights), 2)
