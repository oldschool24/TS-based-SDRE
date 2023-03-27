function [u, fHat, hatB, errorFlag] = tsBasedControl( ...
    x, extendedModel, sysName, dt, known, Q, R)
% this function finds SDRE-control based on tsModel, x - state vector
    arguments
        x
        extendedModel
        sysName
        dt double {mustBePositive}
        known = []
        Q = []
        R = []
    end

    % 1. Set default values
    errorFlag = false;
    if strcmp(sysName, 'motorLink')
        wrapped_x = x;
        n = 2;
        r = 1;
    elseif strcmp(sysName, 'invPend') 
        wrapped_x = sys.invPendWrapper(x);
        n = 4;
        r = 1;
    elseif strcmp(sysName, 'flex2link')
        wrapped_x = sys.flex2linkWrapper(x);
        n = 8;
        r = 2;
    end
    tsModel = extendedModel.model;
    thenParams = extendedModel.thenParams;    
    modelRange = extendedModel.range;
    normC = extendedModel.normC;
    normS = extendedModel.normS;
    if ~isempty(normC) && ~isempty(normS)
        wrapped_x = normalize(wrapped_x', 2, ...
            'center', normC(1:n), 'scale', normS(1:n));
        wrapped_x = wrapped_x';
    end
    if isempty(Q)
        Q = 10 * eye(n);
    end
    if isempty(R)
        R = 5 * eye(r);
    end
    nRules = length(tsModel.Rules);

    % 2. Calculate waveA, waveB: 
    %    tsModel(x(k), u(k)) = waveA * x(k) + waveB * u(k)
    waveA = zeros(n, n);
    waveB = zeros(n, r);
    [~, ~, ~, ~, ruleFiring] = utils.evalProjection( ...
        tsModel, [wrapped_x; zeros(r, 1)], modelRange);
    
    ruleFiring = ruleFiring / sum(ruleFiring);
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

    % 3. Calculate hatA, hatB: estimates of A(x), B(x)
    hatA = 1/dt * (waveA - eye(n));
    hatB = 1/dt * waveB;
    if ~isempty(known)
        A = sys.get_A(x, sysName);
        B = sys.get_B(x, sysName);  
        hatA(known, :) = A(known, :);
        hatB(known, :) = B(known, :);
    end
    fHat = hatA * x;

    % 4. Calculate u = SDRE(hatA, hatB)
    P = icare(hatA, hatB, Q, R);
    if isempty(P)
        disp('Solution of SDRE has not been found')
        u = zeros(r, 1);
        errorFlag = true;
    else
        u = -inv(R) * hatB' * P * x;
    end
end
