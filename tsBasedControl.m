function [u, fHat, hatB, errorFlag] = tsBasedControl( ...
    x, extendedModel, sysName, dt, known, Q, R, isWrap)
% this function finds SDRE-control based on tsModel, x - state vector
    arguments
        x
        extendedModel
        sysName
        dt double {mustBePositive}
        known = []
        Q = []
        R = []
        isWrap = false
    end

    % 1. Set default values
    errorFlag = false;
    if strcmp(sysName, 'motorLink')
        n = 2;
        r = 1;
    elseif strcmp(sysName, 'invPend') 
        n = 4;
        r = 1;
    elseif strcmp(sysName, 'flex2link')
        n = 8;
        r = 2;
    end
    tsModel = extendedModel.model;
    thenParams = extendedModel.thenParams;    
    modelRange = extendedModel.range;
    normC = extendedModel.normC;
    normS = extendedModel.normS;
    if ~isempty(normC) && ~isempty(normS)
        x = normalize(x', 2, ...
            'center', normC(1:n), 'scale', normS(1:n));
        x = x';
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
    [~, ~, ~, ~, ruleFiring, pure, processed] = utils.evalProjection( ...
        tsModel, [x; zeros(r, 1)], modelRange, isWrap, sysName);
    x_pure = pure(1:n);  % before wrapper, after projection
    x_processed = processed(1:n);  % after wrapper and projection
    
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
    if isWrap
        usual_f_hat = 1/dt * (waveA - eye(n)) * x_processed;
        hatA = diag(usual_f_hat ./ x_pure);
    else
        hatA = 1/dt * (waveA - eye(n));
    end
    hatB = 1/dt * waveB;
    if ~isempty(known)
        A = sys.get_A(x_pure, sysName);
        B = sys.get_B(x_pure, sysName);  
        hatA(known, :) = A(known, :);
        hatB(known, :) = B(known, :);
    end
    fHat = hatA * x;  % TODO: x_pure instead x?

    % 4. Calculate u = SDRE(hatA, hatB)
    P = icare(hatA, hatB, Q, R);
    if isempty(P)
        disp('Solution of SDRE has not been found')
        u = zeros(r, 1);
        errorFlag = true;
    else
        u = -inv(R) * hatB' * P * x;  % TODO: x_pure instead x?
    end
end
