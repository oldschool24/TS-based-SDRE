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
    hatB = 1/dt * waveB;
    if isWrap && any(x_pure ~= x_processed)
%         [hatA, P, info] = unwrappedAfromWrapped(waveA, n, x_processed, ...
%                                                 x_pure, Q, R);
        [hatA, P, info] = minimalAndSimilarA(dt, waveA, x_processed, ...
                                             x_pure, sysName, hatB, Q, R);
%         stabilizable(hatA, hatB) && detectable(hatA, sqrtm(Q))  % check
    else
        hatA = 1/dt * (waveA - eye(n));
%         hatA = weightedA(dt, waveA, x_pure, 'perElement');
        [hatA, hatB] = knownChange(sysName, known, x_pure, hatA, hatB);
        [P, ~, ~, info] = icare(hatA, hatB, Q, R);
    end
    fHat = hatA * x_pure;

    % 4. Calculate u = SDRE(hatA, hatB)
    if isempty(P)
        disp('Solution of SDRE has not been found')
%         Dx = diag(info.Sx);
%         P = Dx * info.V * inv(info.U) * Dx;
        u = zeros(r, 1);
        errorFlag = true;
    else
        u = -inv(R) * hatB' * P * x;  % TODO: x_pure instead x?
    end
end

function [hatA, hatB] = knownChange(sysName, known, x_pure, hatA, hatB)
    if ~isempty(known)
        A = sys.get_A(x_pure, sysName);
        B = sys.get_B(x_pure, sysName);  
        hatA(known, :) = A(known, :);
        hatB(known, :) = B(known, :);
    end
end

function [hatA, P, info] = minimalAndSimilarA(dt, waveA, x_processed, ...
                                              x_pure, sysName, hatB, Q, R)
    n = numel(x_processed);
    if strcmp(sysName, 'flex2link')
        pIdxs = 1:2;  % periodic components
        sIdxs = 3:8;  % standard components
    else
        pIdxs = find(x_pure ~= x_processed);  
        sIdxs = find(x_pure == x_processed);  
    end
    
    % for standard: use usual estimate of A
    usualHatA = 1/dt * (waveA - eye(n));
    hatA = zeros(n);
    hatA(:, sIdxs) = usualHatA(:, sIdxs);

% remaining columns are obtained as a solution to the optim. problem
% norm(hatA(:, pIdxs)) -> min, when:
% hatA(:, pIdxs) * x_pure(pIdxs) = usualHatA(:, pIdxs) * x_processed(pIdxs)
    x_wrapped = x_processed(pIdxs);
    x_unwrapped = x_pure(pIdxs);
    usualHatA = usualHatA(:, pIdxs);
    lambda = 2 * (usualHatA*x_wrapped) / sum(x_unwrapped.^2);
    hatA(:, pIdxs) = lambda * x_unwrapped' / 2;  % lambda(i) * unwrapped(k)
%     sum(1/dt*(waveA-eye(n))*x_processed - hatA*x_pure, 'all')  % check 

    [P, ~, ~, info] = icare(hatA, hatB, Q, R);
end

function hatA = weightedA(dt, waveA, x, type)
    n = numel(x);
    fHat = 1/dt * (waveA-eye(n)) * x;
    
    if strcmp(type, 'perMatrix')
        candidates = zeros(n, n, n+1);
        for iCandidate=1:n
            A = zeros(n, n);
            if abs(x(iCandidate)) > 1e-9
                A(:, iCandidate) = fHat ./ x(iCandidate);
            end
            candidates(:, :, iCandidate) = A;
        end
        candidates(:, :, n+1) = 1/dt * (waveA-eye(n));   

        Aeq = ones(1, n+1);  % sum(w) = 1
        beq = 1;
        lb = zeros(n+1, 1);  % 0 <= w <= 1
        ub = ones(n+1, 1);
    elseif strcmp(type, 'perElement')
        candidates = fHat ./ x';  % a_ij = f_i / x_j
        [ii, jj] = find(isinf(candidates));
        linIdxs = sub2ind([n, n], ii, jj);
        candidates(isinf(candidates)) = 10^12;  % zero leads to Inf
        
        Aeq = repmat(eye(n), 1, n);
%         Aeq(:, linIdxs) = 0;
        beq = ones(n, 1);
        lb = zeros(n^2, 1);
        ub = ones(n^2, 1);
%         ub(linIdxs) = 0;
    end
   
    global w_alpha; % [zeros(n, 1); 1]
    wOpt = fmincon(@(w) bestFactorizationObjective(w, candidates, type), ...
                  w_alpha, [], [], Aeq, beq, lb, ub, [], ...
                  optimoptions('fmincon', 'Display', 'off', ...
                  'OptimalityTolerance', 1e-6, 'MaxIterations', 10, ...
                  'ConstraintTolerance', 1e-6));
    w_alpha = wOpt;

    if strcmp(type, 'perMatrix')
        hatA = sum(candidates .* reshape(wOpt, 1, 1, []), 3);
    elseif strcmp(type, 'perElement')
        hatA = candidates .* reshape(wOpt, size(candidates));
    end
%     fHat - hatA * x
end

function value = bestFactorizationObjective(w, candidates, type)
    if strcmp(type, 'perMatrix')
        weighted = sum(candidates .* reshape(w, 1, 1, []), 3);
    elseif strcmp(type, 'perElement')
        matrW = reshape(w, size(candidates));
        weighted = matrW .* candidates;
    end
    value = norm(weighted) * norm(inv(weighted));
end
