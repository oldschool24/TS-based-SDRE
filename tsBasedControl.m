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
%         hatA = estimateA(dt, waveA, x_pure, 'offset', hatB, Q, R);
        
        % TODO: COMMENT!!!!!!
%         hatB = sys.get_B(x, 'flex2link');      
        
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

function hatA = estimateA(dt, waveA, x, type, hatB, Q, R)
    n = numel(x);
    hatA = 1/dt * (waveA-eye(n));
    fHat = hatA * x;
    
    global w_alpha; % [zeros(n, 1); 1]
    global x_old;
    if strcmp(type, 'perMatrix')
        % solution = weighted sum of matrices(candidates)
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
        if isempty(w_alpha)
            w_alpha = [zeros(n, 1); 1];
        end
    elseif strcmp(type, 'perElement')
        % as previous, but candidates = elements, not matrix
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
        if isempty(w_alpha)
            w_alpha = ones(n^2, 1) / n;
        end
    elseif strcmp(type, 'offset')
        % solution = offset that improves the original estimate
                
        % TODO: COMMENT!!!!!!
%         hatA = sys.get_A(x, 'flex2link');

        candidates = hatA;

        % Equality constrain for all i: sum_j(dA_ij * x_j) = 0
        Aeq = zeros(n, n^2);
        for ix=1:n
            Aeq(ix, ix:n:n^2) = x(:);
        end
        beq = zeros(n, 1);
        % -eps < dA < eps
        % NOW CODE IS WORKING IN MAIN MODE!!!!!!!!!!!!!!!!!!
        dA_abs = 5e-3;
%         dA_thr = 1e-3;
%         eps = max(dA_abs * abs(hatA), dA_thr);
        eps = dA_abs * abs(hatA);
        if isempty(w_alpha)
%             eps = 1e-2 * abs(hatA);
            w_alpha = zeros(n^2, 1);
        end
        lb = reshape(-eps, [], 1);
        ub = reshape(eps, [], 1);
    end
   
    eps_x = 1;
    if norm(x - x_old) > eps_x
        wOpt = fmincon(@(w) bestFactorizationObjective(w, candidates, type, hatB, Q, R, x), ...
                      w_alpha, [], [], Aeq, beq, lb, ub, [], ...
                      optimoptions('fmincon', 'Display', 'off', ... % off, iter
                      'OptimalityTolerance', 1e-6, 'MaxIterations', 10, ...
                      'ConstraintTolerance', 1e-6));
        w_alpha = wOpt;
        x_old = x;
    end

    if strcmp(type, 'perMatrix')
        hatA = sum(candidates .* reshape(w_alpha, 1, 1, []), 3);
    elseif strcmp(type, 'perElement')
        hatA = candidates .* reshape(w_alpha, size(candidates));
    elseif strcmp(type, 'offset')
        dA = reshape(w_alpha, size(hatA));
%         disp(all(dA >= -eps, 'all') && all(dA <= eps, 'all'))
        hatA = hatA + dA;
    end
%     disp(max(abs(fHat - hatA * x)))
end

function value = bestFactorizationObjective(w, candidates, type, ...
                                            hatB, Q, R, x)
    if strcmp(type, 'perMatrix')
        % weighted sum of matrix candidates
        A = sum(candidates .* reshape(w, 1, 1, []), 3);
        value = norm(A) * norm(inv(A));
    elseif strcmp(type, 'perElement')
        matrW = reshape(w, size(candidates));
        % weighted sum of element candidates
        A = matrW .* candidates;
        value = norm(A) * norm(inv(A));
    elseif strcmp(type, 'offset')
        A = candidates;
        dA = reshape(w, size(A));
        A = A + dA;

        value = norm(A) * norm(inv(A));
%         P = icare(A, hatB, Q, R);
%         value = x' * P * x;
    end
end
