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
        [hatA, hatB] = knownChange(sysName, known, x_pure, hatA, hatB);
        [P, ~, ~, info] = icare(hatA, hatB, Q, R);
    end
    fHat = hatA * x_pure;

    % 4. Calculate u = SDRE(hatA, hatB)
    if isempty(P)
        disp('Solution of SDRE has not been found')
        Dx = diag(info.Sx);
        P = Dx * info.V * inv(info.U) * Dx;
%         u = zeros(r, 1);
        errorFlag = true;
    end
    u = -inv(R) * hatB' * P * x;  % TODO: x_pure instead x?
end

function choices = getChoices(maxNum, nonzeroIdxs, zeroIdxs)
    % choice = join(permutation(nonzeroIdxs), combination(nonzeroIdxs))
    nZero = numel(zeroIdxs);
    permutations = perms(nonzeroIdxs);  % without repetitions
    nPerm = size(permutations, 1);

    if nZero == 0
        nChoices = min(maxNum, nPerm);
        choices = datasample(1:nPerm, nChoices, 'Replace', false);
        choices = permutations(choices, :);
    else
        combinations = randsample(nonzeroIdxs, 10*nZero, true);
        combinations = reshape(combinations, [10, nZero]);
        nComb = size(combinations, 1);
    
        [permIdxs, combIdxs] = ndgrid(1:nPerm, 1:nComb);
        nChoices = min(maxNum, numel(permIdxs));
        choices = datasample(1:numel(permIdxs), nChoices, 'Replace', false);
        choices = [permutations(permIdxs(choices), :), ...
                   combinations(combIdxs(choices), :)];
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

function isHold = stabilizable(A, B)
    isHold = true;  % has full row rank for all RE(lamda)>=0
    [n, ~] = size(A);
    eigValues = eig(A);
    eigValues(real(eigValues) < 0) = []; % only Re(lambda)>=0
    for iValue=1:length(eigValues)
        lambda = eigValues(iValue);
        if rank([A - lambda*eye(n), B]) ~= n
            isHold = false;
            break
        end
    end
end

function isHold = detectable(A, C)
    isHold = true;  % has full column rank for all RE(lamda)>=0
    [~, n] = size(A);
    eigValues = eig(A);
    eigValues(real(eigValues) < 0) = []; % only Re(lambda)>=0
    for iValue=1:length(eigValues)
        lambda = eigValues(iValue);
        if rank([A - lambda*eye(n); C]) ~= n
            isHold = false;
            break
        end
    end
end

function outputArray = reorderRows(choices, f, x, ...
                                   nNZ, almostZeroIdxs, zeroIdxs)
    arrNorm = zeros(size(choices, 1), 1);
    for iRow = 1:size(choices, 1)
        % Calculate norm of matrix A
        hatA = choice2hatA(choices(iRow, :), f, x, ...
                           nNZ, almostZeroIdxs, zeroIdxs);
        arrNorm(iRow) = norm(hatA);
    end
    
    % Sort the rows based on the number of unique elements
    % if equal, take the row corresponding to the smaller matrix norm 
    [~, sortOrder] = sortrows(arrNorm, 'ascend'); 
    
    % Reorder the rows based on the sorting order
    outputArray = choices(sortOrder, :);  
end

function hatA = choice2hatA(choice, f, x, nNZ, almostZeroIdxs, zeroIdxs)
    n = numel(x);
    nAZ = numel(almostZeroIdxs); 
    hatA = zeros(n);
    for row=1:n
        if row > nNZ
            if row <= nNZ+nAZ
                % A(row,k)*x_k + A(row,s)*x_s = f(row)
                % A(row,k)=-1 => A(row,s) = (f(row)+x_k) / x_s
                column = almostZeroIdxs(row-nNZ);
                f(row) = f(row) + x(column);
            else
                % fill in the elements that:
                % 1. affect the stabilizability and detectability
                % 2. do not affect hatA * x_pure
                column = zeroIdxs(row-nNZ-nAZ);
            end
            hatA(row, column) = -1;
        end
        column = choice(row);
        hatA(row, column) = f(row) / x(column);
    end
end

function [hatA, P, info] = unwrappedAfromWrapped(dt, waveA, x_processed, ...
                                                 x_pure, Q, R)
    n = numel(x_processed);
    usual_f_hat = 1/dt * (waveA - eye(n)) * x_processed;
    % must find such hatA: usual_f_hat = hatA * pure_x
    % n equations, n^2 elements in A => choose n elements of A
    % a_ik is candidate, if pure_x(k) ~= 0
    nonzeroIdxs = find(abs(x_pure) > 1e-6);
    almostZeroIdxs = find(abs(x_pure) < 1e-6 & x_pure);
    zeroIdxs = setdiff(setdiff(1:n, nonzeroIdxs), almostZeroIdxs);
    nNZ = numel(nonzeroIdxs);
    maxNumChoices = 100;
    choices = getChoices(maxNumChoices, nonzeroIdxs, ...
                         union(almostZeroIdxs, zeroIdxs));
    choices = reorderRows(choices, usual_f_hat, x_pure, ...
                          nNZ, almostZeroIdxs, zeroIdxs); 

%         nonzeroND = cell(nNZ, 1);
%         nonzeroND(:) = {nonzeroIdxs};  % every cell contains nonzeroIdxs
%         combs = cell(nNZ, 1);
%         [combs{:}] = ndgrid(nonzeroND{:}); % combs(:, k) = valid choice
    
%         nChoices = min(300, numel(combs{1}));
%         combsArr = zeros(nChoices, nNZ);
%         iArr = 1;
%         for iChoice=datasample(1:numel(combs{1}), nChoices, 'Replace', false)
%             for k=1:nNZ
%                 combsArr(iArr, k) = combs{k}(iChoice);
%             end
%             iArr = iArr + 1;
%         end

    sqrtmQ = sqrtm(Q);
    nChoices = size(choices, 1);
    for iChoice=1:nChoices  % TODO: continuity?
        hatA = choice2hatA(choices(iChoice, :), usual_f_hat, x_pure,...
                           nNZ, almostZeroIdxs, zeroIdxs);
        [hatA, hatB] = knownChange(sysName, known, x_pure, hatA, hatB);
        if stabilizable(hatA, hatB) && detectable(hatA, sqrtmQ)
            % then unique positive semi-definite solution exists
            % sum(usual_f_hat - hatA * x_pure, 'all')
            [P, ~, ~, info] = icare(hatA, hatB, Q, R);
            break
        end
    end
    if iChoice == nChoices
        [P, ~, ~, info] = icare(hatA, hatB, Q, R);
    end

%         selected = randsample(nonzeroIdxs, n, true);
%         hatA = zeros(n);
%         for k=1:n
%             column = selected(k);
%             row = k;
%             hatA(row, column) = usual_f_hat(row) / x_pure(column);
%         end

%         hatA = diag(usual_f_hat ./ x_pure);
%         infIdxs = isinf(hatA);
%         if any(infIdxs, 'all')
%             hatA(infIdxs) = 0;
%             options = optimset('MaxFunEvals', 200*(n^2), ...
%                                'MaxIter', 200*(n^2), ...
%                                'TolFun', 1e-3, 'TolX', 1e-3);
%             [hatA, ~, ~, optInfo] = fminsearch( ...
%                 @(A) norm(usual_f_hat - A*x_pure), hatA, options);
%             disp(optInfo)
%         end
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
