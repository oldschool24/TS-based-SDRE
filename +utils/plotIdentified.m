function plotIdentified(sysName, extendedModel, T, dt, x0)
    if strcmp(sysName, 'motorLink')
        rhs = @sys.rhsMotorLink;
        wrapper = @(x) x;
        r = 1;
    elseif strcmp(sysName, 'invPend')
        rhs = @sys.rhsInvPend;
        wrapper = @sys.invPendWrapper;
        r = 1;
    elseif strcmp(sysName, 'flex2link')
        rhs = @sys.rhsFlex2link;
        wrapper = @sys.flex2linkWrapper;
        r = 2;
    end

    n = length(x0);
    if isempty(extendedModel.normC) || isempty(extendedModel.normS)
        isNormalize = false;
    else
        isNormalize = true;
        xNormC = extendedModel.normC(1:n);
        uNormC = extendedModel.normC(n+1:n+r);
        xNormS = extendedModel.normS(1:n);
        uNormS = extendedModel.normS(n+1:n+r);
    end

    timesteps = 0:dt:T;
    nSteps = length(timesteps);    
    uTest = valFunctions(sysName);
    nTests = length(uTest);
    X_true = zeros(nTests, nSteps, n);
    for iTest=1:nTests
        % use spline approximation of random control
        uList = arrayfun(uTest{iTest}, timesteps, 'UniformOutput', false);
        uList = cell2mat(uList);
        pp = spline(timesteps, uList);  
        u = @(t) ppval(pp, t);  
        % collect true answers
        [~, X] = ode45(@(t, x) rhs(x, u(t)), timesteps, x0);
        for iStep=1:nSteps
            X(iStep, :) = wrapper(X(iStep, :));
        end
        if isNormalize
            X = normalize(X, 'center', xNormC, 'scale', xNormS);
        end
        X_true(iTest, :, :) = X;
    end
    [~, ~, n] = size(X_true);
    
    tsModel = extendedModel.model;
    modelRange = extendedModel.range;
    X_pred = zeros(nTests, nSteps-1, n);     % X_pred(1:2) = x0(1);
    X = zeros(nSteps, n);
    fTrue = zeros(nTests, nSteps, n);
    fPred = zeros(nTests, nSteps, n);
    Btrue = zeros(nTests, nSteps, n, r);
    Bpred = zeros(nTests, nSteps, n, r);
%     warning('off', 'fuzzy:general:warnEvalfis_NoRuleFired')
%     warning('off', 'fuzzy:general:diagEvalfis_OutOfRangeInput')
    for iTest=1:nTests
        u = uTest{iTest};
        X(:, :) = X_true(iTest, :, :);
        for iStep=2:nSteps
            if isNormalize
                uNormalized = normalize(u(timesteps(iStep-1))', ...
                    'center', uNormC, 'scale', uNormS);
                inp = [X(iStep-1, :), uNormalized]';
            else
                inp = [X(iStep-1, :)'; u(timesteps(iStep-1))];
            end
            X_pred(iTest, iStep, :) = utils.evalProjection(tsModel, inp, ...
                                                           modelRange);
        end
        [~, f, fHat, B, Bhat] = utils.logger(sysName, X, r, extendedModel, dt);
        fTrue(iTest, :, :) = f;
        fPred(iTest, :, :) = fHat;
        Btrue(iTest, :, :, :) = B;
        Bpred(iTest, :, :, :) = Bhat;
    end
%     warning('on', 'fuzzy:general:warnEvalfis_NoRuleFired')
%     warning('on', 'fuzzy:general:diagEvalfis_OutOfRangeInput')
    
    for iTest=1:nTests
        utils.plotEstimates('X', squeeze(X_true(iTest, :, :)), ...
                            squeeze(X_pred(iTest, :, :)), n, timesteps)
        utils.plotEstimates('f', squeeze(fTrue(iTest, :, :)), ...
                            squeeze(fPred(iTest, :, :)), n, timesteps)
        utils.plotEstimates('B', squeeze(Btrue(iTest, :, :, :)), ...
                            squeeze(Bpred(iTest, :, :, :)), n, timesteps)
    end
end
