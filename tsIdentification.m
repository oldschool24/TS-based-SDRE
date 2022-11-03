function tsIdentification(isLoad, sysName)
% identification of TS fuzzy model based on IO data
% x(k+1) ~ f(x(k), u(k)); 
% isLoad: 1, if load dataset; 0 if create

    if nargin == 1
        sysName = 'motorLink';
    end
    % u = [u_1 ... u_r]', x = [x_1 ... x_n]'
    if strcmp(sysName, 'motorLink')
        r = 1;      
        n = 2;
        uRange = [-0.5;
                   0.5];
        xRange = [-pi, -6*pi;
                   pi,  6*pi];
        x0 = [0; 0];
    elseif strcmp(sysName, 'invPend')
        r = 1;
        n = 4;
        uRange = [-3;
                   3];
%         uRange = [-97 138];
        xRange = [-6, -pi/2, -12, -10; ...
                   6,  pi/2,  12,  10];
        x0 = [0.1; pi/30; -0.05; -pi/10];
    elseif strcmp(sysName, 'flex2link')
        r = 2;
        n = 8;
%         uRange = [-211 -615; 
%                    254  608];
        uRange = [-1 -3.7;  % 100 times less
                   1  3.7];
        xRange = [-pi, -pi, -30, -70, -5*pi, -5*pi, -150, -700;  % real
                   pi,  pi,  30,  70,  5*pi,  5*pi,  150,  700];
        x0 = zeros(n, 1);
    end
    
    T = 1;
    dt = 0.01;
    method = 'SubtractiveClustering';
    dataName = ['data/' sysName];
    if isLoad == 1
        load(dataName, 'trainData')
    else
        if strcmp(sysName, 'motorLink')
            x0Range = [-pi/2 -2*pi; ...
                        pi/2  2*pi];
            nPoints = 4;
        elseif strcmp(sysName, 'invPend')
            x0Range = zeros(size(xRange));
            xAmp = xRange(2, :) - xRange(1, :);
            x0Range(1, :) = xRange(1, :) + 0.025 * xAmp;    % TODO: retry it for motorLink
            x0Range(2, :) = xRange(2, :) - 0.025 * xAmp;
            nPoints = 200;
        elseif strcmp(sysName, 'flex2link')
            x0Range = [-pi/2, -pi/2, -pi/2, -pi/2, -pi, -pi, -pi, -pi;
                        pi/2,  pi/2,  pi/2,  pi/2,  pi,  pi,  pi,  pi];
            nPoints = 300;
        end
        x0Grid = utils.uniformGrid(x0Range, nPoints);
        trainData = collectData(sysName, x0Grid, xRange, uRange', T, dt, r);
        save(dataName, 'trainData')
    end

    % 1. identify number of rules and antecedents params: x(k+1) ~ f(x(k))
    opt = genfisOptions(method); 
    tsModel = genfis(trainData(:, 1:n), trainData(:, end-n+1:end), opt);
    for iState=1:n
        tsModel.Inputs(iState).Name = ['x_' num2str(iState) '(t)'];
        tsModel.Outputs(iState).Name = ['x_' num2str(iState) '(t+1)'];
    end
    nRules = length(tsModel.Rules);

    % 2. identify consequents of extended model: x(k+1) ~ f(x(k), u(k))
    % Note: u(k) is not included in if-condition
    for iControl=1:r
        tsModel = addInput(tsModel, uRange(:, iControl)', ...
                           'Name', ['u_' num2str(iControl) '(t)']);
    end
    % consequent: coefficients for u are zero -> u doesn't affect output  
%     plotIdentified(tsModel, rhs, method, 20, dt)

    thenParams = bls(tsModel, trainData, r, n);
    thenParams = utils.addBiasNules(thenParams, nRules, n, r);
    thenParams = reshape(thenParams, 1, []);
    [~, out] = getTunableSettings(tsModel);
    tsModel = setTunableValues(tsModel, out, thenParams);

    % 3. Calculate the RMSE
    [nData, ~] = size(trainData); 
    RMSE = 0;
    for iData=1:nData
        pred = evalfis(tsModel, trainData(iData, 1:n+r));
        RMSE = RMSE + norm(trainData(iData, n+r+1:end) - pred) ^ 2;
%         trainData(iData, n+r+1:end) - pred
    end
    RMSE = sqrt(RMSE / nData);
    disp(['RMSE = ', num2str(RMSE)])

    tic
    % 4. plot and compare
    plotIdentified(sysName, tsModel, 10, dt, x0)
    toc

    % 5. save
    modelName = ['models/' sysName];
    writeFIS(tsModel, modelName)
end

function dataset = collectData(sysName, x0Grid, xRange, uRange, T, dt, r)
% create simulated data for TS model identification
    if strcmp(sysName, 'motorLink')
        rhs = @sys.rhsMotorLink;
        wrapper = @(x) x;
        expAmp = -2;
        trainParams = reshape(uRange, 1, [], 2);
    elseif strcmp(sysName, 'invPend')
        rhs = @sys.rhsInvPend;
        wrapper = @sys.invPendWrapper;
        expAmp = -2;
        trainParams = reshape(uRange, 1, [], 2);
    elseif strcmp(sysName, 'flex2link')
        rhs = @sys.rhsFlex2link;
        wrapper = @sys.flex2linkWrapper;
%         expAmp = [-2;
%                   -2];
        expAmp = [];
        trainParams = reshape(uRange, 2, [], 2);
    end

    [nPoints, n] = size(x0Grid);
    timesteps = 0:dt:T;
    nSteps = length(timesteps);
    uTrain = trainFuncs(trainParams, expAmp);    
    nControls = length(uTrain);
    dataset = zeros(nPoints * nControls * (nSteps-1), r + 2*n);
    iData = 0;
    if isempty(xRange)
        options = odeset('Events', []);
    else
        options = odeset('Events', ...
            @(t, x) utils.xRangeEvent(wrapper(x), xRange));
    end
    trajsLen = zeros(nPoints*nControls, 1);
    iTraj = 1;
    for iPoint=1:nPoints            % TODO: change order of for loops
        x0 = x0Grid(iPoint, :)';
        for iControl=1:nControls
            % 1. Integrate with initial condition = x0
            uList = arrayfun(uTrain{iControl}, timesteps, ...
                'UniformOutput', false);
            uList = cell2mat(uList);
            pp = spline(timesteps, uList);  
            u = @(t) ppval(pp, t);
            [t, X] = ode45(@(t, x) rhs(x, u(t)), timesteps, x0, options);
            nSteps = length(t);
            trajsLen(iTraj) = nSteps;
            iTraj = iTraj + 1;
            for iStep=1:nSteps
                X(iStep, :) = wrapper(X(iStep, :));
            end

%             % 2. Plot
%             figure()
%             plot(t, X)
%             legend('$x$','$\theta$', '$\dot{x}$', '$\dot{\theta}$', ...
%                    'Interpreter', 'latex')

            % 3. Save data from simulation
            for iStep=1:nSteps-1
                dataset(iData+iStep, 1:n) = X(iStep, :);            
                dataset(iData+iStep, n+1:n+r) = uList(:, iStep)';
                dataset(iData+iStep, r+n+1:end) = X(iStep+1, :);
            end
            iData = iData + (nSteps - 1);
        end
    end
    dataset(iData + 1:end, :) = [];
end

function res = trainFuncs(uniformInterval, expAlpha)
    amp = uniformInterval(:, :, 2) - uniformInterval(:, :, 1);
    amp = squeeze(amp);
    offset = squeeze(uniformInterval(:, :, 1));
    
    [r, nFuncs, ~] = size(uniformInterval);
    [~, nExp] = size(expAlpha);
    res = cell(nFuncs + nExp, 1);
    for iFunc=1:nFuncs
        res{iFunc} = @(t) amp(:, iFunc).*rand(r, 1) + offset(:, iFunc);
    end
    for iExp=1:nExp
        amp = uniformInterval(:, randi(nFuncs), 2);
        res{nFuncs + iExp} = @(t) amp .* exp(expAlpha(:, iExp) * t);
    end
end

function thenParams = bls(tsModel, dataset, r, n)
% this function find consequents parameters of tsModel
% using bls-algorithm on the dataset
% r = length(u), n = length(x), x - state vector, u - control vector
    [nSamples, ~] = size(dataset);
    nRules = length(tsModel.Rules);
    % 1. Extract ground truth from dataset.
    X = dataset(:, end-n+1:end);
    input = dataset(:, 1:end-n);
    % 2. Calculating the firing strength of each rule.
    firings = zeros(nSamples, nRules);
    for iSample=1:nSamples
        [~, ~, ~, ~, ruleFiring] = evalfis(tsModel, input(iSample, :));
        ruleFiring = ruleFiring / sum(ruleFiring);
        firings(iSample, :) = ruleFiring;
    end
    % 3. Create Phi = [phi(1)'; ...; phi(n_d)']
%     % comment 115 line and uncomment 113, 114 if you need bias parameter:
%     firings = repelem(firings, 1, n+r+1);
%     input(:, end+1) = 1;  % fake input for bias parameter
    firings = repelem(firings, 1, n+r);
    Phi = repmat(input, 1, nRules) .* firings;
    % 4. bls: use mldivide
    thenParams = Phi \ X;
end

function plotIdentified(sysName, tsModel, T, dt, x0)
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

    timesteps = 0:dt:T;
    nSteps = length(timesteps);
    n = length(x0);
    uTest = testFunctions(sysName);
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
        X_true(iTest, :, :) = X;
    end
    [~, ~, n] = size(X_true);
    
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
            X_pred(iTest, iStep, :) = evalfis(tsModel, ...
                                              [X(iStep-1, :)'; ...
                                              u(timesteps(iStep-1))]);
        end
        [~, f, fHat, B, Bhat] = utils.logger(sysName, X, r, tsModel, dt);
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
