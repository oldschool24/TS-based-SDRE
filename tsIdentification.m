function tsIdentification(isLoad, sysName, dt)
% identification of TS fuzzy model based on IO data
% x(k+1) ~ f(x(k), u(k)); 
% isLoad: 1, if load dataset; 0 if create
    arguments
        isLoad {mustBeNumericOrLogical}
        sysName
        dt double {mustBePositive}
    end
    
    utils.setDefaultVars;
    modelName = [sysName '(not-dt-' num2str(dt) ...
                 '_T-' num2str(T) ...
                 '_N-' num2str(nPoints) ...
                 '_reduct-' num2str(reduction) ...
                 ').mat'];

    dataName = ['data/train/' modelName];
    if isLoad
        load(dataName, 'trainData')
    else
        x0Grid = utils.uniformGrid(x0Range, nPoints);
        if withPI   
%             % trajectory -> ss-model -> PI-control -> dataset 
%             trainData = collectData(sysName, x0Grid(1, :), xRange, ...
%                 uRange'/10, T, dt, r, []);
%             % automatic search of Kp, Ki (not reliable)
%             [Kp, Ki] = utils.buildPI(trainData, xIdxPI, uIdxPI, dt, n);

            pidCoefs = struct('Kp', Kp, 'Ki', Ki, 'Kd', Kd, ...
                'xIdxPI', xIdxPI, 'uIdxPI', uIdxPI);
            trainData = collectData(sysName, x0Grid(2:end, :), xRange, ...
                uRange, T, dt, r, pidCoefs);
        else
            trainData = collectData(sysName, x0Grid, xRange, ...
                uRange', T, dt, r, []);
        end
        if isNormalize
            [trainData, normC, normS] = normalize( ...
                trainData, 'range', [-1, 1]);
        end
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

    thenParams = bls(tsModel, trainData, r, n);
    thenParams = utils.addBiasNules(thenParams, nRules, n, r);
    thenParams = reshape(thenParams, 1, []);
    [~, out] = getTunableSettings(tsModel);
    tsModel = setTunableValues(tsModel, out, thenParams);
    extendedModel.thenParams = thenParams;

    % 3. save
    modelName = ['models/' modelName];
    extendedModel.model = tsModel;
    extendedModel.range = utils.getTsRange(tsModel, n);
    if isNormalize
        extendedModel.normC = normC(1:n+r);
        extendedModel.normS = normS(1:n+r);
    else
        extendedModel.normC = [];
        extendedModel.normS = [];
    end
    save(modelName, "extendedModel")
    % writeFIS(extendedModel, modelName)

    % 4. plot and compare
    tic
    utils.plotIdentified(sysName, extendedModel, 5, dt, x0)
    toc
end

function dataset = collectData(sysName, x0Grid, xRange, uRange, T, ...
    dt, r, pidCoefs)
% create simulated data for TS model identification
    if strcmp(sysName, 'motorLink')
        rhs = @sys.rhsMotorLink;
        wrapper = @(x) x;
        expAlpha = -2;
        trainParams = reshape(uRange', 1, [], 2);
    elseif strcmp(sysName, 'invPend')
        rhs = @sys.rhsInvPend;
        wrapper = @sys.invPendWrapper;
        expAlpha = -2;
        trainParams = reshape(uRange', 1, [], 2);
    elseif strcmp(sysName, 'flex2link')
        rhs = @sys.rhsFlex2link;
        wrapper = @sys.flex2linkWrapper;
%         expAlpha = [-2;
%                     -2];
        expAlpha = [];
        trainParams = reshape(uRange', 2, [], 2);
    end

    if ~isempty(pidCoefs)
        xIdxPI = pidCoefs.xIdxPI;
        uIdxPI = pidCoefs.uIdxPI;
        Kp = pidCoefs.Kp;
        Ki = pidCoefs.Ki;
        Kd = pidCoefs.Kd;
    end

    uTrain = trainFuncs(trainParams, expAlpha);  
    nControls = length(uTrain);
    [nPoints, n] = size(x0Grid);
    timesteps = 0:dt:T;
    nSteps = length(timesteps);
    dataset = zeros(nPoints * (nControls+1) * (nSteps-1), r + 2*n);
    iData = 0;
    if isempty(xRange)
        options = odeset('Events', []);
    else
        options = odeset('Events', ...
            @(t, x) utils.xRangeEvent(wrapper(x), xRange));
        if ~isempty(pidCoefs)
            integralRange = [-Inf(1, length(xIdxPI));
                              Inf(1, length(xIdxPI))];
            xRange = [xRange integralRange];
        end
        options_extended = odeset('Events', ...
            @(t, x) utils.xRangeEvent(wrapper(x), xRange));
    end
    trajsLen = zeros(nPoints, nControls+1);
    outOfRangePI = [];
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
            trajsLen(iPoint, iControl) = nSteps;
            for iStep=1:nSteps
                X(iStep, :) = wrapper(X(iStep, :));
            end

            % 2. Save data from simulation
            for iStep=1:nSteps-1
                dataset(iData+iStep, 1:n) = X(iStep, :);            
                dataset(iData+iStep, n+1:n+r) = uList(:, iStep)';
                dataset(iData+iStep, n+r+1:end) = X(iStep+1, :);
            end
            iData = iData + (nSteps - 1);
        end
        if ~isempty(pidCoefs)
            % 1. Integrate with initial condition = x0
            [t, X, ~, ~, ie] = ode45( ...
                @(t, x) utils.rhsWithPI(x, sysName, xIdxPI, uIdxPI, Kp, Ki, Kd), ...
                timesteps, [x0; zeros(r, 1)], options_extended);
            outOfRangePI = [outOfRangePI; ie];

            % 2. Data postprocessing
            nSteps = length(t);
            trajsLen(iPoint, end) = nSteps;
            uList = zeros(nSteps, r);
            for iStep=1:nSteps
                for k=1:r
                    err = X(iStep, xIdxPI(k));   % reference = 0
                    I = X(iStep, n+k);
                    uList(iStep, uIdxPI(k)) = Kp(k)*err + Ki(k)*I;
                end
                X(iStep, :) = wrapper(X(iStep, :));
            end
            X(:, end-r+1 : end) = [];

            % 3. Save data from simulation
            for iStep=1:nSteps-1
                dataset(iData+iStep, 1:n) = X(iStep, :);  
                dataset(iData+iStep, n+1:n+r) = uList(iStep, :);
                dataset(iData+iStep, n+r+1:end) = X(iStep+1, :);
            end
            iData = iData + (nSteps - 1);
        end
    end
    dataset(iData + 1:end, :) = [];
    
    if ~isempty(outOfRangePI)
        [numberOfTimes, violatedComponent] = groupcounts(outOfRangePI);
        tableOfViolations = table(violatedComponent, numberOfTimes);
        disp('going beyond xRange under PI-control')
        disp(tableOfViolations)
    end

    if isempty(expAlpha)
        disp(['mean trajectory length: random control -- ' ...
            num2str(mean(trajsLen(:, 1))) ...       
            ', PI control -- ' ...
            num2str(mean(trajsLen(:, 2)))])     % TODO: case nFuncs > 1
    else
        disp(['mean trajectory length: random control -- ' ...
              num2str(mean(trajsLen(:, 1))) ...
              ', sin*exp control -- ' ...
              num2str(mean(trajsLen(:, 2))) ...  % TODO: case nFuncs > 1
              ', PI control -- ' ...
              num2str(mean(trajsLen(:, 3)))])
    end
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
        alpha = expAlpha(:, iExp);
        res{nFuncs + iExp} = @(t) amp .* sin(5* alpha * t) .* exp(alpha * t);
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
