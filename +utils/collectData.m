function dataset = collectData(sysName, x0Grid, xRange, uRange, T, ...
    dt, r, pidCoefs, isWrap)
% create simulated data for TS model identification
    
    if strcmp(sysName, 'motorLink')
        rhs = @sys.rhsMotorLink;
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
    if ~isWrap
        wrapper = @(x) x;
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
                X(iStep, :) = wrapper(X(iStep, :)')';
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
                X(iStep, :) = wrapper(X(iStep, :)')';
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
