function tsIdentification(isLoad, sysName)
% identification of TS fuzzy model based on IO data
% x(k+1) ~ f(x(k), u(k)); 
% isLoad: 1, if load dataset; 0 if create

    if nargin == 1
        sysName = 'motor_link';
    end
    % u = [u_1 ... u_m]', x = [x_1 ... x_n]'
    if strcmp(sysName, 'motor_link')
        rhs = @sys.rhsMotorLink;
        m = 1;      
        n = 2;
    else
        rhs = @sys.rhsInvPend;
        m = 1;
        n = 4;
    end
    
    T = 1;
    learnStep = 0.01;
    method = 'SubtractiveClustering';
    dataName = ['data/' sysName];
    if isLoad == 1
        load(dataName, 'dataset')
    else
        if strcmp(sysName, 'motor_link')
            x0Ranges = [-pi/2 pi/2; -2*pi 2*pi];
        else
            x0Ranges = [-0.1 0.1; 0 2*pi; -1 1; -pi pi];
        end
        x0Grid = uniformGrid(x0Ranges, 4);
        dataset = collectData(rhs, x0Grid, T, learnStep, m);
        save(dataName, 'dataset')
    end

    % 1. identify number of rules and antecedents params: x(k+1) ~ f(x(k))
    opt = genfisOptions(method); 
    tsModel = genfis(dataset(:, 1:n), dataset(:, end-n+1:end), opt);
    for k=1:n
        tsModel.Inputs(k).Name = ['x_' num2str(k) '(t)'];
        tsModel.Outputs(k).Name = ['x_' num2str(k) '(t+1)'];
    end
    nRules = length(tsModel.Rules);

    % 2. identify consequents of extended model: x(k+1) ~ f(x(k), u(k))
    % Note: u(k) is not included in if-condition
    tsModel = addInput(tsModel, [-0.5 0.5], 'Name', 'u(t)');
    % consequent: coefficients for u are zero -> u doesn't affect output  
%     plotIdentified(tsModel, rhs, method, 20, learnStep)

    thenParams = bls(tsModel, dataset, m, n);
    thenParams = utils.addBiasNules(thenParams, nRules, n, m);
    thenParams = reshape(thenParams, 1, []);
    [~, out] = getTunableSettings(tsModel);
    tsModel = setTunableValues(tsModel, out, thenParams);

    % 3. plot and compare
    plotIdentified(tsModel, rhs, method, 20, learnStep)

    % 4. save
    modelName = ['models/' sysName];
    writeFIS(tsModel, modelName)
end

function x0Grid = uniformGrid(ranges, nPoints)
% create grid of initial points from uniform distribution 
% ranges(i, :) -- range of i-th component
% nDivisions -- the number of segments into which the range is divided
    [n, ~] = size(ranges);  % n -- length(state vector)
    x0Grid = zeros(nPoints, n);
    for iComponent=1:n
        scale = ranges(iComponent, 2) - ranges(iComponent, 1);
        bias = ranges(iComponent, 1);
        x0Grid(:, iComponent) = scale * rand(nPoints, 1) + bias;
    end
end

function dataset = collectData(rhs, x0Grid, T, learnStep, m)
% create simulated data for TS model identification
    [nPoints, n] = size(x0Grid);
    timesteps = 0:learnStep:T;
    nSteps = length(timesteps);
    dataset = zeros(nPoints * (nSteps-1), m + 2*n);
    for iPoint=1:nPoints
        % 1. Simulate with x0 = x0Grid(iPoint, :)'
        uList = rand(nSteps, m) - 0.5;
        pp = spline(timesteps, uList); % future: problems with m > 1
        uRand = @(t) ppval(pp, t);  % random uniform control on [-0.5, 0.5]
        [~, X] = ode45(@(t, x) rhs(x, uRand(t)), timesteps, x0Grid(iPoint, :)');
        
        % 2. Plot
%         hold on
%         for k=1:n         % future: comment
%             plot(X(:, k)) 
%         end
%         legend('x', 'theta', 'x_dot', 'theta_dot')
%         hold off

        % 3. Save data from simulation
        iData = (iPoint-1) * (nSteps-1);
        for iStep=1:nSteps-1
            dataset(iData+iStep, 1:n) = X(iStep, :);            
            dataset(iData+iStep, n+1:n+m) = uList(iStep, :);
            dataset(iData+iStep, m+n+1:end) = X(iStep+1, :);
        end
    end
end

function thenParams = bls(tsModel, dataset, m, n)
% this function find consequents parameters of tsModel
% using bls-algorithm on the dataset
% m = length(u), n = length(x), x - state vector, u - control vector
    [nSamples, ~] = size(dataset);
    nRules = length(tsModel.Rules);
    % 1. Extract ground truth from dataset.
    X = dataset(:, end-n+1:end);
    input = dataset(:, 1:end-n);
    % 2. Calculating the firing strength of each rule.
    firings = zeros(nSamples, nRules);
    for k=1:nSamples
        [~, ~, ~, ~, ruleFiring] = evalfis(tsModel, input(k, :));
        ruleFiring = ruleFiring / sum(ruleFiring);
        firings(k, :) = ruleFiring;
    end
    % 3. Create Phi = [phi(1)'; ...; phi(n_d)']
%     % comment 115 line and uncomment 113, 114 if you need bias parameter:
%     firings = repelem(firings, 1, n+m+1);
%     input(:, end+1) = 1;  % fake input for bias parameter
    firings = repelem(firings, 1, n+m);
    Phi = repmat(input, 1, nRules) .* firings;
    % 4. bls: use mldivide
    thenParams = Phi \ X;
end

function plotIdentified(tsModel, rhs, method, T, learnStep)
    timesteps = 0:learnStep:T;
    uTest = @(t) 0.02*sin(0.1*pi*t) + 0.15*sin(pi*t) + ...
                 0.2*sin(10*pi*t) + 0.2*sin(100*pi*t);
%     uTest = @(t) 0.5;
%     uTest = @(t) -0.5;
%     uTest = @(t) 0.2 * sin(t);
%     uTest = @(t) -0.2 * sin(t);
%     uTest = @(t) 0.01 * t;
    [~, X_true] = ode45(@(t, x) rhs(x, uTest(t)), timesteps, [0; 0]);
    
    [~, n] = size(X_true);
    X_pred = zeros(length(timesteps), n);     % X_pred(1:2) = x0(1);
    for k=2:length(timesteps)
        X_pred(k, :) = evalfis(tsModel, ...
                               [X_true(k-1, :)'; uTest(timesteps(k-1))]);
    end

    figure()
    title(method)
    nLines = ceil(n/2);
    nColumns = 2;
    for k=1:n
        subplot(nLines, nColumns, k)
        plot(timesteps, X_true(:, k), timesteps, X_pred(:, k))
        legend('true', 'identified')
        title(['x_' num2str(k)])
    end
%     % uncomment in case when you want predict only x1
%     figure()
%     plot(timesteps, X_true(:, 1), timesteps, X_pred(:, 1))
%     legend('true', 'identified')
%     title(method)
end
