function tsIdentification(isLoad, mode)
% identification of TS fuzzy model based on IO data
% isLoad: 1, if load dataset; 0 if create
% mode: 'current' x(k+1) ~ f(x(k), u(k)); 
%       'current_previous' x(k+1) ~ f(x(k), x(k-1), u(k))

    T = 1;
    learnStep = 0.01;

    if nargin == 1
        mode = 'current';     
    end
    rhs = @rhsMotorLink;
    m = 1;      % u = [u_1 ... u_m]';
    n = 2;      % x = [x_1 ... x_n]';
    x0 = [0; -2*pi];
    method = 'SubtractiveClustering';
    if isLoad == 1
        if strcmp(mode, 'current')
            load data/motor_data_current.mat dataset
        else
            load data/motor_data_current_previous.mat dataset
        end
    else
        dataset = collectData(rhs, x0, mode, T, learnStep, m, n);
        dataName = ['data/motor_data_' mode];
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
%     plotIdentified(tsModel, rhs, method, mode, 20, learnStep)
    thenParams = bls(tsModel, dataset, m, n);
    thenParams = addBiasNules(thenParams, nRules, n, m);
    thenParams = reshape(thenParams, 1, []);
    [~, out] = getTunableSettings(tsModel);
    tsModel = setTunableValues(tsModel, out, thenParams);

    % 3. plot and compare
    plotIdentified(tsModel, rhs, method, mode, 20, learnStep)

    % 4. save
    writeFIS(tsModel, ['models/motor_link_ts_' mode])

    % 5. control
    u = tsBasedControl(tsModel, [0; 0], learnStep);
end

function dXdt = rhsMotorLink(x, u)
% right hand size of motor link system
    dXdt = zeros(2, 1);
    dXdt(1) = x(2);
    dXdt(2) = -64*sin(x(1)) - 5*x(2) + 400*u;
end

function dataset = collectData(rhs, x0, mode, T, learnStep, m, n)
% create simulated data for TS model identification
    global uList
    uList = [];
    timesteps = 0:learnStep:T;
    [~, X] = ode45(@(t, x) rhs(x, uRand()), timesteps, x0);
    plot(X(:, 1))
    
    nSteps = length(timesteps);
    if strcmp(mode, 'current')  % x(k+1) ~ f(x(k))
        dataset = zeros(nSteps - 1, m + 2*n);
        for k=1:nSteps-1
            dataset(k, 1:n) = X(k, :);            
            dataset(k, n+1:n+m) = uList(k, :);
            dataset(k, m+n+1:end) = X(k+1, :);
%             % uncomment in case when you want predict only x1
%             dataset(k, 1:n) = X(k, 1);
%             dataset(k, n+1:n+m) = uList(k, :);
%             dataset(k, m+n+1:end) = X(k+1, 1);
        end
    end
end

function u = uRand()
% random uniform control on [-0.5, 0.5]
    u = rand() - 0.5;
    global uList
    uList(end+1, :) = u;       % ? problems with indexing
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

function extParams = addBiasNules(thenParams, nRules, n, m)
% Set bias parameters of tsModel as null.
    extParams = zeros(nRules * (n+m+1), n);
    firstLine = 1;
    for iRule=1:nRules
        extLines = firstLine : firstLine+(n+m-1);
        thenLines = 1+(iRule-1)*(n+m) : iRule*(n+m);
        extParams(extLines, :) = thenParams(thenLines, :);
        firstLine = firstLine + (n+m+1);
    end
end

function res = removeBiasNules(extParams, nRules, n, m)
% Remove bias parameters of tsModel
    linesToRemove = (n+m+1) : (n+m+1) : nRules*(n+m+1);
    extParams(linesToRemove, :) = [];
    res = extParams;
end

function plotIdentified(tsModel, rhs, method, mode, T, learnStep)
    timesteps = 0:learnStep:T;
    uTest = @(t) 0.02*sin(0.1*pi*t) + 0.15*sin(pi*t) + ...
                 0.2*sin(10*pi*t) + 0.2*sin(100*pi*t);
    [~, X_true] = ode45(@(t, x) rhs(x, uTest(t)), timesteps, [0; 0]);
    
    [~, n] = size(X_true);
    X_pred = zeros(length(timesteps), n);     % X_pred(1:2) = x0(1);
    if strcmp(mode, 'current')
        for k=2:length(timesteps)
            X_pred(k, :) = evalfis(tsModel, [X_true(k-1, :)'; ...
                                             uTest(timesteps(k-1))]);
%             % uncomment in case when you want predict only x1
%             X_pred(k, :) = evalfis(tsModel, [X_true(k-1, 1)'; ...
%                                               uTest(timesteps(k-1))]);
        end
    else
        for k=3:length(timesteps)
            X_pred(k, :) = evalfis(tsModel, [X_true(k-1, :)'; ...
                                             X_true(k-2, :)'; ...
                                             uTest(timesteps(k-1))]);
        end
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





