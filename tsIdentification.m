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
        rng(1, 'twister')   % for reproducibility
        x0Grid = utils.uniformGrid(x0Range, nPoints);
        if withPI   
%             % trajectory -> ss-model -> PI-control -> dataset 
%             trainData = utils.collectData(sysName, x0Grid(1, :), xRange, ...
%                 uRange'/10, T, dt, r, []);
%             % automatic search of Kp, Ki (not reliable)
%             [Kp, Ki] = utils.buildPI(trainData, xIdxPI, uIdxPI, dt, n);

            pidCoefs = struct('Kp', Kp, 'Ki', Ki, 'Kd', Kd, ...
                'xIdxPI', xIdxPI, 'uIdxPI', uIdxPI);
            trainData = utils.collectData(sysName, x0Grid(2:end, :), ...
                xRange, uRange, T, dt, r, pidCoefs);
        else
            trainData = utils.collectData(sysName, x0Grid, xRange, ...
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
%     tic
%     utils.plotIdentified(sysName, extendedModel, 5, dt, x0)
%     toc
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
