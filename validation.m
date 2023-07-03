function validation(valConfigPath, isPlot, isFast)
    % 0. load settings
    valConfig = jsondecode(fileread(valConfigPath));
    expPath = fullfile('runs', valConfig.expName);
    copyfile(valConfigPath, fullfile(expPath, 'valData.json'))
    trainConfigPath = fullfile(expPath, 'trainData.json');
    trainConfig = jsondecode(fileread(trainConfigPath));
    n = trainConfig.n;
    r = trainConfig.r;
    sysName = trainConfig.sysName;
    dt = trainConfig.dt;
    isWrap = trainConfig.isWrap;

    if isempty(valConfig.valDataName)
        T = valConfig.T;
        xRange = valConfig.xRange;
        nPoints = valConfig.nPoints;
        reduction = valConfig.reduction;
        valData = utils.collectValData(sysName, dt, T, xRange, ...
                                       nPoints, reduction, isWrap);
        valPath = fullfile('data', valConfig.expName, 'data');
        mkdir(fullfile('data', valConfig.expName))
        save(valPath, 'valData')
    else
        valPath = fullfile('data', valConfig.valDataName, 'data');
    end
    Q = valConfig.qCriterion * eye(n);
    R = valConfig.rCriterion * eye(r);
    isAnalyze = valConfig.isAnalyze;

    % 1. load model and data
    load(fullfile(expPath, 'model'), "extendedModel")
    load(fullfile(expPath, 'trainData'), 'trainData')
    load(valPath, 'valData')
    if isFast  % use only 25% of data
        trainData = trainData(1:4:end, :);  
        valData = valData(1:4:end, :);
%         nSamples = min(size(trainData, 1), size(valData, 1));
%         nSamples = min(nSamples, 100000);
%         trainData = trainData(1:nSamples, :);
%         valData = valData(1:nSamples, :);
    end
    xTrain = trainData(:, end-n+1:end);
    xVal = valData(:, end-n+1:end);
    tsModel = extendedModel.model;
    modelRange = extendedModel.range;

    % 2. predict
    % turn off warnings
    warning('off', 'fuzzy:general:warnEvalfis_NoRuleFired')
    warning('off', 'fuzzy:general:diagEvalfis_OutOfRangeInput')
    parfevalOnAll(@warning,0,'off','fuzzy:general:warnEvalfis_NoRuleFired');
    parfevalOnAll(@warning,0,'off','fuzzy:general:diagEvalfis_OutOfRangeInput');
    % predictions
    xTrainPred = modelPrediction(trainData(:, 1:n+r), tsModel, ...
                                 modelRange, isWrap, sysName);
    disp('Prediction on the training data is done')
    xPred = modelPrediction(valData(:, 1:n+r), tsModel, ...
                            modelRange, isWrap, sysName);
    disp('Prediction on the validation data is done')
    [~, f_train, f_trainPred, B_train, B_trainPred] = utils.logger( ...
        sysName, xTrain, r, extendedModel, dt, [], Q, R, isWrap, true);
    disp('Estimates on the training data are calculated')
    [~, f_val, f_pred, B_val, B_pred] = utils.logger( ...
        sysName, xVal, r, extendedModel, dt, [], Q, R, isWrap, true);
    disp('Estimates on the validation data are calculated')
    
    % 3. calculate metrics
    ts_results = calcMetrics(xVal, xPred, xTrain, xTrainPred);
    f_results = calcMetrics(f_val, f_pred, f_train, f_trainPred);
    B_results = cell(r, 1);
    for iControl=1:r
        B_results{iControl} = calcMetrics( ...
            B_val(:, :, iControl), B_pred(:, :, iControl), ...
            B_train(:, :, iControl), B_trainPred(:, :, iControl));
    end
    save(fullfile(expPath, 'ts'), "ts_results")
    save(fullfile(expPath, 'f'), "f_results")
    save(fullfile(expPath, 'B'), "B_results")

    % 4. Plots: error(t)
    if isPlot
        [nSamplesTrain, ~] = size(xTrain);
        [nSamplesVal, ~] = size(xVal);
        utils.plotEstimates('x', xTrain, xTrainPred, n, ...
                            1:nSamplesTrain, 'x on train')
        utils.plotEstimates('x', xVal, xPred, n, ...
                            1:nSamplesVal, 'x on test')
        utils.plotEstimates('f', f_train, f_trainPred, n, ...
                            1:nSamplesTrain, 'f on train')
        utils.plotEstimates('f', f_val, f_pred, n, ...
                            1:nSamplesVal, 'f on test')
        utils.plotEstimates('B', B_train, B_trainPred, n, ...
                            1:nSamplesTrain, 'B on train')
        utils.plotEstimates('B', B_val, B_pred, n, ...
                            1:nSamplesVal, 'B on test')
    end

    % 5. Plots: error(x)
    if isAnalyze
        % 5.1 On train
        analyzeDir = fullfile(expPath, "trainAnalysis");
        if ~exist(analyzeDir, 'dir')
            mkdir(analyzeDir)
        end
        utils.wrapperStats(sysName, xTrain, xTrainPred, f_train, ...
                           f_trainPred, B_train, B_trainPred, analyzeDir)
        
        % 5.2 On validation
        analyzeDir = fullfile(expPath, "valAnalysis");
        if ~exist(analyzeDir, 'dir')
            mkdir(analyzeDir)
        end
        utils.wrapperStats(sysName, xVal, xPred, f_val, f_pred, ...
                           B_val, B_pred, analyzeDir)
    end
end

function yPred = modelPrediction(y, tsModel, modelRange, isWrap, sysName)
    [nSamples, ~] = size(y);
    [~, n] = size(modelRange);
    yPred = zeros(nSamples, n);
    yPred(1, :) = y(1, 1:n);
    parfor iSample=2:nSamples
        yPred(iSample, :) = utils.evalProjection( ...
            tsModel, y(iSample-1, :)', modelRange, isWrap, sysName);
    end
%     tic 
%     yCell = mat2cell(y(1:end-1, :), ones(nSamples-1, 1));
%     yPred = cellfun(@(x) utils.evalProjection(tsModel, x', modelRange, ...
%                                               isWrap, sysName), ...
%                     yCell, 'UniformOutput', false);
%     yPred = cell2mat([{y(1, 1:n)}; yPred]);
%     toc
end

function res = calcMetrics(yVal, yPred, yTrain, yTrainPred)
    res.MASE = metrics.MASE(yVal, yPred, yTrain);
    res.sMAE = metrics.sMAE(yVal, yPred, yTrain);
    res.sMAPE_Train = metrics.sMAPE(yTrain, yTrainPred);
    res.sMAPE_Val = metrics.sMAPE(yVal, yPred);
    res.RMSE_Train = metrics.RMSE(yTrain, yTrainPred);
    res.RMSE_Val = metrics.RMSE(yVal, yPred);
end
