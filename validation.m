function validation(valConfigPath, isPlot)
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

    % 1. load model and data
    load(fullfile(expPath, 'model'), "extendedModel")
    load(fullfile(expPath, 'trainData'), 'trainData')
    xTrain = trainData(:, end-n+1:end);
    load(valPath, 'valData')
    xVal = valData(:, end-n+1:end);
    tsModel = extendedModel.model;
    modelRange = extendedModel.range;

    % 2. predict
    xTrainPred = modelPrediction(trainData(:, 1:n+r), tsModel, ...
                                 modelRange, isWrap, sysName);
    xPred = modelPrediction(valData(:, 1:n+r), tsModel, ...
                            modelRange, isWrap, sysName);
    [~, f_train, f_trainPred, B_train, B_trainPred] = utils.logger( ...
        sysName, xTrain, r, extendedModel, dt, [], Q, R, isWrap);
    [~, f_val, f_pred, B_val, B_pred] = utils.logger( ...
        sysName, xVal, r, extendedModel, dt, [], Q, R, isWrap);
    
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

    % 4. Plot
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
end

function yPred = modelPrediction(y, tsModel, modelRange, isWrap, sysName)
    [nSamples, ~] = size(y);
    [~, n] = size(modelRange);

    yPred = zeros(nSamples, n);
    yPred(1, :) = y(1, 1:n);
    for iSample=2:nSamples
        yPred(iSample, :) = utils.evalProjection( ...
            tsModel, y(iSample-1, :)', modelRange, isWrap, sysName);
    end
end

function res = calcMetrics(yVal, yPred, yTrain, yTrainPred)
    res.MASE = metrics.MASE(yVal, yPred, yTrain);
    res.sMAE = metrics.sMAE(yVal, yPred, yTrain);
    res.sMAPE_Train = metrics.sMAPE(yTrain, yTrainPred);
    res.sMAPE_Val = metrics.sMAPE(yVal, yPred);
    res.RMSE_Train = metrics.RMSE(yTrain, yTrainPred);
    res.RMSE_Val = metrics.RMSE(yVal, yPred);
end
