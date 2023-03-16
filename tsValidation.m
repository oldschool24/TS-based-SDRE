function [MASE, sMAE, sMAPE] = tsValidation(modelPath, trainPath, testPath)
    if contains(trainPath, 'invPend')
        n = 4;
        r = 1;
    elseif contains(trainPath, 'motorLink')
        n = 2;
        r = 1;
    elseif contains(trainPath, 'flex2link')
        n = 8;
        r = 2;
    end

    % 1. load data and model
    load(trainPath, 'trainData')
    xTrain = trainData(:, end-n+1:end);
    load(testPath, 'testData')
    xTest = testData(:, end-n+1:end);
    load(modelPath, "extendedModel")
    tsModel = extendedModel.model;
    modelRange = extendedModel.range;
    
    % 2. predict
    [nSamples, ~] = size(testData);
    xPred = zeros(nSamples, n);
    for iSample=1:nSamples
        xPred(iSample, :) = utils.evalProjection( ...
            tsModel, testData(iSample, 1:n+r)', modelRange);
    end
    
    % 3. calculate metrics
    MASE = metrics.MASE(xTest, xPred, xTrain);
    sMAE = metrics.sMAE(xTest, xPred, xTrain);
    sMAPE = metrics.sMAPE(xTest, xPred);
    disp(['test MASE = ', num2str(MASE)]);
    disp(['test sMAE = ', num2str(sMAE)]);
    disp(['test sMAPE = ', num2str(sMAPE)]);
end