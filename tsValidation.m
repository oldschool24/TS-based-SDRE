function [MASE, sMAE, sMAPE] = tsValidation(sysName)
    arguments
        sysName = 'invPend'
    end

    if strcmp(sysName, 'invPend')
        n = 4;
        r = 1;
    elseif strcmp(sysName, 'motorLink')
        n = 2;
        r = 1;
    end

    % 1. load data and model
    load(['data/', sysName], 'trainData')
    xTrain = trainData(:, end-n+1:end);
    load(['data/', sysName, 'Test'], 'testData')
    xTest = testData(:, end-n+1:end);
    tsModel = readfis(['models/' sysName '.fis']);
    
    % 2. predict
    [nSamples, ~] = size(testData);
    xPred = zeros(nSamples, n);
    for iSample=1:nSamples
        xPred(iSample, :) = evalfis(tsModel, testData(iSample, 1:n+r)');
    end
    
    % 3. calculate metrics
    MASE = metrics.MASE(xTest, xPred, xTrain);
    sMAE = metrics.sMAE(xTest, xPred, xTrain);
    sMAPE = metrics.sMAPE(xTest, xPred);
    disp(['test MASE = ', num2str(MASE)]);
    disp(['test sMAE = ', num2str(sMAE)]);
    disp(['test sMAPE = ', num2str(sMAPE)]);
end