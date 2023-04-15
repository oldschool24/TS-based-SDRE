function estimatesValidation(modelPath, trainPath, testPath, dt, isPlot)
    arguments
        modelPath {mustBeText}
        trainPath {mustBeText}
        testPath {mustBeText}
        dt double {mustBePositive} 
        isPlot = true
    end

    if contains(trainPath, 'invPend')
        n = 4;
        r = 1;
        sysName = 'invPend';
    elseif contains(trainPath, 'motorLink')
        n = 2;
        r = 1;
        sysName = 'motorLink';
    elseif contains(trainPath, 'flex2link')
        n = 8;
        r = 2;
        sysName = 'flex2link';
    end

    % 1. load data and model
    load(trainPath, 'trainData')
    xTrain = trainData(:, end-n+1:end);
    load(testPath, 'testData')
    xTest = testData(:, end-n+1:end);
    load(modelPath, "extendedModel")

    % 2. predict
    [~, fTrain, fTrainPred, B_Train, B_TrainPred] = utils.logger( ...
                                     sysName, xTrain, r, extendedModel, dt);
    [~, fTest, fPred, B_Test, B_Pred] = utils.logger(sysName, xTest, ...
                                                     r, extendedModel, dt);
    
    % 3. calculate metrics
    dispMetrics('f', fTest, fPred, fTrain, fTrainPred)
    for iControl=1:r
        dispMetrics(['B^' num2str(iControl)], ...
            B_Test(:, :, iControl), B_Pred(:, :, iControl), ...
            B_Train(:, :, iControl), B_TrainPred(:, :, iControl))
    end

    % 4. Plot
    [nSamplesTest, ~] = size(xTest);
    [nSamplesTrain, ~] = size(xTrain);
    if isPlot
        utils.plotEstimates('f', fTrain, fTrainPred, n, 1:nSamplesTrain, 'f on train')
        utils.plotEstimates('f', fTest, fPred, n, 1:nSamplesTest, 'f on test')
        utils.plotEstimates('B', B_Train, B_TrainPred, n, 1:nSamplesTrain, 'B on train')
        utils.plotEstimates('B', B_Test, B_Pred, n, 1:nSamplesTest, 'B on test')
    end
end

function dispMetrics(name, yTest, yPred, yTrain, yTrainPred)
    MASE = metrics.MASE(yTest, yPred, yTrain);
    sMAE = metrics.sMAE(yTest, yPred, yTrain);
    sMAPE_Train = metrics.sMAPE(yTrain, yTrainPred);
    sMAPE_Test = metrics.sMAPE(yTest, yPred);
    RMSE_Train = metrics.RMSE(yTrain, yTrainPred);
    RMSE_Test = metrics.RMSE(yTest, yPred);

    disp([name, ' MASE = ', num2str(MASE)]);
    disp([name, ' sMAE = ', num2str(sMAE)]);
    disp([name, ' sMAPE on train = ', num2str(sMAPE_Train)]);
    disp([name, ' sMAPE on test = ', num2str(sMAPE_Test)]);
    disp([name, ' RMSE on train = ', num2str(RMSE_Train)]);
    disp([name, ' RMSE on test = ', num2str(RMSE_Test)]);
end
