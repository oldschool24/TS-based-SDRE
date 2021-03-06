function estimatesValidation(sysName, dt, isPlot)
    arguments
        sysName = 'invPend';
        dt = 0.01;
        isPlot = true;
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
    xTrain = trainData(:, 1:n);
    load(['data/', sysName, 'Test'], 'testData')
    xTest = testData(:, 1:n);
    tsModel = readfis(['models/' sysName '.fis']);

    % 2. predict
    [~, fTrain, fTrainPred, B_Train, B_TrainPred] = utils.logger( ...
                                     sysName, xTrain, r, tsModel, dt);
    [~, fTest, fPred, B_Test, B_Pred] = utils.logger(sysName, xTest, ...
                                                     r, tsModel, dt);
    
    % 3. calculate metrics
    dispMetrics('f', fTest, fPred, fTrain, fTrainPred)
    dispMetrics('B', B_Test, B_Pred, B_Train, B_TrainPred)
%     MASE = metrics.MASE(fTest, fPred, fTrain);
%     sMAE = metrics.sMAE(fTest, fPred, fTrain);
%     sMAPE_Train = metrics.sMAPE(fTrain, fTrainPred);
%     sMAPE_Test = metrics.sMAPE(fTest, fPred);
%     disp(['f MASE = ', num2str(MASE)]);
%     disp(['f sMAE = ', num2str(sMAE)]);
%     disp(['f sMAPE on train = ', num2str(sMAPE_Train)]);
%     disp(['f sMAPE on test = ', num2str(sMAPE_Test)]);

%     MASE = metrics.MASE(B_Test, B_Pred, B_Train);
%     sMAE = metrics.sMAE(B_Test, B_Pred, B_Train);
%     sMAPE_Train = metrics.sMAPE(B_Train, B_TrainPred);
%     sMAPE_Test = metrics.sMAPE(B_Test, B_Pred);
%     disp(['B MASE = ', num2str(MASE)]);
%     disp(['B sMAE = ', num2str(sMAE)]);
%     disp(['B sMAPE on train = ', num2str(sMAPE_Train)]);
%     disp(['B sMAPE on test = ', num2str(sMAPE_Test)]);

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
