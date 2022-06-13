function res = sMAE(yTest, yPred, yTrain)
    % scaled Mean Absolute Error
    % yTest = ground truth test, yPred = predicted on test
    % yTrain = ground truth train
    meanTestError = mean(abs(yTest - yPred));
    meanTrain = mean(abs(yTrain));
    res = 100 * meanTestError ./ meanTrain;  % in percents
end