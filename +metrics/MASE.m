function res = MASE(yTest, yPred, yTrain)
    % Mean Absolute Scaled Error
    % yTest = ground truth test, yPred = predicted on test
    % yTrain = ground truth train
    meanTestError = mean(abs(yTest - yPred));
    naiveMeanTrainError = mean(abs(yTrain(2:end, :) - yTrain(1:end-1, :)));
    res = meanTestError ./ naiveMeanTrainError;
end
