function res = sMAPE(yTest, yPred)
    % symmetric Mean Absolute Percentage Error
    % yTest = ground truth test, yPred = predicted values
    res = abs(yPred - yTest) ./ (abs(yPred) + abs(yTest));
    res = 100 * mean(res);        % in percents
end