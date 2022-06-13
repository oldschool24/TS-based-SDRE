function res = RMSE(yTest, yPred)
    res = sqrt(mean((yTest-yPred) .^ 2));
end