function q = analyzeDataset(dataPath)
    load(dataPath, 'valData')
    valData(:, 9:end) = [];
    
    q = quantile(valData, [0, 0.25, 0.5, 0.75, 1]);

    valData = valData(1:20:end, :);
    [coeff, score, ~, ~, explained, ~] = pca(zscore(valData, 0, 1));
    
    figure('Name', '2D PCA')
    plot(score(:, 1), score(:, 2), '.')
    title([num2str(sum(explained(1:2))) ' explained'])
    
    figure('Name', '3D PCA')
    scatter3(score(:, 1), score(:, 2), score(:, 3))
    title([num2str(sum(explained(1:3))) ' explained'])
    
    disp('Which variables explain the clustering')
    disp(coeff(:, 1:3)')
end