function [qP, qR, refYellow] = compareDatasets(prototypePath, referencePath)
    prototype = importdata(prototypePath);
    prototype(:, 9:end) = [];
    step = 10;

    qP = quantile(prototype, [0, 0.25, 0.5, 0.75, 1]);
    prototype = prototype(1:step:end, :);
    [~, mu, sigma] = zscore(prototype, 0, 1);
    [coeff, score, ~, ~, explained, ~] = pca(zscore(prototype, 0, 1));
    
    if exist("referencePath", "var")
        reference = importdata(referencePath);
        reference(:, 9:end) = [];
        reference(750051:758126-1, :) = []; % case test-good_1 
        reference = reference(1:step:end, :);
        qR = quantile(reference, [0, 0.25, 0.5, 0.75, 1]); 

        scoreRef = (reference-mu)./sigma * coeff;

        figure('Name', '2D PCA')
        plot(score(:, 1), score(:, 2), 'b.', scoreRef(:, 1), scoreRef(:, 2), 'r.')
        title([num2str(sum(explained(1:2))) ' explained'])

        figure('Name', '3D PCA')
        scatter3(score(:, 1), score(:, 2), score(:, 3), 'blue')
        hold on
        scatter3(scoreRef(:, 1), scoreRef(:, 2), scoreRef(:, 3), 'red')
        title([num2str(sum(explained(1:3))) ' explained'])
    else
        figure('Name', '2D PCA')
        plot(score(:, 1), score(:, 2), '.')
        title([num2str(sum(explained(1:2))) ' explained'])
        
        figure('Name', '3D PCA')
        scatter3(score(:, 1), score(:, 2), score(:, 3))
        title([num2str(sum(explained(1:3))) ' explained'])
        
        disp('Which variables explain the clustering')
        disp(coeff(:, 1:3)')
    end
end