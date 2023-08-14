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
%         reference(17289:end, :) = [];  % case 16__train 
        reference(750051:758126-1, :) = []; % case test(good) 
        reference = reference(1:step:end, :);
        qR = quantile(reference, [0, 0.25, 0.5, 0.75, 1]); 

%         scoreRef = zscore(reference, 0, 1) * coeff;
        scoreRef = (reference-mu)./sigma * coeff;

%         figure('Name', '2D PCA Prot')
%         plot(score(:, 1), score(:, 2), 'b.')
%         title([num2str(sum(explained(1:2))) ' explained'])
% 
%         figure('Name', '2D PCA Ref')
%         plot(scoreRef(:, 1), scoreRef(:, 2), 'r.')

        figure('Name', '2D PCA')
        plot(score(:, 1), score(:, 2), 'b.', scoreRef(:, 1), scoreRef(:, 2), 'r.')
        title([num2str(sum(explained(1:2))) ' explained'])

        black_score =  [ 3.89, -2.97,  2.99;
                         2.12, -3.20,  3.89;
                        -0.64, -3.55,  3.41;
                        -5.43, -0.89,  0.08;
                        -5.98, -1.16,  0.70;
                         4.70,  0.42, -1.68;
                         5.35,  0.61, -1.06;
                        -0.11,  3.66, -1.21;
                        -2.52, -2.73, -0.18; % new
                         1.53,  2.83, -0.64;
                        -1.34, -1.60, -2.36;
                         2.02, -3.62,  2.28;
                        -1.31, -1.47, -2.27;
                         2.61, -1.92, 2.11;
                         0.79, 2.72, -0.39];
%         [~, idxs] = ismembertol(black_score, scoreRef(:, 1:3), 2e-2, ...
%                                 "ByRows", true);
%         yellow_score = scoreRef(idxs, :);        
%         % get the coordinates of the yellow points in the original space
%         refYellow = reference(idxs, :);

        figure('Name', '3D PCA')
        scatter3(score(:, 1), score(:, 2), score(:, 3), 'blue')
        hold on
        scatter3(scoreRef(:, 1), scoreRef(:, 2), scoreRef(:, 3), 'red')
        scatter3(black_score(:, 1), black_score(:, 2), ...
                 black_score(:, 3), 'black', 'filled')
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