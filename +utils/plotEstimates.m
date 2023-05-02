function plotEstimates(name, data, pred, n, timesteps, figName, imgDir)
    if nargin == 6
        estFigure = figure('Name', figName);
    else
        estFigure = figure();
    end
    nLines = ceil(n/2);
    nColumns = 2;
    for k=1:n
        subplot(nLines, nColumns, k)
        plot(timesteps, data(:, k), timesteps, pred(:, k))
        legend('true', 'identified', 'FontSize', 14)
        title([name '_' num2str(k)], 'FontSize', 20)
        ax = gca;
        ax.FontSize = 18;
    end
    if nargin == 7
        estFigure.Position = [0, 0, 1920, 1080];
        imgPath = fullfile(imgDir, [figName '.png']);
        exportgraphics(estFigure, imgPath, 'Resolution', 300);
    end
end