function plotEstimates(name, data, pred, n, timesteps, figName)
    if nargin == 6
        figure('Name', figName)
    else
        figure()
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
end