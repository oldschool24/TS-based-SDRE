function savePlots(data, imgDir, figName)
    fig = figure('name', figName);
    scatter(1:length(data), data, 'filled')

    iComponent = figName(2:end);
    if isempty(iComponent)
        varName = figName;
    else
        varName = strrep(figName, iComponent, ['_{.' iComponent '}']);
    end

    xlabel('Trial number')
    ylabel(sprintf('$\\frac{ ||%s(x) - \\hat %s(x)|| }{ ||%s(x)|| }$', ...
                   varName, varName, varName), ...
           'Interpreter', 'latex')
    xlim([0, length(data)+1])
    yRange = max(data) - min(data);
    ylim([min(data) - 0.05*yRange, max(data) + 0.05*yRange])
    ax = gca;
    ax.FontSize = 18;

    fig.Position = [0, 0, 1920, 1080];
    imgPath = fullfile(imgDir, [figName '_error.png']);
    exportgraphics(fig, imgPath, 'Resolution', 300);
end
