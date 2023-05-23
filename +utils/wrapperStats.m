function wrapperStats(sysName, tsX, predX, f_true, f_pred, ...
                      B_true, B_pred, imgDir)
    % 1. x density plot
    if strcmp(sysName, 'flex2link')
        idxs = [1, 2];
        r = 2;
        intervals = [-inf, -5*pi; -5*pi -3*pi; -3*pi, -pi; 
                     -pi, pi; pi, 3*pi; 3*pi, 5*pi; 5*pi, inf];
        intervals(:, :, 2) = [-inf, -5*pi; -5*pi -3*pi; -3*pi, -pi; 
                              -pi, pi; pi, 3*pi; 3*pi, 5*pi; 5*pi, inf];
    elseif strcmp(sysName, 'invPend')
        idxs = 1;
        r = 1;
    end
    counts = zeros(size(intervals, 1), length(idxs));
    labels = cell(size(intervals, 1), length(idxs));
    for iComponent=idxs
        data = tsX(:, iComponent);
        for iInterval=1:size(intervals, 1)
            low = intervals(iInterval, 1, iComponent);
            high = intervals(iInterval, 2, iComponent);
            counts(iInterval, iComponent) = sum(data > low & data < high);
            labels{iInterval, iComponent} = sprintf('(%.2f, %.2f)', low, high);
        end
    end
    figure('visible', 'off');
    bar(counts);
    set(gca, 'XTickLabel', labels);
    ylabel('Number of occurrences');
    title('Density of wrapped components');
    exportgraphics(gca, fullfile(imgDir, 'x_density.png'))

    % 2. Plots: dependence of the model error on the x_k, k in idxs
    modelError = abs(tsX(2:end, idxs) - predX(2:end, idxs));
    tsX(end, :) = [];  % x(t=k) <-> abs(x(t=k+1)-pred(t=k+1))
    meanError(tsX, modelError, idxs, 'ts_', imgDir);

    % 3. Same for f
    f_error = abs(f_true(2:end, idxs) - f_pred(2:end, idxs));
    meanError(tsX, f_error, idxs, 'f_', imgDir); 
    
    % 4. Same for B
    for k=1:r
        B_error = abs(B_true(2:end, idxs, k) - B_pred(2:end, idxs, k));
        meanError(tsX, B_error, idxs, ['B_' num2str(k)], imgDir);
    end
end

function meanError(data, error, idxs, figName, imgDir)
    % pred = prediction
    % idxs = indexes of interesting components of state vector
    for iComponent = 1:numel(idxs)
        componentIdx = idxs(iComponent);
        componentName = ['x' num2str(componentIdx)];
%         states = unique(data(:, iComponent)); 
%         meanError = zeros(size(states));
%         for iUnique = 1:numel(states)
%             state = states(iUnique);
%             stateIdxs = (data(:, iComponent) == state);
%             meanError(iUnique) = mean(error(stateIdxs, iComponent));
%         end
        [states, stateIdxs, ~] = unique(data(:, iComponent)); 
        componentError = error(stateIdxs, iComponent);
        figure('Name', [figName num2str(componentIdx) ' error'], ...
               'Visible', 'off')
%         plot(states, meanError)
        plot(states, componentError)
        xlabel(componentName)
        imgName = [figName num2str(componentIdx) '.png'];
        exportgraphics(gca, fullfile(imgDir, imgName))
    end
end
