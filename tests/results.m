function res = results(mat_path)
    if nargin == 0
        mat_path = '../results/testInvPend(5-5)';
    end
    load(mat_path, 'criterion')

    res = zeros(5, 1);
    nPoints = length(criterion);
    indFailed = all(criterion > -1, 2);
    criterion = criterion(indFailed, :);
    res(1) = 100 * length(criterion) / nPoints; % success rate
    relative = 100 * (criterion(:, 1) - criterion(:, 2)) ./ criterion(:, 2);
    res(2) = min(relative);
    res(3) = max(relative);
    res(4) = mean(relative);
    res(5) = std(relative);
end