function res = results(mat_path)
    if nargin == 0
        mat_path = '../results/testInvPend(5-5)';
    end
    load(mat_path, 'criterion')

    nPoints = length(criterion);
    indSucceed = all(criterion(:, end-1:end) > -1, 2);
    newCriterion = criterion(indSucceed, :);
    nSucceed = sum(indSucceed);
    new2old = zeros(nSucceed, 1);
    nFails = 0;
    for k=1:nSucceed
        if ~indSucceed(k)
            nFails = nFails + 1;
        end
        new2old(k) = k + nFails;
    end

    res = zeros(5, 1);
    res(1) = 100 * nSucceed / nPoints; % success rate
    relative = 100 * (newCriterion(:, end-1) - newCriterion(:, end)) ...
        ./ newCriterion(:, end);
    
    [res(2), ind] = min(relative);
    disp(['best point for TS-based: ' ...
        num2str(criterion(new2old(ind), 1:end-2))])

    [res(3), ind] = max(relative);
    disp(['best point for SDRE: ' ...
        num2str(criterion(new2old(ind), 1:end-2))])

    res(4) = mean(relative);
    res(5) = std(relative);
end