function res = results(mat_path)
    if nargin == 0
        mat_path = '../results/testInvPend(5-5)';
    end
    load(mat_path, 'criterion')
    nPoints = length(criterion);

    % TODO: new tests for invPend, motorLink
    % success, if integration with TS-control stopped after SDRE
    indSucceed = criterion(:, 11) >= criterion(:, 12);
    
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
    % criterion(:, end-5) contains tsCriterion values, end-4: sdreCriterion
    % criterion(:, 9) contains tsCriterion values, 10: sdreCriterion
    % TODO: change for invPend, motorLink
    relative = 100 * (newCriterion(:, 9) - newCriterion(:, 10)) ...
        ./ newCriterion(:, 10);
    
    [res(2), ind] = min(relative);
    disp(['best point for TS-based: ' ...
        num2str(criterion(new2old(ind), 1:end-6))])

    [res(3), ind] = max(relative);
    disp(['best point for SDRE: ' ...
        num2str(criterion(new2old(ind), 1:end-6))])

    res(4) = mean(relative);
    res(5) = std(relative);
end
