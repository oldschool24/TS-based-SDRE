function res = printResults(testStats)
    if isa(testStats, "char")    % if a path is given
        load(testStats, 'testStats')
    end
    x0Grid = testStats('x0');
    [nTests, ~] = size(x0Grid);

    % TODO: new tests for invPend, motorLink
    % success, if integration with TS-control stopped after SDRE
    indSucceed = testStats('tsTime') >= testStats('sdreTime');
    nSucceed = sum(indSucceed);

    removeBadRuns(testStats, indSucceed);
    new2old = zeros(nSucceed, 1);
    iSucceed = 0;
    for iTest=1:nTests
        if indSucceed(iTest)
            iSucceed = iSucceed + 1;
            new2old(iSucceed) = iTest;
        end
    end

    res = zeros(5, 1);
    res(1) = 100 * nSucceed / nTests; % success rate
    % TODO: change for invPend, motorLink
    relative = (testStats('tsCriterion') - testStats('sdreCriterion')) ...
        ./ testStats('sdreCriterion') * 100;
    
    [res(2), ind] = min(relative);
    disp(['best point for TS-based: ' num2str(x0Grid(new2old(ind), :))])
    [res(3), ind] = max(relative);
    disp(['best point for SDRE: ' num2str(x0Grid(new2old(ind), :))])

    res(4) = mean(relative);
    res(5) = std(relative);
end

function removeBadRuns(map, idxs)
    keySet = keys(map);
    for iKey=1:length(keySet)
        key = keySet{iKey};
        column = map(key);
        map(key) = column(idxs, :);
    end
end
