function plotGT(sysName, T, dt, x0)
    % использовал для нахождения параметров testFunctions
    if strcmp(sysName, 'flex2link')
        rhs = @sys.rhsFlex2link;
        wrapper = @sys.flex2linkWrapper;
        interest_idxs = [7, 8];
    end

    timesteps = 0:dt:T;
    nSteps = length(timesteps);
    n = length(x0);
    uTest = testFunctions(sysName);
    nTests = length(uTest);
    X_true = zeros(nTests, nSteps, n);
    for iTest=1:nTests
        % use spline approximation of random control
        uList = arrayfun(uTest{iTest}, timesteps, 'UniformOutput', false);
        uList = cell2mat(uList);
        pp = spline(timesteps, uList);  
        u = @(t) ppval(pp, t);  
        % collect true answers
        [~, X] = ode45(@(t, x) rhs(x, u(t)), timesteps, x0);
        for iStep=1:nSteps
            X(iStep, :) = wrapper(X(iStep, :));
        end
        X_true(iTest, :, :) = X;

        figure()
        plot(X(:, interest_idxs))
    end
    
end