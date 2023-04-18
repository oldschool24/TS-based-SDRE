function testData = collectValData(sysName, dt, T, xRange, nPoints, reduction)
    arguments
        sysName
        dt double {mustBePositive}
        T = 10
        xRange = [-6, -pi/2, -12, -10; ...
                   6,  pi/2,  12,  10]
        nPoints = 10
        reduction = 0.05
    end
    
    if strcmp(sysName, 'motorLink')
        r = 1;
        rhs = @sys.rhsMotorLink;
        wrapper = @(x) x;
    elseif strcmp(sysName, 'invPend')
        r = 1;
        rhs = @sys.rhsInvPend;
        wrapper = @sys.invPendWrapper;
    elseif strcmp(sysName, 'flex2link')
        r = 2;
        rhs = @sys.rhsFlex2link;
        wrapper = @sys.flex2linkWrapper;
    end

    uTest = testFunctions(sysName); % test control functions
    nTests = length(uTest);

    [~, n] = size(xRange);
    x0Range = utils.reduceRange(xRange, reduction);
    x0Grid = utils.uniformGrid(x0Range, nPoints); % set of initial conditions
    
    timesteps = 0:dt:T;
    nSteps = length(timesteps);
    odeOpt = odeset('Events', ...
                    @(t, x) utils.xRangeEvent(wrapper(x), xRange));
    testData = zeros(nTests * nPoints * (nSteps-1), 2*n + r);
    iData = 0;
    for iTest=1:nTests
        % use spline approximation of random control
        uList = arrayfun(uTest{iTest}, timesteps, 'UniformOutput', false);
        uList = cell2mat(uList);
        pp = spline(timesteps, uList);  
        u = @(t) ppval(pp, t);  
        % collect true answers       
        for iPoint=1:nPoints
            % integrate with initial condition = x0
            x0 = x0Grid(iPoint, :)';
            [t, X] = ode45(@(t, x) rhs(x, u(t)), timesteps, x0, odeOpt);
            
            % save data from simulation
            nSteps = length(t);
            for iStep=2:nSteps
                testData(iData+iStep-1, 1:n) = wrapper(X(iStep-1, :)');
                testData(iData+iStep-1, n+1:n+r) = uList(:, iStep-1);
                testData(iData+iStep-1, r+n+1:end) = wrapper(X(iStep, :)');
            end
            iData = iData + (nSteps - 1);
        end
    end
    testData(iData+1:end, :) = [];
end
