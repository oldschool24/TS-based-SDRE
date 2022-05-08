function mainSim(sysRhs, learnStep, T, tsModel)
    arguments
        sysRhs = @sys.rhsMotorLink;
        learnStep = 0.01;
        T = 20;
        tsModel = readfis('models/motor_link.fis'); 
    end

    u = @(x) tsBasedControl(x, tsModel, learnStep);
    timesteps = 0:learnStep:T;
    x0 = [-pi/2; -3*pi];  % change it sometime
    global fTrue fPred Btrue Bpred
    [~, X_true] = ode45(@(t, x) sysRhs(x, u(x)), timesteps, x0);
    
    nSteps = length(timesteps);
    m = length(u(X_true(1, :)'));
    uList = zeros(nSteps, m);
    for iStep=1:nSteps
        uList(nSteps, :) = u(X_true(1, :)');
    end
    plot(timesteps, uList)

    figure()
    plot(timesteps, X_true)
    legend('x1', 'x2')

    L = length(fTrue);
    figure()
    title('f1')
    plot(1:L, fTrue(:, 1), 1:L, fPred(:, 1))
    legend('f_1', 'f_1Pred')

    figure()
    title('f2')
    plot(1:L, fTrue(:, 2), 1:L, fPred(:, 2))
    legend('f_2', 'f_2Pred')

    figure()
    title('B1')
    plot(1:L, Btrue(:, 1), 1:L, Bpred(:, 1))
    legend('B_1', 'B_1Pred')

    figure()
    title('B2')
    plot(1:L, Btrue(:, 2), 1:L, Bpred(:, 2))
    legend('B_2', 'B_2Pred')
end
