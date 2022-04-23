function mainSim(sysRhs, learnStep, T, tsModel)
    arguments
        sysRhs = @rhsMotorLink;
        learnStep = 0.01;
        T = 20;
        tsModel = readfis('models/motor_link_ts_current.fis'); 
    end

    u = @(x) tsBasedControl(x, tsModel, learnStep);
    timesteps = 0:learnStep:T;
    x0 = [-0.4; 2];  % change it sometime
    [~, X_true] = ode45(@(t, x) sysRhs(x, u(x)), timesteps, x0);
    
    figure()
    plot(X_true')
end

function dXdt = rhsMotorLink(x, u)
% right hand size of motor link system
    dXdt = zeros(2, 1);
    dXdt(1) = x(2);
    dXdt(2) = -64*sin(x(1)) - 5*x(2) + 400*u;
end

