function mainSim(sysRhs, learnStep, T, tsModel, x0)
    arguments
        sysRhs = @sys.rhsMotorLink;
        learnStep = 0.01;
        T = 20;
        tsModel = readfis('models/motor_link.fis'); 
        x0 = [-pi/2; -3*pi];  % change it sometime
    end

    if isequal(sysRhs, @sys.rhsMotorLink)
        sysName = 'motorLink';
        wrapper = @(x) x;
    elseif isequal(sysRhs, @sys.rhsInvPend)
        sysName = 'invPend';
        wrapper = @sys.invPendWrapper;
    end

    % 1.1 integrate
    u = @(x) tsBasedControl(x, sysName, tsModel, learnStep);
    timesteps = 0:learnStep:T;
    [~, X_true] = ode45(@(t, x) sysRhs(x, u(x)), timesteps, x0);
    % 1.2 process data for plot
    nSteps = length(timesteps);
    for iStep=1:nSteps
        X_true(iStep, :) = wrapper(X_true(iStep, :));
    end

    % 2. Calculate u and estimates(f, B) at timesteps
    n = length(x0);
    r = length(u(X_true(1, :)'));
    [uList, fTrue, fPred, Btrue, Bpred] = utils.logger(sysName, X_true, ...
                                                       nSteps, n, r, ...
                                                       tsModel, learnStep);
%     nSteps = length(timesteps);
%     n = length(x0);
%     r = length(u(X_true(1, :)'));
%     uList = zeros(nSteps, r);
%     fTrue = zeros(nSteps, n);
%     fPred = zeros(nSteps, n);
%     Btrue = zeros(nSteps, n);
%     Bpred = zeros(nSteps, n);
%     for iStep=1:nSteps
%         X_true(iStep, :) = wrapper(X_true(iStep, :));
%         x = X_true(iStep, :)';
%         [u, fHat, Bhat] = tsBasedControl(x, sysName, tsModel, learnStep);
%         uList(iStep, :) = u;
%         fPred(iStep, :) = fHat;
%         Bpred(iStep, :) = Bhat;
%         % calculate true f, B
%         if strcmp(sysName, 'motorLink')
%             fTrue(iStep, :) = [x(2), -64*sin(x(1)) - 5*x(2)];
%             Btrue(iStep, :) = [0, 400];
%         elseif strcmp(sysName, 'invPend')
%             denominator = (M + m)*I + (m*l*sin(x(2)))^2 + M*m*l^2;
%             fTrue(iStep, :) = [x(3), x(4), ...
%                 (-b*(I + m*l^2)*x(3) - m^2*l^3*sin(x(2))*x(4)^2 + ...
%                 g*m^2*l^2*sin(x(2))*cos(x(2)) - I*m*l*sin(x(2))*x(4)^2) ...
%                 / denominator, ...
%                 -m*l*(m*l*sin(x(2))*cos(x(2))*x4^2 + b*cos(x(2))*x(3) - ...
%                 (M+m)*g*sin(x(2))) / denominator];
%             Btrue(iStep, :) = [0, 0, ...
%                 (I + m * l^2) / denominator, ...
%                 -m*l*cos(x(2)) / denominator];
%         end
%     end

    % 3. Plot u, estimates and trajectories
    plot(timesteps, uList)
    labels = cell(n, 1);
    for k=1:n
        figure()
        title(['f' num2str(k)])
        plot(timesteps, fTrue(:, k), timesteps, fPred(:, k))
        legend(['f' num2str(k)], ['f' num2str(k) 'Pred'])

        figure()
        title(['B' num2str(k)])
        plot(timesteps, Btrue(:, k), timesteps, Bpred(:, k))
        legend(['B' num2str(k)], ['B' num2str(k) 'Pred'])

        labels{k} = ['x' num2str(k)];
    end
    figure()
    plot(timesteps, X_true)
    legend(labels)
end
