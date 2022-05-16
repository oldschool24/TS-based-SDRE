function [uList, fTrue, fPred, Btrue, Bpred] = logger(sysName, X, nSteps, ...
                                                      n, r, tsModel, learnStep)
    if strcmp(sysName, 'invPend')
        M = 0.5;
        m = 0.2;  % PROBLEM: same as len(u)
        b = 0.1;
        l = 0.3;
        I = 0.006;
        g = 9.8;
    end

    uList = zeros(nSteps, r);
    fTrue = zeros(nSteps, n);
    fPred = zeros(nSteps, n);
    Btrue = zeros(nSteps, n);
    Bpred = zeros(nSteps, n);
    for iStep=1:nSteps
        x = X(iStep, :)';
        [u, fHat, Bhat] = tsBasedControl(x, sysName, tsModel, learnStep);
        uList(iStep, :) = u;
        fPred(iStep, :) = fHat;
        Bpred(iStep, :) = Bhat;
        % calculate true f, B
        if strcmp(sysName, 'motorLink')
            fTrue(iStep, :) = [x(2), -64*sin(x(1)) - 5*x(2)];
            Btrue(iStep, :) = [0, 400];
        elseif strcmp(sysName, 'invPend')
            denominator = (M + m)*I + (m*l*sin(x(2)))^2 + M*m*l^2;
            fTrue(iStep, :) = [x(3), x(4), ...
                (-b*(I + m*l^2)*x(3) - m^2*l^3*sin(x(2))*x(4)^2 + ...
                g*m^2*l^2*sin(x(2))*cos(x(2)) - I*m*l*sin(x(2))*x(4)^2) ...
                / denominator, ...
                -m*l*(m*l*sin(x(2))*cos(x(2))*x4^2 + b*cos(x(2))*x(3) - ...
                (M+m)*g*sin(x(2))) / denominator];
            Btrue(iStep, :) = [0, 0, ...
                (I + m * l^2) / denominator, ...
                -m*l*cos(x(2)) / denominator];
        end
    end
end
