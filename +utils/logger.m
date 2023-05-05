function [uList, fTrue, fPred, Btrue, Bpred] = logger( ...
    sysName, X, r, extendedModel, dt, known, Q, R, isWrap)

    [nSteps, n] = size(X);
    uList = zeros(nSteps, r);
    fTrue = zeros(nSteps, n);
    fPred = zeros(nSteps, n);
    Btrue = zeros(nSteps, n, r);
    Bpred = zeros(nSteps, n, r);
    for iStep=1:nSteps
        x = X(iStep, :)';
        [u, fHat, Bhat] = tsBasedControl(x, extendedModel, sysName, ...
                                         dt, known, Q, R, isWrap);
        uList(iStep, :) = u;
        fPred(iStep, :) = fHat;
        Bpred(iStep, :, :) = Bhat;
        fTrue(iStep, :) = sys.get_f(x, sysName);
        Btrue(iStep, :, :) = sys.get_B(x, sysName);
    end
end
