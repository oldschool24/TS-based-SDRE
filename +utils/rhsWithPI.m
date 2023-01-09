function dXdt = rhsWithPI(x, sysName, xIdxPI, uIdxPI, Kp, Ki)
    % x = [state; Integral of error for PIDs]
    dXdt = zeros(size(x));
    r = length(xIdxPI); 
    u = zeros(r, 1);
    I = x(end-r+1 : end);
    for k=1:r
        err = x(xIdxPI(k));   % reference = 0 (stabilisation)
        u(uIdxPI(k)) = Kp(k)*err + Ki(k)*I(k);
        dXdt(end - r + k) = err;
    end

    x = x(1:end-r);   % x = state
    if strcmp(sysName, 'motorLink')
        dXdt(1:end-r) = sys.rhsMotorLink(x, u);
    elseif strcmp(sysName, 'invPend')
        dXdt(1:end-r) = sys.rhsInvPend(x, u);
    elseif strcmp(sysName, 'flex2link')
        dXdt(1:end-r) = sys.rhsFlex2link(x, u);
    end
end
