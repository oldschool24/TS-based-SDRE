function dXdt = rhsWithPI(x, sysName, xIdxPI, uIdxPI, Kp, Ki, Kd)
    % x = [state; Integral of error for PIDs]
    dXdt = zeros(size(x));
    r = length(xIdxPI); 
    
    if strcmp(sysName, 'flex2link')
        xdotIdx = [5, 6, 7, 8, -1, -1, -1, -1]; % x1_dot = x5, ...
    else
        xdotIdx = [];
    end

    u = zeros(r, 1);
    I = x(end-r+1 : end);
    for k=1:r
        componentNum = xIdxPI(k);
        err = x(componentNum);   % reference = 0 (stabilisation)
        u(uIdxPI(k)) = Kp(k)*err + Ki(k)*I(k);
        dXdt(end - r + k) = err;
        if ~isempty(xdotIdx)
            dotComponentNum = xdotIdx(componentNum);
            if dotComponentNum ~= -1
                errDot = x(dotComponentNum);
                u(uIdxPI(k)) = u(uIdxPI(k)) + Kd(k) * errDot;
            end
        end
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
