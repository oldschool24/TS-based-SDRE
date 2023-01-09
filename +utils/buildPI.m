function [Kp, Ki] = buildPI(data, xIdxPI, uIdxPI, dt, n, order, isPlot)
    arguments
        data
        xIdxPI
        uIdxPI
        dt {mustBePositive}
        n {mustBePositive}
        order {mustBePositive} = 1
        isPlot {mustBeNumericOrLogical} = false
    end

    nRegulators = length(xIdxPI); 
    Kp = zeros(nRegulators, 1);
    Ki = zeros(nRegulators, 1);
    for iRegulator=1:nRegulators
        y = data(:, xIdxPI(iRegulator));
        u = data(:, n+uIdxPI(iRegulator));

        estimationData = iddata(y, u, dt);
        sys_id = ssest(estimationData, order);
        [~, fit, ~] = compare(estimationData, sys_id); 
        disp(['Goodness of fit (ss for u_' num2str(iRegulator) ') ' ...
            'value: ' num2str(fit)])
        
        % PID tuning algorithm for linear plant model
        [pidControl, ~] = pidtune(sys_id, 'PI');
        Kp(iRegulator) = pidControl.Kp;
        Ki(iRegulator) = pidControl.Ki;

        if isPlot
            Response = getPIDLoopResponse(pidControl, sys_id, 'closed-loop');
            stepplot(Response)
            title('Step Plot: Reference tracking')
            grid on
        end
    end
end
