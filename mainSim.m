function [tsCriterion, sdreCriterion] = mainSim(sysName, dt, T, tsModel, ...
                                                x0, Q, R)
    arguments
        sysName = 'motorLink';
        dt = 0.01;
        T = 0.5;
        tsModel = readfis('models/motorLink.fis');
        x0 = [-pi/2; -3*pi];  % change it sometime
        Q = 10 * eye(2);
        R = 5;
    end

    % 0.1 set vars for model
    if strcmp(sysName, 'motorLink')
        wrapper = @(x) x;
        n = 2;
        r = 1;
        known = [];
    elseif strcmp(sysName, 'invPend')
        wrapper = @sys.invPendWrapper;
        n = 4;
        r = 1;
        known = [];
    elseif strcmp(sysName, 'flex2link')
        wrapper = @sys.flex2linkWrapper;
        n = 8;
        r = 2;
        known = [];
    end
    tsFailed = false;  % var for event detection
    
    % 0.2 add calculation of criterion to rhs
    function dXdt = rhsWithCriterion(x, uName, sysName, Q, R, tsModel, ...
                                     dt, known)
        if strcmp(uName, 'tsBased')
            [u, ~, ~, errorFlag] = tsBasedControl(x, sysName, tsModel, ...
                                                  dt, known, Q, R);
            tsFailed = or(tsFailed, errorFlag); 
        elseif strcmp(uName, 'SDRE')
            u = sdre(x, sysName, Q, R);
        end
        
        if strcmp(sysName, 'motorLink')
            dXdt = sys.rhsMotorLink(x, u);
        elseif strcmp(sysName, 'invPend')
            dXdt = sys.rhsInvPend(x, u);
        elseif strcmp(sysName, 'flex2link')
            dXdt = sys.rhsFlex2link(x, u);
        end
        dXdt(end + 1) = x'*Q*x + u'*R*u;
    end

    % 0.3 event detection, event = there is no solution of SDRE
    function [condition, isTerminal, direction] = nonvalidTS()
        condition = 1 - double(tsFailed);
        isTerminal = 1;
        direction = 0;
    end

    % 1.1 integrate
    timesteps = 0:dt:T;
    rhs = @(x) rhsWithCriterion(x, 'tsBased', sysName, Q, R, tsModel, ...
                                dt, known);
    opt = odeset('Events', @(t, x) nonvalidTS);
    [~, tsX] = ode45(@(t, x) rhs(x(1:end-1)), timesteps, [x0; 0], opt);

    if tsFailed
        disp(['TS-based SDRE does not work. x:' num2str(x0')])
        tsCriterion = -1;
        sdreCriterion = -1;
    else
        rhs = @(x) rhsWithCriterion(x, 'SDRE', sysName, Q, R);
        [~, sdreX] = ode45(@(t, x) rhs(x(1:end-1)), timesteps, [x0; 0]);
        tsCriterion = tsX(end, end);
        disp(['Criterion value of new method: ' num2str(tsCriterion)])
        tsX(:, end) = [];       % delete column with criterion values
        sdreCriterion = sdreX(end, end);
        disp(['Criterion value of classic method: ' num2str(sdreCriterion)])
        % 1.2 process data for plot
        sdreX(:, end) = [];     % delete column with criterion values
        nSteps = length(timesteps);
        for iStep=1:nSteps
            tsX(iStep, :) = wrapper(tsX(iStep, :));
            sdreX(iStep, :) = wrapper(sdreX(iStep, :));
        end
    
        % 2. Calculate u and estimates(f, B) at timesteps
        [uList, fTrue, fPred, Btrue, Bpred] = utils.logger(sysName, tsX, ...
                                                           r, tsModel, dt);
        sdreList = zeros(nSteps, r);
        for iStep=1:nSteps
            sdreList(iStep, :) = sdre(sdreX(iStep, :)', sysName, Q, R);
        end
    
        % 3. Plot u, estimates and trajectories
        figure('Name', 'Controls')
        plot(timesteps, sdreList, timesteps, uList)
        legend('SDRE', 'tsBased', 'FontSize', 24)
        ax = gca;
        ax.FontSize = 18;
    
    %     for k=1:n
    %         figure()
    %         title(['f' num2str(k)])
    %         plot(timesteps, fTrue(:, k), timesteps, fPred(:, k))
    %         legend(['f' num2str(k)], ['f' num2str(k) 'Pred'], 'FontSize', 18)
    % 
    %         figure()
    %         title(['B' num2str(k)])
    %         plot(timesteps, Btrue(:, k), timesteps, Bpred(:, k))
    %         legend(['B' num2str(k)], ['B' num2str(k) 'Pred'], 'FontSize', 18)
    %     end
        figure('Name', 'Trajectories')
        lineSpec = ["r-", "r--"; "g-", "g--"; "b-", "b--"; "k-", "k--"];
        labels = cell(2*n, 1);
        hold on
        for k=1:n
            plot(timesteps, sdreX(:, k), lineSpec(k, 1), ...
                 timesteps, tsX(:, k), lineSpec(k, 2))
            labels{2*k - 1} = ['x' num2str(k) '(SDRE)'];
            labels{2*k} = '';
        end
        legend(labels, 'FontSize', 24)
        hold off
        ax = gca;
        ax.FontSize = 18;
    end
end
