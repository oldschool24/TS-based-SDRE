function [tsCriterion, sdreCriterion] = mainSim(modelPath, sysName, ...
                                                dt, T, x0, Q, R, ode)
    arguments
        modelPath = 'motorLink.fis';
        sysName = 'motorLink';
        dt = 0.01;
        T = 0.5;
        x0 = [-pi/2; -3*pi];  % change it sometime
        Q = 10 * eye(2);
        R = 5;
        ode = @ode45;  % flex2link: ode23s faster
    end

    % 0.1 Set default values
    if contains(modelPath, 'fis')
        tsModel = readfis(modelPath);
        extendedModel.model = tsModel;
        extendedModel.normC = [];
        extendedModel.normS = [];
    else
        load(modelPath, "extendedModel")
    end
    if strcmp(sysName, 'motorLink')
        wrapper = @(x) x;
        n = 2;
        r = 1;
        known = [];
        tsOpt = odeset('RelTol', 1e-3, 'AbsTol', 1e-6); 
        sdreOpt = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);
    elseif strcmp(sysName, 'invPend')
        wrapper = @sys.invPendWrapper;
        n = 4;
        r = 1;
        known = [];
        tsOpt = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);
        sdreOpt = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);
    elseif strcmp(sysName, 'flex2link')
        wrapper = @sys.flex2linkWrapper;
        n = 8;
        r = 2;
        known = [];
        tsOpt = odeset('RelTol', 1e-1, 'AbsTol', 5e-3); % default: 1e-3, 1e-6
        sdreOpt = odeset('RelTol', 1e-1, 'AbsTol', 5e-3);
    end
    tsFailed = false;  % var for event detection

%     tsTime = [];
%     sdreTime = [];
%     intTime = [];
    
    % 0.2 add calculation of criterion to rhs
    function dXdt = rhsWithCriterion(x, uName, sysName, Q, R, ...
                                     extendedModel, dt, known)
        if strcmp(uName, 'tsBased')
%             tStart = tic;
            [u, ~, ~, errorFlag] = tsBasedControl(x, extendedModel, ...
                                                  sysName, dt, known, Q, R);
%             tEnd = toc(tStart);
%             tsTime = [tsTime; tEnd];
            tsFailed = or(tsFailed, errorFlag); 
        elseif strcmp(uName, 'SDRE')
%             tStart = tic;
            u = sdre(x, sysName, Q, R);
%             tEnd = toc(tStart);
%             sdreTime = [sdreTime; tEnd];
        end
        
%         tStart = tic;
        if strcmp(sysName, 'motorLink')
            dXdt = sys.rhsMotorLink(x, u);
        elseif strcmp(sysName, 'invPend')
            dXdt = sys.rhsInvPend(x, u);
        elseif strcmp(sysName, 'flex2link')
            dXdt = sys.rhsFlex2link(x, u);
        end
        dXdt(end + 1) = x'*Q*x + u'*R*u;
%         tEnd = toc(tStart);
%         intTime = [intTime; tEnd];
    end

    % 0.3 event detection, event = there is no solution of SDRE
    function [condition, isTerminal, direction] = nonvalidTS()
        condition = 1 - double(tsFailed);
        isTerminal = 1;
        direction = 0;
    end

    % 1.1 integrate
    timesteps = 0:dt:T;
    rhs = @(x) rhsWithCriterion(x, 'tsBased', sysName, Q, R, ...
                                extendedModel, dt, known);
    tsOpt = odeset(tsOpt, 'Events', @(t, x) nonvalidTS);
    [t, tsX] = ode(@(t, x) rhs(x(1:end-1)), timesteps, [x0; 0], tsOpt);

    if tsFailed
        disp(['TS-based SDRE does not work. x:' num2str(tsX(end, :))])
        disp(['t:' num2str(t(end))])
        tsCriterion = -1;
        sdreCriterion = -1;
    else
%         dispTime('TS-based control', tsTime)
        tsCriterion = tsX(end, end);
        disp(['Criterion value of new method: ' num2str(tsCriterion)])
        tsX(:, end) = [];       % delete column with criterion values
%         dispTime('rhs evaluation', intTime)
        
        rhs = @(x) rhsWithCriterion(x, 'SDRE', sysName, Q, R);
        [~, sdreX] = ode(@(t, x) rhs(x(1:end-1)), timesteps, [x0; 0], sdreOpt);
%         dispTime('SDRE control', sdreTime)
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
        [uList, fTrue, fPred, Btrue, Bpred] = utils.logger( ...
            sysName, tsX, r, extendedModel, dt);
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
    
        for k=1:n
            figure()
            title(['f' num2str(k)])
            plot(timesteps, fTrue(:, k), timesteps, fPred(:, k))
            legend(['f' num2str(k)], ['f' num2str(k) 'Pred'], 'FontSize', 18)
    
            figure()
            title(['B' num2str(k)])
            plot(timesteps, Btrue(:, k), timesteps, Bpred(:, k))
            legend(['B' num2str(k)], ['B' num2str(k) 'Pred'], 'FontSize', 18)
        end
        figure('Name', 'Trajectories')
        % TODO: AUTOMATIC COLOURS!
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

function dispTime(name, time)
    minT = round(min(time), 4);
    maxT = round(max(time), 4);
    meanT = round(mean(time), 4);

    output = ['Time of ' name ':'];
    output = [output ' min=' num2str(minT)];
    output = [output ', max=' num2str(maxT)];
    output = [output ', mean=' num2str(meanT)];
    disp(output)
end
