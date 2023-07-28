function simStats = mainSim(modelPath, sysName, dt, T, x0, Q, R, stopType, ...
                            xRange, ode, isWrap, imgDir, known, ...
                            isAnalyze, isWarn, verbose)
    arguments
        modelPath
        sysName
        dt double {mustBePositive}
        T = 0.5
        x0 = [-pi/2; -3*pi]  % change it sometime
        Q = 10 * eye(2)
        R = 5
        stopType = 'trajectory'
        xRange = [-inf, -10*pi; inf, 10*pi]
        ode = @ode45         % flex2link: ode15s faster
        isWrap = false
        imgDir = ''
        known = []
        isAnalyze = false
        isWarn = false
        verbose = true
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
    tsOpt = odeset('RelTol', 1e-3, 'AbsTol', 1e-6); 
    sdreOpt = odeset('RelTol', 1e-3, 'AbsTol', 1e-6);
    if strcmp(sysName, 'motorLink')
        wrapper = @(x) x;
        n = 2;
        r = 1;
    elseif strcmp(sysName, 'invPend')
        wrapper = @sys.invPendWrapper;
        n = 4;
        r = 1;
    elseif strcmp(sysName, 'flex2link')
        wrapper = @sys.flex2linkWrapper;
        n = 8;
        r = 2;
        tsOpt = odeset('RelTol', 5e-3, 'AbsTol', 5e-6); % ode23s (long)
        sdreOpt = odeset('RelTol', 5e-3, 'AbsTol', 5e-6);
    end
    
    % 0.2 add calculation of criterion to rhs
    function dXdt = rhsWithCriterion(x, uName, sysName, Q, R, ...
                                     extendedModel, dt, known)
        if strcmp(uName, 'tsBased')
            u = tsBasedControl(x, extendedModel, sysName, dt, known, ...
                               Q, R, isWrap);
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
    function [condition, isTerminal, direction] = nonvalidTS(x)
        [u, ~, ~, errorFlag] = tsBasedControl( ...
            x, extendedModel, sysName, dt, known, Q, R, isWrap);
        if strcmp(stopType, 'trajectory')
            [condition, isTerminal, direction] = stopDueTraj(x, xRange);
        elseif strcmp(stopType, 'control')
            [condition, isTerminal, direction] = stopDueControl(u, u0);
        end
        condition = [condition; 1-double(errorFlag)];
        isTerminal = [isTerminal; 1];
        direction = [direction; 0];
    end

    if ~isWarn
        warning('off', 'fuzzy:general:warnEvalfis_NoRuleFired')
        warning('off', 'fuzzy:general:diagEvalfis_OutOfRangeInput')
    end

    % 1.1 integrate
    timesteps = 0:dt:T;
    rhs = @(x) rhsWithCriterion(x, 'SDRE', sysName, Q, R);
    global w_alpha;
    global x_old;
    w_alpha = [];
    x_old = x0;
    
    sdreWallTime = tic;
    [t, sdreX] = ode(@(t, x) rhs(x(1:end-1)), timesteps, [x0; 0], sdreOpt);
    sdreWallTime = toc(sdreWallTime);
    simStats.sdreTime = t(end);
    sdreCriterion = sdreX(end, end);
    sdreX(:, end) = [];   

    if strcmp(stopType, 'control')
        u0 = tsBasedControl(x0, extendedModel, sysName, dt, known, Q, R, isWrap);
    end
    rhs = @(x) rhsWithCriterion(x, 'tsBased', sysName, Q, R, ...
                                extendedModel, dt, known);
    tsOpt = odeset(tsOpt, 'Events', @(t, x) nonvalidTS(x(1:end-1)));
    tsWallTime = tic;
    [t, tsX] = ode(@(t, x) rhs(x(1:end-1)), timesteps, [x0; 0], tsOpt);
    tsWallTime = toc(tsWallTime);
    simStats.stopPoint = tsX(end, 1:end-1);
    simStats.tsTime = t(end);
    tsCriterion = tsX(end, end);
    tsX(:, end) = [];   % delete column with criterion values

    simStats.insideEpsTube = insideEpsTube(tsX, sdreX);
    simStats.tsWallTime = tsWallTime;
    simStats.tsCriterion = tsCriterion;
    simStats.sdreWallTime = sdreWallTime;
    simStats.sdreCriterion = sdreCriterion;

    if verbose
        if simStats.tsTime < simStats.sdreTime
            disp(['TS-based SDRE does not work. x:' num2str(simStats.stopPoint)])
            disp(['t:' num2str(simStats.tsTime)])    
        end
        disp(['Criterion value of new method: ' num2str(tsCriterion)])
        disp(['Criterion value of classic method: ' num2str(sdreCriterion)])
    end

    if or(~isempty(imgDir), isAnalyze)
        % 2.1 Calculate u, estimates(x, f, B) at timesteps
        [uList, f_true, f_pred, B_true, B_pred] = utils.logger( ...
            sysName, tsX, r, extendedModel, dt, known, Q, R, isWrap, false);
        [nSteps, ~] = size(tsX);
        sdreList = zeros(nSteps, r);    % SDRE values at timesteps
        for iStep=1:nSteps
            sdreList(iStep, :) = sdre(sdreX(iStep, :)', sysName, Q, R);
        end

        % 2.2 process trajectories for plot
        timesteps(nSteps+1 : end) = [];
        predX = zeros(nSteps, n);   % tsModel preds (sys under ts-control)
        tsModel = extendedModel.model;
        modelRange = extendedModel.range;
        for iStep=1:nSteps
            if iStep > 1
                predX(iStep, :) = utils.evalProjection( ...
                    tsModel, [tsX(iStep-1, :), uList(iStep-1, :)]', ...
                    modelRange, isWrap, sysName);
            end
        end
        predX(1, :) = tsX(1, :);
    end

    if isAnalyze
        simStats.tsX = tsX;
        simStats.predX = predX;
        simStats.f_true = f_true;
        simStats.f_pred = f_pred; 
        simStats.B_true = B_true; 
        simStats.B_pred = B_pred;
    end

    if ~isempty(imgDir)
        % 2.3 Plot u, estimates and trajectories
        plotComparison('Controls', 'SDRE', 0, timesteps, sdreList, ...
                       uList, imgDir)
        utils.plotEstimates('f', f_true, f_pred, n, timesteps, ...
                            'f', imgDir)
        for k=1:r
            utils.plotEstimates(['B^' num2str(k)], B_true(:, :, k), ...
                                B_pred(:, :, k), n, timesteps, ...
                                ['B_' num2str(k)], imgDir)
        end
        if isWrap
            for iStep=1:nSteps
                tsX(iStep, :) = wrapper(tsX(iStep, :)');
                sdreX(iStep, :) = wrapper(sdreX(iStep, :)');
            end
        end
        for k=1:2:n
            plotComparison(['Trajectories-', num2str(k)], 'SDRE', ...
                           k-1, timesteps, sdreX(1:nSteps, [k, k+1]), ...
                           tsX(:, [k, k+1]), imgDir)
            plotComparison(['Trajectories-Preds-' num2str(k)], 'TS', ...
                           k-1, timesteps, tsX(:, [k, k+1]), ...
                           predX(:, [k, k+1]), imgDir)
        end
    end
end

function [condition, isTerminal, direction] = stopDueTraj(x, xRange)
    [condition, isTerminal, direction] = utils.xRangeEvent(x, xRange);
end

function [condition, isTerminal, direction] = stopDueControl(u, u0)
    if norm(u) > 1.03 * norm(u0)
        condition = 0;
    else
        condition = 1;
    end
    isTerminal = 1;
    direction = 0;
end

function isInside = insideEpsTube(tsX, sdreX)
    eps = max(abs(sdreX(ceil(end/2):end, :)));
    
    % light version
    if all(tsX(end, :) > -eps) && all(tsX(end, :) < eps)
        isInside = true;
    else
        isInside = false;
    end
end

function plotComparison(figName, lineName, kStart, timesteps, ...
                        sdreData, tsData, imgDir)
    figure('Name', figName)
    if strcmp(figName, 'Controls')
        varName = 'u';
    else
        varName = 'x';
    end

    [~, nComponents] = size(sdreData);
    labels = cell(2*nComponents, 1);
    hold on
    for k=1:nComponents
        sdreLine = plot(timesteps, sdreData(:, k), 'LineWidth', 1.5);
        labels{2*k - 1} = [varName num2str(k + kStart) '(' lineName ')'];
        color = get(sdreLine, 'Color');
        plot(timesteps, tsData(:, k), '--', 'Color', color, 'LineWidth', 1.5)
        labels{2*k} = '';
    end
    legend(labels, 'FontSize', 24)
    hold off
    ax = gca;
    ax.FontSize = 18;
    exportgraphics(gca, fullfile(imgDir, [figName '.png']))
end
