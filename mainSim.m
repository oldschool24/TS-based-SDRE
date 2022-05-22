function [tsCriterion, sdreCriterion] = mainSim(sysRhs, dt, T, tsModel, x0)
    arguments
        sysRhs = @sys.rhsMotorLink;
        dt = 0.01;
        T = 0.5;
        tsModel = readfis('models/motorLink.fis'); 
        x0 = [-pi/2; -3*pi];  % change it sometime
    end

    if isequal(sysRhs, @sys.rhsMotorLink)
        sysName = 'motorLink';
        wrapper = @(x) x;
        n = 2;
        r = 1;
        known = [];
    elseif isequal(sysRhs, @sys.rhsInvPend)
        sysName = 'invPend';
        wrapper = @sys.invPendWrapper;
        n = 4;
        r = 1;
        known = [];
    end
    Q = 10 * eye(n);
    R = 5 * eye(r);

    % 1.1 integrate
    u = @(x) tsBasedControl(x, sysName, tsModel, dt, known);
    rhs = @(x) [sysRhs(x, u(x)); x'*Q*x + u(x)'*R*u(x)];
    timesteps = 0:dt:T;
    [~, tsX] = ode45(@(t, x) rhs(x(1:end-1)), timesteps, [x0; 0]);
    u = @(x) sdre(x, sysName, Q, R);
    rhs = @(x) [sysRhs(x, u(x)); x'*Q*x + u(x)'*R*u(x)];
    [~, sdreX] = ode45(@(t, x) rhs(x(1:end-1)), timesteps, [x0; 0]);

    % 1.2 process data for plot
    tsCriterion = tsX(end, end);
    disp(['Criterion value of new method: ' num2str(tsCriterion)])
    tsX(:, end) = [];       % delete column with criterion values
    sdreCriterion = sdreX(end, end);
    disp(['Criterion value of classic method: ' num2str(sdreCriterion)])
%     sdreX(:, end) = [];     % delete column with criterion values
%     nSteps = length(timesteps);
%     for iStep=1:nSteps
%         tsX(iStep, :) = wrapper(tsX(iStep, :));
%         sdreX(iStep, :) = wrapper(sdreX(iStep, :));
%     end
% 
%     % 2. Calculate u and estimates(f, B) at timesteps
%     n = length(x0);
%     r = length(u(tsX(1, :)'));
%     [uList, fTrue, fPred, Btrue, Bpred] = utils.logger(sysName, tsX, ...
%                                                        nSteps, n, r, ...
%                                                        tsModel, dt);
%     sdreList = zeros(nSteps, r);
%     for iStep=1:nSteps
%         sdreList(iStep, :) = sdre(sdreX(iStep, :)', sysName, Q, R);
%     end
% 
%     % 3. Plot u, estimates and trajectories
%     figure('Name', 'Controls')
%     plot(timesteps, sdreList, timesteps, uList)
%     legend('SDRE', 'tsBased', 'FontSize', 18)
% 
% %     for k=1:n
% %         figure()
% %         title(['f' num2str(k)])
% %         plot(timesteps, fTrue(:, k), timesteps, fPred(:, k))
% %         legend(['f' num2str(k)], ['f' num2str(k) 'Pred'], 'FontSize', 18)
% % 
% %         figure()
% %         title(['B' num2str(k)])
% %         plot(timesteps, Btrue(:, k), timesteps, Bpred(:, k))
% %         legend(['B' num2str(k)], ['B' num2str(k) 'Pred'], 'FontSize', 18)
% %     end
%     figure('Name', 'Trajectories')
%     lineSpec = ["r-", "r--"; "g-", "g--"; "b-", "b--"; "k-", "k--"];
%     labels = cell(2*n, 1);
%     hold on
%     for k=1:n
%         plot(timesteps, sdreX(:, k), lineSpec(k, 1), ...
%              timesteps, tsX(:, k), lineSpec(k, 2))
%         labels{2*k - 1} = ['x' num2str(k) '(SDRE)'];
%         labels{2*k} = '';
%     end
%     legend(labels, 'FontSize', 14)
%     hold off
end

function u = sdre(x, sysName, Q, R, A, B)
    if nargin == 4
        if strcmp(sysName, 'motorLink')
            A = [0, 1; -64*sin(x(1))/x(1), -5];
            B = [0; 400];
        elseif strcmp(sysName, 'invPend')
            M = 0.5;
            m = 0.2;  
            b = 0.1;
            l = 0.3;
            I = 0.006;
            g = 9.8;
            
            denominator = (M + m)*I + (m*l*sin(x(2)))^2 + M*m*l^2;
            A = [0, 0, 1, 0;
                 0, 0, 0, 1;
                 0, g*m^2*l^2*sin(x(2))*cos(x(2))/x(2), -b*(I + m*l^2), ...
                 -m*l*(I + m*l^2)*sin(x(2))*x(4);
                 0, m*l*(M+m)*g*sin(x(2))/x(2), -m*l*b*cos(x(2)), ...
                 -m^2*l^2*sin(x(2))*cos(x(2))*x(4)];
            A(3:4, :) = A(3:4, :) / denominator;
            B = [0; 0; ...
                (I + m * l^2) / denominator; ...
                m*l*cos(x(2)) / denominator];
        end
    end
    P = icare(A, B, Q, R);
    u = -inv(R) * B' * P * x;
end
