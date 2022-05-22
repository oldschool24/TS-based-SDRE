function [u, fHat, hatB] = tsBasedControl(x, sysName, tsModel, ...
                                          dt, known)
% this function finds SDRE-control based on tsModel, x - state vector
    if nargin < 5
        known = [];
    end
    if strcmp(sysName, 'invPend') 
        x = sys.invPendWrapper(x);
    end
    nRules = length(tsModel.Rules);
    n = length(x);
    nInputs = length(tsModel.Inputs);
    r = nInputs - n;

    % 1. Calculate waveA, waveB: 
    %    tsModel(x(k), u(k)) = waveA * x(k) + waveB * u(k)
    waveA = zeros(n, n);
    waveB = zeros(n, r);
    [~, ~, ~, ~, ruleFiring] = evalfis(tsModel, [x; zeros(r, 1)]);
    
    if sum(ruleFiring) < 1e-12
        disp(['The rules do not work. x:' num2str(x')])
        u = zeros(r, 1);
        fHat = zeros(n, 1);
        hatB = zeros(n, 1);
        return
    end
%     assert(sum(ruleFiring)>1e-3, ['The rules do not work. x:' num2str(x')])
    
    ruleFiring = ruleFiring / sum(ruleFiring);
    [~, out] = getTunableSettings(tsModel);
    thenParams = getTunableValues(tsModel, out);
    % column <-> component of state vector:
    thenParams = reshape(thenParams, [], n);
    thenParams = utils.removeBiasNules(thenParams, nRules, n, r);
%     theta = zeros(nRules, n, n+r);
%     % theta = reshaped thenParams: rule <-> matrix of parameters
    for iRule=1:nRules
        temp = thenParams(1+(iRule-1)*(n+r) : iRule*(n+r), :)';
        A = temp(:, 1:n);
        B = temp(:, n+1:end);
        waveA = waveA + ruleFiring(iRule) * A;
        waveB = waveB + ruleFiring(iRule) * B;
%         theta(iRule, :, :) = thenParams(1+(iRule-1)*(n+r):iRule*(n+r), :)';
    end

    % 2. Calculate hatA, hatB: estimates of A(x), B(x)
    hatA = 1/dt * (waveA - eye(n));
    hatB = 1/dt * waveB;
    if ~isempty(known)
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
        hatA(known, :) = A(known, :);
        hatB(known, :) = B(known, :);
    end
    fHat = hatA * x;

    % 3. Calculate u = SDRE(hatA, hatB)
    Q = 10 * eye(n);
    R = 5 * eye(r);
    P = icare(hatA, hatB, Q, R);
    u = -inv(R) * hatB' * P * x;
end

%     M = 0.5;
%     m = 0.2;  
%     b = 0.1;
%     l = 0.3;
%     I = 0.006;
%     g = 9.8;
%     denominator = (M + m)*I + (m*l*sin(x(2)))^2 + M*m*l^2;
%     fTrue = [x(3); x(4); ...
%         (-b*(I + m*l^2)*x(3) - m^2*l^3*sin(x(2))*x(4)^2 + ...
%         g*m^2*l^2*sin(x(2))*cos(x(2)) - I*m*l*sin(x(2))*x(4)^2) ...
%         / denominator; ...
%         -m*l*(m*l*sin(x(2))*cos(x(2))*x(4)^2 + b*cos(x(2))*x(3) - ...
%         (M+m)*g*sin(x(2))) / denominator];
%     Btrue = [0; 0; ...
%         (I + m * l^2) / denominator; ...
%         m*l*cos(x(2)) / denominator];
%     fHat - fTrue
%     hatB - Btrue

    % analysis of tsModel coefs 
%     [~, out] = getTunableSettings(tsModel);
%     theta = getTunableValues(tsModel, out);
% uWeights = theta(5:6:360);
% uWeights = reshape(uWeights, 4, []);
% sum(abs(uWeights), 2), max(abs(uWeights), [], 2), mean(abs(uWeights), 2)