function u = tsBasedControl(x, tsModel, learnStep)
% this function finds SDRE-control based on tsModel, x - state vector
    nRules = length(tsModel.Rules);
    n = length(x);
    nInputs = length(tsModel.Inputs);
    m = nInputs - n;

    % 1. Calculate waveA, waveB: 
    %    tsModel(x(k), u(k)) = waveA * x(k) + waveB * u(k)
    waveA = zeros(n, n);
    waveB = zeros(n, m);
    [~, ~, ~, ~, ruleFiring] = evalfis(tsModel, [x; zeros(m, 1)]);
    
    assert(sum(ruleFiring)>1e-3, ['The rules do not work. x:' num2str(x')])
    
    ruleFiring = ruleFiring / sum(ruleFiring);
    [~, out] = getTunableSettings(tsModel);
    thenParams = getTunableValues(tsModel, out);
    % column <-> component of state vector:
    thenParams = reshape(thenParams, [], n);
    thenParams = removeBiasNules(thenParams, nRules, n, m);
%     theta = zeros(nRules, n, n+m);
%     % theta = reshaped thenParams: rule <-> matrix of parameters
    for iRule=1:nRules
        temp = thenParams(1+(iRule-1)*(n+m) : iRule*(n+m), :)';
        A = temp(:, 1:n);
        B = temp(:, n+1:end);
        waveA = waveA + ruleFiring(iRule) * A;
        waveB = waveB + ruleFiring(iRule) * B;
%         theta(iRule, :, :) = thenParams(1+(iRule-1)*(n+m):iRule*(n+m), :)';
    end

    % 2. Calculate hatA, hatB: estimates of A(x), B(x)
    hatA = 1/learnStep * (waveA - eye(n));
    hatB = 1/learnStep * waveB;
    % 2.1 compare with ground truth
    f = [x(2); -64*sin(x(1)) - 5*x(2)]
    hatA * x
    hatB
    B = [0; 400]


    % 3. Calculate u = SDRE(hatA, hatB)
    Q = 10;
    R = 5;
    P = icare(hatA, hatB, Q, R);
    u = -inv(R) * hatB' * P * x;
end

function res = removeBiasNules(extParams, nRules, n, m)
% Remove bias parameters of tsModel
    linesToRemove = (n+m+1) : (n+m+1) : nRules*(n+m+1);
    extParams(linesToRemove, :) = [];
    res = extParams;
end