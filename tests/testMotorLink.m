addpath('../')
thetaRange = -pi/2:pi/9:pi/2;
thetadotRange = -pi/2:pi/9:pi/2;
Q = 5*eye(2);
R = 10;

nTests = length(thetaRange) * length(thetadotRange);
criterion = zeros(nTests, 2);  % 1st column - ts-based SDRE, 2nd - SDRE
k = 1;

warning('off', 'fuzzy:general:warnEvalfis_NoRuleFired')
warning('off', 'fuzzy:general:diagEvalfis_OutOfRangeInput')
for theta=thetaRange
    for thetadot=thetadotRange
        tic
        [tsCriterion, sdreCriterion] = mainSim(@sys.rhsMotorLink, ...
            0.01, 10, readfis('../models/motorLink.fis'), ...
            [theta; thetadot], Q, R);
        toc
        criterion(k, 1) = tsCriterion;
        criterion(k, 2) = sdreCriterion;
        k = k + 1;
        save('../results/testMotorLink(5-10)', 'criterion')
    end
end
% save('../data/testMotorLink(10-5)', 'criterion')
