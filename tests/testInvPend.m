addpath('../')
dt = 0.01;
xRange = -1:0.5:1;
thetaRange = -pi/2.5:pi/6:pi/2.5;
xdotRange = [-0.5, 0.5];
thetadotRange = [-pi/4, pi/4];
Q = 5*eye(4);
R = 5;

nTests = length(xRange) * length(thetaRange) * length(xdotRange) * length(thetadotRange);
criterion = zeros(nTests, 2);  % 1st column - ts-based SDRE, 2nd - SDRE
k = 1;
warning('off', 'fuzzy:general:warnEvalfis_NoRuleFired')
warning('off', 'fuzzy:general:diagEvalfis_OutOfRangeInput')
for x=xRange
    for theta=thetaRange
        for xdot=xdotRange
            for thetadot=thetadotRange
                tic
                [tsCriterion, sdreCriterion] = mainSim(@sys.rhsInvPend, ...
                    dt, 10, readfis('../models/invPend.fis'), ...
                    [x; theta; xdot; thetadot], Q, R);
                toc
                criterion(k, 1) = tsCriterion;
                criterion(k, 2) = sdreCriterion;
                k = k + 1;
            end
        end
    end
end
save('../results/testInvPend(5-5)', 'criterion')
