addpath('../')
q1Range = -pi:pi/2:pi;
q2Range = -pi:pi/2:pi;
z1Range = [-350, 350];
z2Range = [-250, 250];
Q = eye(8);
R = 5*eye(2);

nTests = length(q1Range) * length(q2Range) * length(z1Range) * length(z2Range);
criterion = zeros(nTests, 2);  % 1st column - ts-based SDRE, 2nd - SDRE
k = 1;
warning('off', 'fuzzy:general:warnEvalfis_NoRuleFired')
warning('off', 'fuzzy:general:diagEvalfis_OutOfRangeInput')
for q1=q1Range
    for q2=q2Range
        for z1=z1Range
            for z2=z2Range
                x0 = [q1; q2; z1; z2; 0; 0; 0; 0];
                tic
                [tsCriterion, sdreCriterion] = mainSim( ...
                    '../models/flex2link.mat', 'flex2link', 0.01, 40, ...
                    x0, Q, R, @ode23s);
                toc
                criterion(k, 1) = tsCriterion;
                criterion(k, 2) = sdreCriterion;
                k = k + 1;
            end
        end
    end
end
