function testFlex2link(modelPath, q, r, testT)
% q, r - parameters of the control criterion
    [dt, T, nPoints, reduction] = extractParams(modelPath);

    addpath('../')
    q1Range = -pi:pi/2:pi; % there is nothing in common between q and q1,q2
    q2Range = -pi:pi/2:pi; % q1, q2 - angles; 
    z1Range = [-350, 350];
    z2Range = [-250, 250];
    [q1, q2, z1, z2] = ndgrid(q1Range, q2Range, z1Range, z2Range);
    nTests = numel(q1);
    
    Q = q * eye(8);
    R = r * eye(2);
    criterion = zeros(nTests, 8+6+8);
    warning('off', 'fuzzy:general:warnEvalfis_NoRuleFired')
    warning('off', 'fuzzy:general:diagEvalfis_OutOfRangeInput')
    tic
    parfor k=1:nTests
        x0 = [q1(k); q2(k); z1(k); z2(k); 0; 0; 0; 0];
        simStats = mainSim(modelPath, 'flex2link', dt, testT, ...
                           x0, Q, R, @ode15s);
        simStats = [x0', simStats.tsCriterion, simStats.sdreCriterion, ...
                    simStats.tsTime, simStats.sdreTime, ...
                    simStats.tsWallTime, simStats.sdreWallTime, ...
                    simStats.stopPoint];
        for iStats=1:22
            criterion(k, iStats) = simStats(iStats);
        end
    end
    toc
    expName = ['../results/testFlex2link(dt-' num2str(dt) ...
               '_T-' num2str(T) ...
               '_N-' num2str(nPoints) ...
               '_reduct-' num2str(reduction) ...
               '_q-' num2str(q) ...
               '_r-' num2str(r) ...
               ').mat'];
    save(expName, 'criterion')
end

function [dt, T, nPoints, reduction] = extractParams(modelPath)
    dtStart = strfind(modelPath, 'dt-') + 3;
    dtEnd = strfind(modelPath, '_T-') - 1;
    dt = str2double(modelPath(dtStart : dtEnd));

    T_Start = dtEnd + 4;
    T_End = strfind(modelPath, '_N-') - 1;
    T = str2num(modelPath(T_Start : T_End));

    nPointsStart = T_End + 4;
    nPointsEnd = strfind(modelPath, '_reduct-') - 1;
    nPoints = str2num(modelPath(nPointsStart : nPointsEnd));
    
    reductStart = nPointsEnd + 9;
    reduction = str2double(modelPath(reductStart : end-1));
end
