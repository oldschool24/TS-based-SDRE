function testFlex2link(modelPath, q, r)
    [dt, T, nPoints, reduction] = extractParams(modelPath);

    addpath('../')
    q1Range = -pi:pi/2:pi;
    q2Range = -pi:pi/2:pi;
    z1Range = [-350, 350];
    z2Range = [-250, 250];
    
    Q = q * eye(8);
    R = r * eye(2);
    nTests = length(q1Range) * length(q2Range) * length(z1Range) * length(z2Range);
    criterion = zeros(nTests, 8+2);
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
                        modelPath, 'flex2link', dt, 15, ...
                        x0, Q, R, @ode23s);
                    toc
                    criterion(k, 1:8) = x0';
                    criterion(k, 9) = tsCriterion;
                    criterion(k, 10) = sdreCriterion;
                    k = k + 1;
                end
            end
        end
    end
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
