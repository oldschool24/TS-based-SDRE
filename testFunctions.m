function res = testFunctions(sysName)
    res = {};
    if strcmp(sysName, 'motorLink')
        r = 1;
        T = 20;
        freq = [50, 5, 0.5, 0.05];
        delay = zeros(size(freq));
        uniformInterval = [-0.2 0.2; -0.2 0.2; 0.1 0.15; -0.02 0.02];
        uniformInterval = reshape(uniformInterval, 1, 4, 2);
%         uniformInterval = [-0.2 0.2; -0.4 0.2; 0.1 0.15; -0.01 0];
    elseif strcmp(sysName, 'invPend')
        r = 1;
        T = 10;
        freq = [50, 5, 0.5, 0.05, 0.5];
        delay = zeros(size(freq));
        uniformInterval = [-2 2; -2 2; 1 1.5; -0.2 0.2; -1.5 -1];
        uniformInterval = reshape(uniformInterval, 1, 5, 2);
    elseif strcmp(sysName, 'flex2link')
        r = 2;
        T = 10;
        freq = [30, 10, 5, 0.5, 5;
                30, 10, 5, 0.5, 5];
%         freq = [30, 5, 0.5, 0.05, 0.5;
%                 30, 5, 0.5, 0.05, 0.5];
        delay = zeros(size(freq));
        uniformInterval = zeros(r, 5, 2);
        uniformInterval(1, :, :) = [-3  3; 
                                    -3  3; 
                                    1.5 2.75; 
                                  -0.25 0.25;
                                  -2.75 -1.5];
        uniformInterval(2, :, :) = [-3  3;
                                    -3  3; 
                                    1.5 2.75; 
                                  -0.25 0.25;
                                  -2.75 -1.5];
%         uniformInterval(1, :, :) = [-5  5; 
%                                     -5  5; 
%                                     2.5 3.75; 
%                                    -0.5 0.5;
%                                   -3.75 -2.5];
%         uniformInterval(2, :, :) = [-5  5;
%                                     -5  5; 
%                                     2.5 3.75; 
%                                    -0.5 0.5;
%                                   -3.75 -2.5];
    end
    amp = (uniformInterval(:, :, 2) - uniformInterval(:, :, 1)) / 2;
    amp = squeeze(amp);
    offset = (uniformInterval(:, :, 1) + uniformInterval(:, :, 2)) / 2;
    offset = squeeze(offset);

    % 1. Simple periodic
    res{end+1} = @(t) amp(:, 1) .* ...
                      cos(2*pi*freq(:, 1)*t + delay(:, 1)) + offset(:, 1);
    % 2. Multisine
    res{end+1} = @(t) multisine(t, amp, freq, delay, offset);
    % 3. Sine sweep
    res{end+1} = @(t) sweepsine(t, T, amp(:, 1), freq, offset(:, 1));
%     % 4. Uniform
%     res{end+1} = @(t) amp(:, 1).*rand(r, 1) + offset(:, 1);
%     % 5. Growing uniform
%     res{end+1} = @(t) growing(t, T, amp(:, 1), offset(:, 1), r);
%     % 6. Spiky
%     res{end+1} = @(t) spiky(t, T, uniformInterval, r);
% 
%     % uncomment if you need plots
%     timesteps = 0:0.01:T;
%     figure()
%     nColumns = 2;
%     nLines = ceil(length(res)/nColumns);
%     for k=1:length(res)
%         subplot(nLines, nColumns, k)
%         y = zeros(length(timesteps), r);
%         for iStep=1:length(timesteps)
%             y(iStep, :) = res{k}(timesteps(iStep));
%         end
%         size(res{k}(timesteps(iStep)))
%         plot(timesteps, y)
%         ax = gca;
%         ax.FontSize = 18;
%     end
end

function value = multisine(t, amp, freq, delay, offset)
    [r, nSins] = size(amp);
    value = zeros(r, 1);
    for iSin=1:nSins
        value = value ... 
              + amp(:, iSin) .* sin(2*pi*freq(:, iSin)*t + delay(:, iSin));
        value = value + offset(:, iSin);
    end
end

function value = sweepsine(t, T, amp, freq, offset)
    [~, nSins] = size(freq);
    lenPhase = T / (nSins-1);
    if abs(t - 0) < 1e-10
        iPhase = 1;
    else
        iPhase = min(ceil(t / lenPhase), nSins);
    end
    if iPhase < nSins
        diffPhase = (freq(:, iPhase+1) - freq(:, iPhase)) / (2 * lenPhase);
    else
        diffPhase = 0;
    end
    value = amp .* sin(2*pi*(freq(:, iPhase)*t + diffPhase*t^2)) + offset;
end

function value = growing(t, T, amp, offset, r)
    currentAmp = amp/T * t;
    % value must be in [offset - currentAmp/2, offset + currentAmp/2]
    value = offset + currentAmp .* rand(r, 1) - currentAmp/2;
end

function value = spiky(t, T, uniformInterval, r)
% spiky control signal with 3 peaks
% divide [0, T] into 36 equal parts, every 9th corresponds to the peak
    if mod(ceil(t/T * 36), 9) == 0 
        phase_type = 3; % peak
    else 
        phase_type = 4; % plateau
    end
    left = uniformInterval(:, phase_type, 1);
    right = uniformInterval(:, phase_type, 2);
    value = left + (right - left) .* rand(r, 1);
end
