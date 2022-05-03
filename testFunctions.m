function res = testFunctions(T, freq, delay, uniformInterval)
    res = {};
    % here can be problem with ranges like [0.2, 1]
    if nargin == 0
        T = 10;
        freq = [1, 20, 5];
        delay = [pi/6, pi/4, pi/3];
        uniformInterval = [0.3, 0.4; -0.1, 0.1];
    end

    amp = uniformInterval(:, 2) - uniformInterval(:, 1);
    offset = uniformInterval(:, 1);
    % 1. Simple periodic
    res{1} = @(t) amp(1) * cos(2*pi*freq(1)*t + delay(1)) + offset(1);
    % 2. Multisine
    res{2} = @(t) multisine(t, amp, freq(1), delay, offset);
    % 3. Sine sweep
    res{3} = @(t) sweepsine(t, T, amp(1), freq, offset(1));
    % 4. Uniform
    res{4} = @(t) amp(1)*rand() + offset(1);
    % 5. Growing uniform
    res{5} = @(t) growing(t, T, amp(1), offset(1));
    % 6. Spiky
    res{6} = @(t) spiky(t, T, uniformInterval);

    timesteps = 0:0.01:T;
    for k=1:6
        y = zeros(length(timesteps), 1);
        for iStep=1:length(timesteps)
            y(iStep) = res{k}(timesteps(iStep));
        end
        figure()
        plot(timesteps, y)
    end
end

function value = multisine(t, amp, freq, delay, offset)
    value = 0;
    nSin = length(amp);
    for iSin=1:nSin
        value = value + amp(iSin) * sin(2*pi*iSin*freq*t + delay(iSin));
        value = value + offset(iSin);
    end
end

function value = sweepsine(t, T, amp, freq, offset)
    nSin = length(freq);
    lenPhase = T / (nSin-1);
    if abs(t - 0) < 1e-10
        iPhase = 1;
    else
        iPhase = ceil(t / lenPhase);
    end
    diffPhase = (freq(iPhase+1) - freq(iPhase)) / (2 * lenPhase);
    value = amp * sin(2*pi*(freq(iPhase)*t + diffPhase*t^2)) + offset;
end

function value = growing(t, T, amp, offset)
    currentAmp = amp/T * t;
    value = currentAmp * rand() + offset;
end

function value = spiky(t, T, uniformInterval)
% spiky control signal with 3 peaks
% divide [0, T] into 12 equal parts, every 4th corresponds to the peak
    if mod(ceil(t/T * 12), 4) == 0 
        phase_type = 1; % peak
    else 
        phase_type = 2; 
    end
    left = uniformInterval(phase_type, 1);
    right = uniformInterval(phase_type, 2);
    value = left + (right - left) * rand();
end
