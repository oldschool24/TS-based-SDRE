function [clipped, periods] = flex2linkWrapper(x)
    % clip x(1), x(2) to [-pi, pi]
    [clipped, rotations] = sys.clipAngle(x(1:2));
    clipped = [clipped; x(3:end)];
    periods = (2*pi) * [rotations; zeros(6, 1)];
end
