function [clipped, rotations] = clipAngle(angles)
% clip angle to [-pi, pi]
    clipped = zeros(length(angles), 1);
    rotations = zeros(length(angles), 1);
    for iAngles=1:length(angles)
        angle = angles(iAngles);
        if (angle > pi) || (angle < -pi)
            distance = abs(angle - sign(angle)*pi);
            rotation = sign(angle) * ceil(distance / (2*pi));
            newValue = angle - (2*pi)*rotation;
        else
            newValue = angle;
            rotation = 0;
        end
        clipped(iAngles) = newValue;
        rotations(iAngles) = rotation;
    end
end
