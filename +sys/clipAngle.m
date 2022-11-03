function clipped = clipAngle(angle)
% clip angle to [-pi, pi]
    if (angle > pi) || (angle < -pi)
        clipped = angle - 2*pi*floor(angle/(2*pi));
        if abs(angle) > pi
            clipped = angle - sign(angle)*(2*pi);
        end
    else
        clipped = angle;
    end
end
