function clipped = clipAngle(angle)
% clip angle to [-pi, pi]
    if (angle > pi) || (angle < -pi)
        clipped = angle - 2*pi*floor(angle/(2*pi));
        if abs(clipped) > pi
            clipped = clipped - sign(clipped)*(2*pi);
        end
    else
        clipped = angle;
    end
end
