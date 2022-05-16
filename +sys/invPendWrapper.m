function x = invPendWrapper(x)
    % clip angle to [-pi, pi]
    if (x(2) > pi) || (x(2) < -pi)
        x(2) = x(2) - 2*pi*floor(x(2)/(2*pi));
        if abs(x(2)) > pi
            x(2) = x(2) - sign(x(2))*(2*pi);
        end
    end
end
