function x = invPendWrapper(x)
    % clip x(2) to [-pi, pi]
    x(2) = sys.clipAngle(x(2));
end
