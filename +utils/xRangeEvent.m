function [value, isterminal, direction] = xRangeEvent(x, xRange)
    value = [x - xRange(1, :)'; xRange(2, :)' - x];
    n = length(x);
    isterminal = ones(2*n, 1);
    direction = [-1*ones(n, 1); -1*ones(n, 1)];
end