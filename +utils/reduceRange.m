function reduced = reduceRange(range, percent)
    reduced = range;
    amp = range(2, :) - range(1, :);
    finite = ~isinf(range(1, :));
    reduced(1, finite) = range(1, finite) + percent/2 * amp(finite);
    finite = ~isinf(range(2, :));
    reduced(2, finite) = range(2, finite) - percent/2 * amp(finite);
end