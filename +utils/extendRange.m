function extended = extendRange(range, percent)
    extended = zeros(size(range));
    amp = range(2, :) - range(1, :);
    extended(1, :) = range(1, :) + percent/2 * amp;
    extended(2, :) = range(2, :) - percent/2 * amp;
end