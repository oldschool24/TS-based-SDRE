function x0Grid = uniformGrid(ranges, nPoints)
% create grid of initial points from uniform distribution 
% ranges(i, :) -- range of i-th component
% nDivisions -- the number of segments into which the range is divided
    [~, n] = size(ranges);  % n -- length(state vector)
    x0Grid = zeros(nPoints, n);
    for iComponent=1:n
        scale = ranges(2, iComponent) - ranges(1, iComponent);
        bias = ranges(1, iComponent);
        x0Grid(:, iComponent) = scale * rand(nPoints, 1) + bias;
    end
end