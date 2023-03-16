function modelRange = getTsRange(tsModel, n)
% determine min and max values at which the model works   
    minFiring = 1e-318 ^ (1/n);
    modelRange = zeros(2, n);
    for k=1:n
        params = [tsModel.Inputs(k).MembershipFunctions.Parameters]; 
        means = params(2:2:end);    % mu in gaussian 
        [~, lowIdx] = min(means);   % left fuzzy set
        [~, highIdx] = max(means);  % right fuzzy set
        
        sigma = params(2*lowIdx-1);
        mu = params(2*lowIdx);
        low = invGaussian(minFiring, mu, sigma, 'left');
        modelRange(1, k) = low;

        sigma = params(2*highIdx-1);
        mu = params(2*highIdx);
        high = invGaussian(minFiring, mu, sigma, 'right');
        modelRange(2, k) = high;
    end
end

function x = invGaussian(y, mu, sigma, rootType)
% find root of equation: y = gaussmf(x) = exp(-(x-mu)^2 / (2*sigma^2))
    if strcmp(rootType, 'left')
        x = mu - sqrt(-2 * sigma^2 * log(y));
    elseif strcmp(rootType, 'right')
        x = mu + sqrt(-2 * sigma^2 * log(y));
    end
end
