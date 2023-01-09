T = 1;
method = 'SubtractiveClustering';

% u = [u_1 ... u_r]', x = [x_1 ... x_n]'
if strcmp(sysName, 'motorLink')
    r = 1;      
    n = 2;
    uRange = [-0.5;
               0.5];
    xRange = [-pi, -6*pi;
               pi,  6*pi];
    x0 = [0; 0];
    isNormalize = false;
    withPI = false;
    x0Range = [-pi/2 -2*pi; ...
                pi/2  2*pi];
    nPoints = 4;
elseif strcmp(sysName, 'invPend')
    r = 1;
    n = 4;
    uRange = [-3;
               3];
%         uRange = [-97 138];
    xRange = [-6, -pi/2, -12, -10; 
               6,  pi/2,  12,  10];
    x0 = [0.1; pi/30; -0.05; -pi/10];
    isNormalize = false;
    withPI = false;
    x0Range = zeros(size(xRange));
    xAmp = xRange(2, :) - xRange(1, :);
    x0Range(1, :) = xRange(1, :) + 0.025 * xAmp;    % TODO: retry it for motorLink
    x0Range(2, :) = xRange(2, :) - 0.025 * xAmp;
    nPoints = 200;
elseif strcmp(sysName, 'flex2link')
    r = 2;
    n = 8;
%         uRange = [-211 -615; 
%                    254  608];
    uRange = [-1 -3.7;  % 100 times less
               1  3.7];
%     xRange = [-pi, -pi, -30, -70, -5*pi, -5*pi, -150, -700;  % real
%                pi,  pi,  30,  70,  5*pi,  5*pi,  150,  700];
    xRange = [-pi, -pi, -700, -500, -7*pi, -9*pi, -4*10^4, -4*10^4; 
               pi,  pi,  700,  500,  7*pi,  9*pi,  4*10^4,  4*10^4];
    x0 = zeros(n, 1);
    isNormalize = false;


    x0Range = [-pi/2, -pi/2, -pi/2, -pi/2, -pi, -pi, -pi, -pi; 
                pi/2,  pi/2,  pi/2,  pi/2,  pi,  pi,  pi,  pi];
    nPoints = 300;

    withPI = true;
    xIdxPI = [1, 2];
    uIdxPI = [1, 2];
    Kp = [-220; -120];
    Ki = [7; 4];
end
