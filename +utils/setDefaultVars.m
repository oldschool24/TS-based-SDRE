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
    reduction = 0.5;
    T = 1;
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
    reduction = 0.05;
    T = 1;
    nPoints = 200;
elseif strcmp(sysName, 'flex2link')
    r = 2;
    n = 8;
%         uRange = [-211 -615; 
%                    254  608];
%     xRange = [-pi, -pi, -30, -70, -5*pi, -5*pi, -150, -700;  % real
%                pi,  pi,  30,  70,  5*pi,  5*pi,  150,  700];
    xRange = [-20, -20, -700, -500, -7*pi, -9*pi, -4*10^4, -4*10^4; 
               20,  20,  700,  500,  7*pi,  9*pi,  4*10^4,  4*10^4];
    x0 = zeros(n, 1);
    isNormalize = false;

    uRange = [-20 -50;  
               20  50];
%     x0Range = [-pi/2, -pi/2, -pi/2, -pi/2, -pi, -pi, -pi, -pi; 
%                 pi/2,  pi/2,  pi/2,  pi/2,  pi,  pi,  pi,  pi];
    reduction = 0.7;
    T = 1;         % dt = 0.01 => 1
    nPoints = 300;   % dt = 0.01 => 300

    withPI = true;
    xIdxPI = [1, 2];
    uIdxPI = [1, 2];
    Kp = [-220; -120];
    Ki = [7; 4];
    Kd = [-0.6; 0.01]; 
end

x0Range = utils.reduceRange(xRange, reduction);
