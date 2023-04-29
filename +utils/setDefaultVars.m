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
    config = jsondecode(fileread(configPath));
    nPoints = config.nPoints;
    T = config.T;
    uRange = config.uRange;
    reduction = config.reduction;
    xRange = config.xRange;
    withPI = config.withPI;
    xIdxPI = config.xIdxPI;
    uIdxPI = config.uIdxPI;
    Kp = config.Kp;
    Ki = config.Ki;
    Kd = config.Kd;
    isNormalize = config.isNormalize;
    r = config.r;
    n = config.n;

    x0 = zeros(n, 1);
end

x0Range = utils.reduceRange(xRange, reduction);
