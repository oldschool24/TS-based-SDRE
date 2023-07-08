config = jsondecode(fileread(configPath));

expName = config.expName;
isLoad = config.isLoad;
sysName = config.sysName;
dt = config.dt;
isWrap = config.isWrap;
nPeriod = config.nPeriod;
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
method = config.method;
x0 = config.x0;
r = config.r;
n = config.n;

if isWrap && nPeriod>0
    error('for wrapped version nPeriod should be equal 0.')
end

x0Range = utils.reduceRange(xRange, reduction);
