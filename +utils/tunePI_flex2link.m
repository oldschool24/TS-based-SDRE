nPoints = 5;

sysName = 'flex2link';
wrapper = @sys.flex2linkWrapper;
n = 8;
r = 2;
xIdxPI = [1, 2];
uIdxPI = [1, 2];


% xRange = [-pi, -pi, -30, -70, -5*pi, -5*pi, -150, -700; 
%            pi,  pi,  30,  70,  5*pi,  5*pi,  150,  700];
xRange = [-pi, -pi, -30, -70, -5*pi, -5*pi, -8*10^4, -8*10^4; 
    pi,  pi,  30,  70,  5*pi,  5*pi,  8*10^4,  8*10^4];
x0Range = [-pi/2, -pi/2, -pi/2, -pi/2, -pi, -pi, -pi, -pi; 
            pi/2,  pi/2,  pi/2,  pi/2,  pi,  pi,  pi,  pi];
x0Grid = utils.uniformGrid(x0Range, nPoints);

integralRange = [zeros(1, length(xIdxPI));
                 Inf(1, length(xIdxPI))];
xRange = [xRange integralRange];
options_extended = odeset('Events', ...
    @(t, x) utils.xRangeEvent(wrapper(x), xRange));

% 1. настроены по отдельности
% Kp = [-400; -200];
% Ki = [20; 15];

% 2. настроены совместно
Kp = [-220; -120];
Ki = [7; 4];

timesteps = 0:0.1:160;

for iPoint=1:nPoints
    x0 = x0Grid(iPoint, :)';
%     x0 = [2.8; 2.9; -20; 66; 14; -0.5; 90; -500];
%     x0 = [-pi/50, -pi/50, 0.03, 0.07, -0.005*pi, 0.005*pi, -0.15, 0.7]';

    % 1. Integrate with initial condition = x0
    [t, X] = ode45( ...
        @(t, x) utils.rhsWithPI(x, sysName, xIdxPI, uIdxPI, Kp, Ki), ...
        timesteps, [x0; zeros(r, 1)]);
            
    % 2. Data postprocessing
    nSteps = length(t);
    uList = zeros(nSteps, r);
    for iStep=1:nSteps
        for k=1:r
            err = X(iStep, xIdxPI(k));   % reference = 0
            I = X(iStep, n+k);
            uList(iStep, uIdxPI(k)) = Kp(k)*err + Ki(k)*I;
        end
        X(iStep, :) = wrapper(X(iStep, :));
    end
    X(:, end-r+1 : end) = [];

%     figure()
%     plot(t, X(:, 1:2))
%     lgd = legend('x_1', 'x_2');
%     lgd.FontSize = 14;
% 
%     figure()
%     plot(t, uList)
%     lgd = legend('u_1', 'u_2');
%     lgd.FontSize = 14;

%     figure()
%     plot(t, X(:, 7:end))
%     lgd = legend('x_7', 'x_8');
%     lgd.FontSize = 14;
end
