function dXdt = rhsInvPend(x, u, M, m, b, l, I)
% x1 = x, x2 = theta
% theta - angle of deviation from the vertical
    arguments
        x (4, 1)
        u (1, 1)
        M = 0.5;
        m = 0.2;
        b = 0.1;
        l = 0.3;
        I = 0.006;
    end
    g = 9.8;
    dXdt = zeros(4, 1);
    dXdt(1) = x(3);
    dXdt(2) = x(4);
    Mat = [M + m, -m*l*cos(x(2)); -m*l*cos(x(2)), I + m * l^2];
    vec = [u - b*x(3) - m*l * x(4)^2 * sin(x(2)); m*g*l*sin(x(2))];
    dXdt(3:4) = Mat \ vec;
end