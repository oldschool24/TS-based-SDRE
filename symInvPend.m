M = sym('M');
m = sym('m');
l = sym('l');
x = sym('x', [4, 1]);
I = sym('I');
g = sym('g');
u = sym('u');
b = sym('b');

Mat = [M + m, -m*l*cos(x(2)); -m*l*cos(x(2)), I + m * l^2];
vec = [u - b*x(3) - m*l * x(4)^2 * sin(x(2)); m*g*l*sin(x(2))];

dXdt = Mat \ vec;
collect(dXdt, u)
