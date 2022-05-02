function dXdt = rhsMotorLink(x, u)
% right hand size of motor link system
    dXdt = zeros(2, 1);
    dXdt(1) = x(2);
    dXdt(2) = -64*sin(x(1)) - 5*x(2) + 400*u;
end