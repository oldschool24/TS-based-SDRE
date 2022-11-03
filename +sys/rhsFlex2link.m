function dXdt = rhsFlex2link(x, u)
    arguments
        x (8, 1)   % [q_1; q_2; z_1; z_2; dot(q_1); ...; dot(z_2)]
        u (2, 1)
    end
    k = 10^(3);   % default: 10^4
    mu = 1/k;

    M = [9.77 + 2.02*cos(x(2)), 1.26 + 1.01*cos(x(2));
         1.26 + 1.01*cos(x(2)), 1.12];
    invM = inv(M);
    C = [-1.01*x(6)*sin(x(2)), -1.01 * (x(5)+x(6)) * sin(x(2));
         1.01*x(5)*sin(x(2)), 0];
    G = [79.38*sin(x(1)) + 11.074*sin(x(1)+x(2));
         11.074*sin(x(1)+x(2))];
    J = diag([0.1, 0.1]);
    invJ = inv(J);

    Z = [x(3); x(4)];
    qDot = [x(5); x(6)];
    
    dXdt = zeros(8, 1);
    dXdt(1:4) = x(5:8);
    dXdt(5:6) = -invM * (C*qDot + G + Z);
    dXdt(7:8) = -1/mu * ((invM + invJ)*Z + invM*(C*qDot + G) + invJ*u);
end