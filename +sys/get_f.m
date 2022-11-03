function f = get_f(x, sysName)
    if strcmp(sysName, 'motorLink')
        f = [x(2), -64*sin(x(1)) - 5*x(2)];
    elseif strcmp(sysName, 'invPend')
        M = 0.5;
        m = 0.2;  
        b = 0.1;
        l = 0.3;
        I = 0.006;
        g = 9.8;
        
        denominator = (M + m)*I + (m*l*sin(x(2)))^2 + M*m*l^2;
        f = [x(3); x(4); ...
            (-b*(I + m*l^2)*x(3) - m^2*l^3*sin(x(2))*x(4)^2 + ...
            g*m^2*l^2*sin(x(2))*cos(x(2)) - I*m*l*sin(x(2))*x(4)^2) ...
            / denominator; ...
            -m*l*(m*l*sin(x(2))*cos(x(2))*x(4)^2 + b*cos(x(2))*x(3) - ...
            (M+m)*g*sin(x(2))) / denominator];
    elseif strcmp(sysName, 'flex2link')
        f = [x(5);
             x(6);
             x(7);
             x(8);
             zeros(4, 1)];
        denominator = 5*(2828*cos(x(2)) + 10201*cos(x(2))^2 - 93548);
        f(5) = -1/denominator * ( ...
             63000*x(4) - 56000*x(3) - 4445280*sin(x(1)) + 63630*x(5)^2*sin(x(2)) ...
             + 56560*x(6)^2*sin(x(2)) + 77518*cos(x(1))*sin(x(2)) + 77518*cos(x(2))*sin(x(1)) ...
             + 559237*cos(x(2))^2*sin(x(1)) + 50500*x(4)*cos(x(2)) ...
             + 559237*cos(x(1))*cos(x(2))*sin(x(2)) + 113120*x(5)*x(6)*sin(x(2)) ...
             + 51005*x(5)^2*cos(x(2))*sin(x(2)) ...
            );
        f(6) = 1/denominator * ( ...
             488500*x(4) - 63000*x(3) - 5000940*sin(x(1)) + 493385*x(5)^2*sin(x(2)) ...
             + 63630*x(6)^2*sin(x(2)) + 4711987*cos(x(1))*sin(x(2)) ...
             + 703297*cos(x(2))*sin(x(1)) + 559237*cos(x(2))^2*sin(x(1)) ...
             - 50500*x(3)*cos(x(2)) + 101000*x(4)*cos(x(2)) ...
             + 559237*cos(x(1))*cos(x(2))*sin(x(2)) + 127260*x(5)*x(6)*sin(x(2)) ...
             + 102010*x(5)^2*cos(x(2))*sin(x(2)) + 51005*x(6)^2*cos(x(2))*sin(x(2)) ...
             + 102010*x(5)*x(6)*cos(x(2))*sin(x(2)) ...
            );
        f(7) = -10000/denominator * ( ...
             63000*x(4) - 4733400*x(3) - 4445280*sin(x(1)) + 510050*x(3)*cos(x(2))^2 ...
             + 63630*x(5)^2*sin(x(2)) + 56560*x(6)^2*sin(x(2)) ...
             + 77518*cos(x(1))*sin(x(2)) + 77518*cos(x(2))*sin(x(1)) ...
             + 559237*cos(x(2))^2*sin(x(1)) + 141400*x(3)*cos(x(2)) ...
             + 50500*x(4)*cos(x(2)) + 559237*cos(x(1))*cos(x(2))*sin(x(2)) ...
             + 113120*x(5)*x(6)*sin(x(2)) + 51005*x(5)^2*cos(x(2))*sin(x(2)) ...
            );
        f(8) = 10000/denominator * ( ...
             - 63000*x(3) + 5165900*x(4) - 5000940*sin(x(1)) ...
             - 510050*x(4)*cos(x(2))^2 + 493385*x(5)^2*sin(x(2)) ...
             + 63630*x(6)^2*sin(x(2)) + 4711987*cos(x(1))*sin(x(2)) ...
             + 703297*cos(x(2))*sin(x(1)) + 559237*cos(x(2))^2*sin(x(1)) ...
             - 50500*x(3)*cos(x(2)) - 40400*x(4)*cos(x(2)) ...
             + 559237*cos(x(1))*cos(x(2))*sin(x(2)) + 127260*x(5)*x(6)*sin(x(2)) ...
             + 102010*x(5)^2*cos(x(2))*sin(x(2)) + 51005*x(6)^2*cos(x(2))*sin(x(2)) ...
             + 102010*x(5)*x(6)*cos(x(2))*sin(x(2)) ...
            );
    end
end