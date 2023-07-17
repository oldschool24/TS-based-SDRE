function A = get_A(x, sysName)
    if strcmp(sysName, 'motorLink')
        if abs(x(1)) < 1e-20
            A = [0, 1; -64, -5];
        else
            A = [0, 1; -64*sin(x(1))/x(1), -5];
        end
    elseif strcmp(sysName, 'invPend')
        M = 0.5;
        m = 0.2;  
        b = 0.1;
        l = 0.3;
        I = 0.006;
        g = 9.8;  
        
        denominator = (M + m)*I + (m*l*sin(x(2)))^2 + M*m*l^2;
        if abs(x(2)) < 1e-20    
            A = [0, 0, 1, 0;
                0, 0, 0, 1;
                0, g*m^2*l^2*cos(x(2)), ...
                -b*(I + m*l^2), -m*l*(I + m*l^2)*sin(x(2))*x(4);
                0, m*l*(M+m)*g, -m*l*b*cos(x(2)), ...
                -m^2*l^2*sin(x(2))*cos(x(2))*x(4)];
        else
            A = [0, 0, 1, 0;
                0, 0, 0, 1;
                0, g*m^2*l^2*sin(x(2))*cos(x(2))/x(2), ...
                -b*(I + m*l^2), -m*l*(I + m*l^2)*sin(x(2))*x(4);
                0, m*l*(M+m)*g*sin(x(2))/x(2), -m*l*b*cos(x(2)), ...
                -m^2*l^2*sin(x(2))*cos(x(2))*x(4)];
        end
        A(3:4, :) = A(3:4, :) / denominator;
    elseif strcmp(sysName, 'flex2link')
        A = [0 0 0 0 1 0 0 0;
             0 0 0 0 0 1 0 0;
             0 0 0 0 0 0 1 0;
             0 0 0 0 0 0 0 1;
             zeros(4, 8)]; 
        if abs(x(1)) < 1e-20
            A(5, 1) = -4445280 + 77518*cos(x(2)) + 559237*cos(x(2))^2;  
            A(6, 1) = -5000940 + 703297*cos(x(2)) + 559237*cos(x(2))^2;
            A(7, 1) = -4445280 + 77518*cos(x(2)) + 559237*cos(x(2))^2;
            A(8, 1) = -5000940 + 703297*cos(x(2)) + 559237*cos(x(2))^2;
        else
            A(5, 1) = (-4445280*sin(x(1)) + 77518*cos(x(2))*sin(x(1)) ...
                + 559237*cos(x(2))^2*sin(x(1))) / x(1);  
            A(6, 1) = (-5000940*sin(x(1)) + 703297*cos(x(2))*sin(x(1)) ...
                + 559237*cos(x(2))^2*sin(x(1))) / x(1);
            A(7, 1) = (-4445280*sin(x(1)) + 77518*cos(x(2))*sin(x(1)) ...
                + 559237*cos(x(2))^2*sin(x(1)))/x(1);
            A(8, 1) = (-5000940*sin(x(1)) + 703297*cos(x(2))*sin(x(1)) ...
                + 559237*cos(x(2))^2*sin(x(1))) / x(1);
        end
        if abs(x(2)) < 1e-20
            A(5, 2) = 77518*cos(x(1)) + 559237*cos(x(1))*cos(x(2));
            A(6, 2) = 4711987*cos(x(1)) + 559237*cos(x(1))*cos(x(2));
            A(7, 2) = 77518*cos(x(1)) + 559237*cos(x(1))*cos(x(2));
            A(8, 2) = 4711987*cos(x(1)) + 559237*cos(x(1))*cos(x(2)) ;
        else
            A(5, 2) = (77518*cos(x(1))*sin(x(2)) ...
                + 559237*cos(x(1))*cos(x(2))*sin(x(2))) / x(2);
            A(6, 2) = (4711987*cos(x(1))*sin(x(2)) ...
                + 559237*cos(x(1))*cos(x(2))*sin(x(2))) / x(2);
            A(7, 2) = (77518*cos(x(1))*sin(x(2)) ... 
                + 559237*cos(x(1))*cos(x(2))*sin(x(2))) / x(2);
            A(8, 2) = (4711987*cos(x(1))*sin(x(2)) ...
                + 559237*cos(x(1))*cos(x(2))*sin(x(2))) / x(2);
        end
           
        A(5, 3) = -56000;
        A(5, 4) = 63000 + 50500*cos(x(2));
        A(5, 5) = 63630*x(5)*sin(x(2)) + 113120*x(6)*sin(x(2)) ...
                + 51005*x(5)*cos(x(2))*sin(x(2));
        A(5, 6) = 56560*x(6)*sin(x(2));
        A(5, :) = -A(5, :) ./ (5*(2828*cos(x(2)) + 10201*cos(x(2))^2 - 93548));

        A(6, 3) = -63000 - 50500*cos(x(2));
        A(6, 4) = 488500 + 101000*cos(x(2));
        A(6, 5) = 493385*x(5)*sin(x(2)) + 127260*x(6)*sin(x(2)) ...
                + 102010*x(5)*cos(x(2))*sin(x(2));
        A(6, 6) = 63630*x(6)*sin(x(2)) + 51005*x(6)*cos(x(2))*sin(x(2)) ...
                + 102010*x(5)*cos(x(2))*sin(x(2));
        A(6, :) = A(6, :) / (5*(2828*cos(x(2)) + 10201*cos(x(2))^2 - 93548));  
        
        A(7, 3) = -4733400 + 510050*cos(x(2))^2 + 141400*cos(x(2));
        A(7, 4) = 63000 + 50500*cos(x(2));
        A(7, 5) = 63630*x(5)*sin(x(2)) + 113120*x(6)*sin(x(2)) ...
                + 51005*x(5)*cos(x(2))*sin(x(2));
        A(7, 6) = 56560*x(6)*sin(x(2));

        A(8, 3) = -63000 - 50500*cos(x(2));
        A(8, 4) = 5165900 - 510050*cos(x(2))^2 - 40400*cos(x(2));
        A(8, 5) = 493385*x(5)*sin(x(2)) + 127260*x(6)*sin(x(2)) ...
                + 102010*x(5)*cos(x(2))*sin(x(2));
        A(8, 6) = 63630*x(6)*sin(x(2)) + 51005*x(6)*cos(x(2))*sin(x(2)) ...
                + 102010*x(5)*cos(x(2))*sin(x(2));
        
        % case k=10^4
        A(7, :) = -2000 * A(7, :) ./ (2828*cos(x(2)) + 10201*cos(x(2))^2 - 93548);
        A(8, :) = 2000 * A(8, :) ./ (2828*cos(x(2)) + 10201*cos(x(2))^2 - 93548);
%         % case k=10^3
%         A(7, :) = A(7, :) / 10;
%         A(8, :) = A(8, :) / 10;
        % case k=10^2
        A(7, :) = A(7, :) / 100;
        A(8, :) = A(8, :) / 100;
    end
end