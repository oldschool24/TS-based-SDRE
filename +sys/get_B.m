function B = get_B(x, sysName)
    if strcmp(sysName, 'motorLink')
        B = [0; 400];
    elseif strcmp(sysName, 'invPend')
        M = 0.5;
        m = 0.2;  
        l = 0.3;
        I = 0.006;        
        
        denominator = (M + m)*I + (m*l*sin(x(2)))^2 + M*m*l^2;
        B = [0; 0; ...
             (I + m * l^2) / denominator; m*l*cos(x(2)) / denominator];
    elseif strcmp(sysName, 'flex2link')
        B = zeros(8, 2);
        
        % case k=10^4
        B(7, 1) = -100000;
        B(8, 2) = -100000;
%         % case k=10^3
%         B(7, 1) = B(7, 1) / 10;
%         B(8, 2) = B(8, 2) / 10;       
        % case k=10^2
        B(7, 1) = B(7, 1) / 100;
        B(8, 2) = B(8, 2) / 100; 
    end
end