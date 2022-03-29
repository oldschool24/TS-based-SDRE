function motor_link_stab()
    

end


function u = ts_based_control(ts_model, x)
    [~, ~, ~, ~, ruleFiring] = evalfis(ts_model, x);
    [~, out] = getTunableSettings(ts_model);

    u = 0;
end

function dXdt = rhs_motor_link(x, u)
% right hand size of motor link system
    dXdt = zeros(2, 1);
    dXdt(1) = x(2);
    dXdt(2) = -64*sin(x(1)) - 5*x(2) + 400*u;
end
