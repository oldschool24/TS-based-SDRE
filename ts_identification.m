function ts_identification(load_data, mode)
% identification of TS fuzzy model based on IO data
    global T learn_step
    T = 1;
    learn_step = 0.01;

    if nargin == 1
        mode = 'current_previous';     
    end
    rhs = @rhs_motor_link;
    x0 = [0; -2*pi];
    method = 'SubtractiveClustering';
    if load_data == 1
        if strcmp(mode, 'current')
            load motor_data_current.mat dataset
        else
            load motor_data_current_previous.mat dataset
        end
    else
        dataset = collect_data(rhs, x0, mode);
        data_name = ['motor_data_' mode];
        save(data_name, 'dataset')
    end

    % 1. identify
    opt = genfisOptions(method);
    ts_model = genfis(dataset(:, 1:end-1), dataset(:, end), opt);

    % 2. plot and compare
    plot_identified(20, ts_model, rhs, method, mode)

    writeFIS(ts_model, ['motor_link_ts_' mode])
end

function dXdt = rhs_motor_link(x, u)
% right hand size of motor link system
    dXdt = zeros(2, 1);
    dXdt(1) = x(2);
    dXdt(2) = -64*sin(x(1)) - 5*x(2) + 400*u;
end

function dataset = collect_data(rhs, x0, mode)
% create simulated data for TS model identification
    global T learn_step u_list
    u_list = [];
    timesteps = 0:learn_step:T;
    [~, X] = ode45(@(t, x) rhs(x, u_rand()), timesteps, x0);
    plot(X(:, 1))
    
    n_steps = length(timesteps);
    if strcmp(mode, 'current')  % x(k+1) ~ f(x(k), u(k))
        dataset = zeros(n_steps-1, 3);
        for k=1:n_steps-1
            dataset(k, 1) = u_list(k);
            dataset(k, 2) = X(k, 1);
            dataset(k, 3) = X(k+1, 1);
        end
    else                        % x(k+1) ~ f(x(k), x(k-1), u(k))
        dataset = zeros(n_steps-2, 4);
        for k=1:n_steps-2
            dataset(k, 1) = u_list(k+1);    % here can be problem for indexing
            dataset(k, 2) = X(k+1, 1);
            dataset(k, 3) = X(k, 1);
            dataset(k, 4) = X(k+2, 1);
        end
    end
end

function u = u_rand()
    u = rand() - 0.5;
    global u_list
    u_list = [u_list; u];
end

function plot_identified(T, ts_model, rhs, method, mode)
    global learn_step
    timesteps = 0:learn_step:T;
    u_test = @(t) 0.02*sin(0.1*pi*t) + 0.15*sin(pi*t) + ...
                  0.2*sin(10*pi*t) + 0.2*sin(100*pi*t);
    [~, X_true] = ode45(@(t, x) rhs(x, u_test(t)), timesteps, [0; 0]);
    
    
    X_pred = zeros(length(timesteps), 1);     % X_pred(1:2) = x0(1);
    if strcmp(mode, 'current')
        for k=2:length(timesteps)
            X_pred(k) = evalfis(ts_model, [u_test(timesteps(k-1)); ...
                                           X_true(k-1)]);
        end
    else
        for k=3:length(timesteps)
            X_pred(k) = evalfis(ts_model, [u_test(timesteps(k-1)); ...
                                           X_true(k-1); ...
                                           X_true(k-2)]);
        end
    end
    plot(timesteps, X_true(:, 1), timesteps, X_pred)
    legend('true', 'identified')
    title(method)
end

function u = ts_based_control(ts_model, x)
    [~, ~, ~, ~, ruleFiring] = evalfis(ts_model, x);
    ruleFiring = ruleFiring / sum(ruleFiring);
    [~, out] = getTunableSettings(ts_model);
    

    u = 0;
end
