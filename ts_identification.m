function ts_identification(load_data, mode)
% identification of TS fuzzy model based on IO data
% load_data: 1, if load dataset; 0 if create
% mode: 'current' x(k+1) ~ f(x(k), u(k)); 
%       'current_previous' x(k+1) ~ f(x(k), x(k-1), u(k))

    T = 1;
    learn_step = 0.01;

    if nargin == 1
        mode = 'current';     
    end
    rhs = @rhs_motor_link;
    m = 1;      % u = [u_1 ... u_m]';
    n = 2;      % x = [x_1 ... x_n]';
    x0 = [0; -2*pi];
    method = 'SubtractiveClustering';
    if load_data == 1
        if strcmp(mode, 'current')
            load motor_data_current.mat dataset
        else
            load motor_data_current_previous.mat dataset
        end
    else
        dataset = collect_data(rhs, x0, mode, T, learn_step);
        data_name = ['motor_data_' mode];
        save(data_name, 'dataset')
    end

    % 1. identify number of rules and antecedents params: x(k+1) ~ f(x(k))
    opt = genfisOptions(method); 
    % !!! attention: indexing 2=m+1, m=dim(u)
    ts_model = genfis(dataset(:, m+1:end-1), dataset(:, end), opt);  

    % 2. identify consequents of extended model: x(k+1) ~ f(x(k), u(k))
    % Note: u(k) is not included in if-condition
    ts_model = addInput(ts_model, [-0.5 0.5], 'Name', 'u');
    % consequent: coefficients for u are zero -> u doesn't affect output  
    plot_identified(ts_model, rhs, method, mode, 20, learn_step)
    then_params = RLS(ts_model, dataset, m, n);
    [~, out] = getTunableSettings(ts_model);
    ts_model = setTunableValues(ts_model, out, then_params');
    % ts_model.consequents = RLS(ts_model, dataset, m, n)

    % 3. plot and compare
    plot_identified(ts_model, rhs, method, mode, 20, learn_step)

    % 4. save
    writeFIS(ts_model, ['motor_link_ts_' mode])
end

function dXdt = rhs_motor_link(x, u)
% right hand size of motor link system
    dXdt = zeros(2, 1);
    dXdt(1) = x(2);
    dXdt(2) = -64*sin(x(1)) - 5*x(2) + 400*u;
end

% function dataset = collect_data(rhs, x0, mode, T, learn_step)
% % create simulated data for TS model identification
%     global u_list
%     u_list = [];
%     timesteps = 0:learn_step:T;
%     [~, X] = ode45(@(t, x) rhs(t, x, u_rand()), timesteps, x0);
%     plot(X(:, 1))
%     
%     n_steps = length(timesteps);
%     if strcmp(mode, 'current')  % x(k+1) ~ f(x(k), u(k))
%         dataset = zeros(n_steps-1, 3);
%         for k=1:n_steps-1
%             dataset(k, 1) = u_list(k);
%             dataset(k, 2) = X(k, 1);
%             dataset(k, 3) = X(k+1, 1);
%         end
%     else                        % x(k+1) ~ f(x(k), x(k-1), u(k))
%         dataset = zeros(n_steps-2, 4);
%         for k=1:n_steps-2
%             dataset(k, 1) = u_list(k+1);    % here can be problem for indexing
%             dataset(k, 2) = X(k+1, 1);
%             dataset(k, 3) = X(k, 1);
%             dataset(k, 4) = X(k+2, 1);
%         end
%     end
% end

function dataset = collect_data(rhs, x0, mode, T, learn_step)
% create simulated data for TS model identification
    global u_list
    timesteps = 0:learn_step:T;
    [~, X] = ode45(@(t, x) rhs(x, u_rand()), timesteps, x0);
    plot(X(:, 1))
    
    n_steps = length(timesteps);
    if strcmp(mode, 'current')  % x(k+1) ~ f(x(k))
        dataset = zeros(n_steps-1, 2);
        for k=1:n_steps-1
            dataset(k, 1) = u_list(k);
            dataset(k, 2) = X(k, 1);
            dataset(k, 3) = X(k+1, 1);
        end
    end
end

function u = u_rand()
    u = rand() - 0.5;
    global u_list
    u_list = [u_list; u];       % here can be problems with indexing
end

function then_params = RLS(ts_model, dataset, m, n)
% this function find consequents parameters of ts_model
% using RLS-algorithm on the dataset
% m = length(u), n = length(x), x - state vector, u - control vector
    [n_samples, ~] = size(dataset);
    n_rules = length(ts_model.Rules);
    % 1. Extract ground truth from dataset
    X = dataset(:, end-n+1:end);
    input = dataset(:, 1:end-n);
    % 2. Calculating the firing strength of each rule
    firings = zeros(n_samples, n_rules);
    for k=1:n_samples
        [~, ~, ~, ~, ruleFiring] = evalfis(ts_model, input(k, :));
        ruleFiring = ruleFiring / sum(ruleFiring);
        firings(k, :) = ruleFiring;
    end
    % 3. Create Phi = [phi(1)'; ...; phi(n_d)']
    firings = repelem(firings, 1, n+m+1);
    input(:, end+1) = 1;  % fake input for bias parameter
    Phi = repmat(input, 1, n_rules) .* firings;
    % 4. RLS: use mldivide
    then_params = Phi \ X;
end

function plot_identified(ts_model, rhs, method, mode, T, learn_step)
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
% this function finds SDRE-control based on ts_model, x - state vector
    n_rules = length(ts_model.Rules);
    n = length(x);
    % 1. Calculate A_wave, B_wave: 
    % ts_model(x(k), u(k)) = A_wave * x(k) + B_wave * u(k)
    [~, ~, ~, ~, ruleFiring] = evalfis(ts_model, x);
    ruleFiring = ruleFiring / sum(ruleFiring);
    [~, out] = getTunableSettings(ts_model);
    then_params = getTunableValues(ts_model, out);
    then_params = reshape(then_params, n_rules, []);
    then_params = 
    A_wave = then_params(:, 1:n);
    B_wave = then_params(:, n+1:end);

    % 2. Calculate A_hat, B_hat: estimates of A(x), B(x)


    % 3. Calculate SDRE 
end
