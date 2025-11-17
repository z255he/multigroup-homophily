function H = homophily_matrix(h)
    n = length(h);
    H = zeros(n); % Initialize n×n matrix

    for i = 1:n
        H(i, :) = (1 - h(i)); % Fill row i with (1 - h_i)
        H(i, i) = 1 + (n - 1) * h(i); % Set diagonal element
    end
    H = 1/n * H;
end

% Final time and uniform time vector
tf = 2000;
t_uniform = linspace(0, tf, 10000)';

% Initial conditions
y0 = [0.01; 1e-6; 0.8; 0.69; 0.95; 0.2];

% Fixed parameters
param_base.h = [0.5, 0.8, 0.8];
param_base.pop = [8e4, 1e4, 1e4];
%options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);  % standard
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);  % grid/heatmap

% Parameter range for rho
rho1_values = 0.5:0.01:0.9;
rho2_values = 0.1:0.02:0.9;
num_sims = length(rho1_values)*length(rho2_values);

% Preallocate result storage (cell array in case y1 has multiple columns)
results = cell(num_sims, 1);

for i = 1:length(rho1_values)
    for j = 1:length(rho2_values)
        % Update parameter
        param = param_base;
        param.pop(1) = 1e5*rho1_values(i);
        param.pop(2) = (1e5 - param.pop(1))*rho2_values(j);
        param.pop(3) = 1e5 - param.pop(1) - param.pop(2);
        param.H = homophily_matrix(param.h);  % Assuming this returns a 3x3 matrix

        % Solve the ODE
        [t_sol, y_sol] = ode45(@(t, y) ode3(t, y, param), [0, tf], y0, options);

        % Interpolate to uniform time grid
        y_interp = interp1(t_sol, y_sol, t_uniform, 'linear');

        % Store result
        idx = (i - 1)*length(rho2_values) + j;
        results{idx} = struct('rho1', rho1_values(i), 'rho2', rho2_values(j), 't', t_uniform, 'y', y_interp);
    end
end


% Save files
num_vars = size(results{1}.y, 2);         % number of variables in y_interp
num_time = length(results{1}.t);          % number of time points
num_sims = length(results);               % number of parameter values

for var_idx = 1:num_vars
    % Initialize matrix: first column is time
    data_matrix = zeros(num_time, num_sims + 1);
    data_matrix(:, 1) = results{1}.t;     % assume same time vector for all
    labels = cell(1, num_sims);

    % Fill in solution data
    for i = 1:num_sims
        data_matrix(:, i+1) = results{i}.y(:, var_idx);
        labels{i} = sprintf('rho1_%.2f_rho2_%.2f', results{i}.rho1, results{i}.rho2);
    end

    % Combine with 'time' as the first header
    header = [{'time'}, labels];

    % Write to CSV
    filename = sprintf('solution_var%d.csv', var_idx);
    fid = fopen(filename, 'w');
    fprintf(fid, '%s,', header{1:end-1});
    fprintf(fid, '%s\n', header{end});
    fclose(fid);
    dlmwrite(filename, data_matrix, '-append');
end

