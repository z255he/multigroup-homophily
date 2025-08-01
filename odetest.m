function H = homophily_matrix(h)
    n = length(h);
    H = zeros(n); % Initialize n×n matrix

    for i = 1:n
        H(i, :) = (1 - h(i)); % Fill row i with (1 - h_i)
        H(i, i) = 1 + (n - 1) * h(i); % Set diagonal element
    end
    H = 1/n * H;
end

% Final time
tf = 100;

% Initial conditions
y0 = [0.01; 1e-6; 0.8; 0.67; 0.95; 0.2];

% Parameters
param.h = [0.5, 0.5, 0.85];
param.pop = [8e4, 1e4, 1e4];
param.H = homophily_matrix(param.h);
param.alpha = 0.1;

% Simulation
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
%[t1, y1] = ode45(@(t, y) ode3(t, y, param), [0, tf], y0, options);
[t1, y1] = ode45(@(t, y) ode3tv(t, y, param), [0, tf], y0, options);

S = y1(:, 1);
I = y1(:, 2);
V = y1(:, 3);
X1 = y1(:, 4);
X2 = y1(:, 5);
X3 = y1(:, 6);

% Plot
figure(1);
hold on;
plot(t1, X1, t1, X2, t1, X3, 'LineWidth',2);
figure(2);
hold on;
plot(t1, I, 'LineWidth',2);