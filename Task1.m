% Davyn Timbang
% 10/10/25
% ASEN 3802

clc;
clear;
close all;


% Constants
H = 91.09;            % W/m^2 or as given
T0 = 17.04;           % °C
x = 13.585 / 100;    % convert cm to m
alpha = 130/(2810*960);      % m^2/s (Aluminum 7075-T651)
L = 0.149225;             % rod length [m] (adjust if known)

% Define lambda_n and b_n
n_max = 10;
n = 1:n_max;
lambda_n = (2*n - 1) * pi / (2*L);

% Compute b_n using the analytical solution to the integral:
% ∫ x sin(λx) dx = (sin(λL) - λL*cos(λL)) / λ^2
%b_n = - (2*H./L) .* ((sin(lambda_n*L) - lambda_n*L.*cos(lambda_n*L)) ./ (lambda_n.^2));

for n = 1:n_max
b_n(n) = (8*H*L)/(pi^2) * ((-1.^n)/(2*n-1)^2);
end

% Define time values
t_values = [1, 1000]; % seconds

% Preallocate
u_t = zeros(length(t_values), 1);

% Compute temperature for each time
for i = 1:length(t_values)
    t = t_values(i);
    sum_term = sum(b_n .* sin(lambda_n * x) .* exp(-lambda_n.^2 * alpha * t));
    u_t(i) = T0 + H*x + sum_term;
end

% Display results
fprintf('Temperature at x = %.3f m\n', x);
fprintf('At t = 1 s:    u = %.3f °C\n', u_t(1));
fprintf('At t = 1000 s: u = %.3f °C\n', u_t(2));

% Convergence plot: u(x,t) vs n
n_series = 0:10;
u_converge_1 = zeros(size(n_series));
u_converge_1000 = zeros(size(n_series));

for k = 0:10
    if k == 0
        u_converge_1(k+1) = T0 + H*x;
        u_converge_1000(k+1) = T0 + H*x;
    else
        partial_sum_1 = sum(b_n(1:k) .* sin(lambda_n(1:k) * x) .* exp(-lambda_n(1:k).^2 * alpha * 1));
        partial_sum_1000 = sum(b_n(1:k) .* sin(lambda_n(1:k) * x) .* exp(-lambda_n(1:k).^2 * alpha * 1000));
        u_converge_1(k+1) = T0 + H*x + partial_sum_1;
        u_converge_1000(k+1) = T0 + H*x + partial_sum_1000;
    end
end

% Plot results
figure;
plot(n_series, u_converge_1, 'o-', 'LineWidth', 1.5);
hold on;
plot(n_series, u_converge_1000, 's-', 'LineWidth', 1.5);
xlabel('Number of terms n');
ylabel('Temperature u(x,t) [°C]');
title('Convergence of the Fourier Series Solution');
legend('t = 1 s', 't = 1000 s', 'Location', 'best');
grid on;