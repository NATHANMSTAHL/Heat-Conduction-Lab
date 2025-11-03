clc;
clear;
close all;

function [u_model, time, T_exp, x_values] = u_x_t (filename, H, T0, k, rho, c_p)

data = readmatrix(filename);

data(end-1:end, :) = [];

time = data(:,1);
CH1 = data(:,2);
CH2 = data(:,3);
CH3 = data(:,4);
CH4 = data(:,5);
CH5 = data(:,6);
CH6 = data(:,7);
CH7 = data(:,8);
CH8 = data(:,9);

T_exp = [CH1 CH2 CH3 CH4 CH5 CH6 CH7 CH8];

x1 = 0.034925;
delta_x = 0.0127;
x_values = [x1 x1+delta_x x1+2*delta_x x1+3*delta_x x1+4*delta_x x1+5*delta_x x1+6*delta_x x1+7*delta_x];

L = 0.149225;   

alpha = k / (rho * c_p);


% Define lambda_n and b_n
% n_max = 10;
% n = 1:n_max;
n = 1:10;
lambda_n = (2*n - 1) * pi / (2*L); 


% for n = 1:n_max
b_n(n) = (8*H*L)/(pi^2) .* ((-1.^n)./(2*n-1).^2);
% end


Nt = length(time);
u_model = zeros(8, Nt);

for i = 1:8                      % loop TCs
    x = x_values(i);
    for j = 1:Nt                 % loop times
        t = time(j);
        sum_term = sum(b_n .* sin(lambda_n * x) .* exp( - (lambda_n.^2) * alpha * t ));
        u_model(i,j) = T0 + H*x + sum_term;
    end
end



end

% Al 25V
H_an_al_25 = 91.09;
H_ex_al_25 = 53.4;
T0_al_25 = 17.04;
k_al_adj_25 = 105;
rho_al = 2810;
c_p_al_adj_25 = 1300; 
alpha_adj_25 = k_al_adj_25 / (rho_al * c_p_al_adj_25);
% k_al = 130;
% rho_al = 2810;
% c_p_al = 960; 

%Al 30V
H_an_al_30 = 132.1;
H_ex_al_30 = 79.21;
T0_al_30 = 17.18;
k_al_adj_30 = 115;
% rho_al = 2810;
c_p_al_adj_30 = 1300; 
alpha_adj_30 = k_al_adj_30 / (rho_al * c_p_al_adj_30);

%Brass 25V
H_an_br_25 = 101.7;
H_ex_br_25 = 105.6;
T0_br_25 = 16.5;
k_br = 50;
rho_br = 8500;
c_p_br = 380;
% k_br = 115;
% rho_br = 8500;
% c_p_br = 380;

%Brass 30V
H_an_br_30 = 146.7;
H_ex_br_30 = 150.4;
T0_br_30 = 16.75;

%Steel 22V
H_an_st_22 = 544.1;
H_ex_st_22 = 287.3;
T0_st_22 = 15.11;
k_st_adj_22 = 20.25;
rho_st = 8000;
c_p_st_adj_22 = 500;
% k_st = 16.2;
% rho_st = 8000;
% c_p_st = 500;
alpha_adj_22 = k_st_adj_22 / (rho_st * c_p_st_adj_22);

AL_25V_data = 'Aluminum_25V_240mA';  % A
AL_30V_data = 'Aluminum_30V_290mA';  % B 
Brass_25V_data = 'Brass_25V_237mA';  % C
Brass_30V_data = 'Brass_30V_285mA';  % D
Steel_22V_data = 'Steel_22V_203mA';  % E


% experimental
[u_model_ex_A, time_ex_A, T_exp_ex_A, x_values_ex_A] = u_x_t (AL_25V_data, H_ex_al_25, T0_al_25, k_al_adj_25, rho_al, c_p_al_adj_25);
[u_model_ex_B, time_ex_B, T_exp_ex_B, x_values_ex_B] = u_x_t (AL_30V_data, H_ex_al_30, T0_al_30, k_al_adj_30, rho_al, c_p_al_adj_30);
[u_model_ex_C, time_ex_C, T_exp_ex_C, x_values_ex_C] = u_x_t (Brass_25V_data, H_ex_br_25, T0_br_25, k_br, rho_br, c_p_br);
[u_model_ex_D, time_ex_D, T_exp_ex_D, x_values_ex_D] = u_x_t (Brass_30V_data, H_ex_br_30, T0_br_30, k_br, rho_br, c_p_br);
[u_model_ex_E, time_ex_E, T_exp_ex_E, x_values_ex_E] = u_x_t (Steel_22V_data, H_ex_st_22, T0_st_22, k_st_adj_22, rho_st, c_p_st_adj_22);







figure;
t = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
sgtitle('Transient Temperature Profiles - Model III (Exp. Steady-State Slope vs Exp. Data with \alpha_a_d_j)')
% Keep axes handles to manage legend reliably
ax = gobjects(5,1);

% 1) Aluminum 25V
ax(1) = nexttile(1);
hold on;
plot(time_ex_A, u_model_ex_A(1,:), 'b', 'LineWidth', 1.5);
plot(time_ex_A, u_model_ex_A(2,:), 'b', 'LineWidth', 1.5);
plot(time_ex_A, u_model_ex_A(3,:), 'b', 'LineWidth', 1.5);
plot(time_ex_A, u_model_ex_A(4,:), 'b', 'LineWidth', 1.5);
plot(time_ex_A, u_model_ex_A(5,:), 'b', 'LineWidth', 1.5);
plot(time_ex_A, u_model_ex_A(6,:), 'b', 'LineWidth', 1.5);
plot(time_ex_A, u_model_ex_A(7,:), 'b', 'LineWidth', 1.5);
plot(time_ex_A, u_model_ex_A(8,:), 'b', 'LineWidth', 1.5);
plot(time_ex_A, T_exp_ex_A, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Temperature u(x,t) [°C]');
title('Aluminum 25V'); 
grid on;

% 2) Aluminum 30V
ax(2) = nexttile(2);
hold on;
plot(time_ex_B, u_model_ex_B(1,:), 'b', 'LineWidth', 1.5);
plot(time_ex_B, u_model_ex_B(2,:), 'b', 'LineWidth', 1.5);
plot(time_ex_B, u_model_ex_B(3,:), 'b', 'LineWidth', 1.5);
plot(time_ex_B, u_model_ex_B(4,:), 'b', 'LineWidth', 1.5);
plot(time_ex_B, u_model_ex_B(5,:), 'b', 'LineWidth', 1.5);
plot(time_ex_B, u_model_ex_B(6,:), 'b', 'LineWidth', 1.5);
plot(time_ex_B, u_model_ex_B(7,:), 'b', 'LineWidth', 1.5);
plot(time_ex_B, u_model_ex_B(8,:), 'b', 'LineWidth', 1.5);
plot(time_ex_B, T_exp_ex_B, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Temperature u(x,t) [°C]');
title('Aluminum 30V'); 
grid on;

% 3) Brass 25V
ax(3) = nexttile(3);
hold on;
plot(time_ex_C, u_model_ex_C(1,:), 'b', 'LineWidth', 1.5);
plot(time_ex_C, u_model_ex_C(2,:), 'b', 'LineWidth', 1.5);
plot(time_ex_C, u_model_ex_C(3,:), 'b', 'LineWidth', 1.5);
plot(time_ex_C, u_model_ex_C(4,:), 'b', 'LineWidth', 1.5);
plot(time_ex_C, u_model_ex_C(5,:), 'b', 'LineWidth', 1.5);
plot(time_ex_C, u_model_ex_C(6,:), 'b', 'LineWidth', 1.5);
plot(time_ex_C, u_model_ex_C(7,:), 'b', 'LineWidth', 1.5);
plot(time_ex_C, u_model_ex_C(8,:), 'b', 'LineWidth', 1.5);
plot(time_ex_C, T_exp_ex_C, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Temperature u(x,t) [°C]');
title('Brass 25V'); 
grid on;

% 4) Brass 30V
ax(4) = nexttile(4);
hold on;
plot(time_ex_D, u_model_ex_D(1,:), 'b', 'LineWidth', 1.5);
plot(time_ex_D, u_model_ex_D(2,:), 'b', 'LineWidth', 1.5);
plot(time_ex_D, u_model_ex_D(3,:), 'b', 'LineWidth', 1.5);
plot(time_ex_D, u_model_ex_D(4,:), 'b', 'LineWidth', 1.5);
plot(time_ex_D, u_model_ex_D(5,:), 'b', 'LineWidth', 1.5);
plot(time_ex_D, u_model_ex_D(6,:), 'b', 'LineWidth', 1.5);
plot(time_ex_D, u_model_ex_D(7,:), 'b', 'LineWidth', 1.5);
plot(time_ex_D, u_model_ex_D(8,:), 'b', 'LineWidth', 1.5);
plot(time_ex_D, T_exp_ex_D, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Temperature u(x,t) [°C]');
title('Brass 30V'); 
grid on;

% 5) Steel 22V (span last row)
ax(5) = nexttile([1 2]);
hold on;
plot(time_ex_E, u_model_ex_E(1,:), 'b', 'LineWidth', 1.5);
plot(time_ex_E, u_model_ex_E(2,:), 'b', 'LineWidth', 1.5);
plot(time_ex_E, u_model_ex_E(3,:), 'b', 'LineWidth', 1.5);
plot(time_ex_E, u_model_ex_E(4,:), 'b', 'LineWidth', 1.5);
plot(time_ex_E, u_model_ex_E(5,:), 'b', 'LineWidth', 1.5);
plot(time_ex_E, u_model_ex_E(6,:), 'b', 'LineWidth', 1.5);
plot(time_ex_E, u_model_ex_E(7,:), 'b', 'LineWidth', 1.5);
plot(time_ex_E, u_model_ex_E(8,:), 'b', 'LineWidth', 1.5);
plot(time_ex_E, T_exp_ex_E, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Temperature u(x,t) [°C]');
title('Steel 22V'); 
grid on;

% Create one legend; attach to first axes (doesn't matter which) then move it
lgd = legend(ax(1), {'Experimental Model', '', '', '', '', '', '', '','Experimental Data'}, 'Orientation','horizontal');

try
    lgd.NumColumns = 2;
    lgd.Layout.Tile = 'south';          % Works in newer MATLAB with tiledlayout
catch
    % Fallback: manual placement under tiles
    set(lgd,'Units','normalized','Position',[0.3 0.01 0.4 0.05]);
end

