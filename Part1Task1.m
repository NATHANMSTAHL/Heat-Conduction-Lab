clear;
clc;
close all;


function[T_0, Hexp, SS_Temps] = SS_slope_T_0 (filename)

data = readmatrix(filename);

time = data(:,1);
CH1 = data(:,2);
CH2 = data(:,3);
CH3 = data(:,4);
CH4 = data(:,5);
CH5 = data(:,6);
CH6 = data(:,7);
CH7 = data(:,8);
CH8 = data(:,9);

ss_time_idx = length(time) - 2;

SS1 = CH1(ss_time_idx);
SS2 = CH2(ss_time_idx);
SS3 = CH3(ss_time_idx);
SS4 = CH4(ss_time_idx);
SS5 = CH5(ss_time_idx);
SS6 = CH6(ss_time_idx);
SS7 = CH7(ss_time_idx);
SS8 = CH8(ss_time_idx);

% CH_Temps = [CH1 CH2 CH3 CH4 CH5 CH6 CH7 CH8];

% x0 = 0;
x1 = 0.034925;
delta_x = 0.0127;
x_values = [x1 x1+delta_x x1+2*delta_x x1+3*delta_x x1+4*delta_x x1+5*delta_x x1+6*delta_x x1+7*delta_x];


SS_Temps = [SS1 SS2 SS3 SS4 SS5 SS6 SS7 SS8];

[polycoeff] = polyfit(x_values, SS_Temps, 1);

Hexp = polycoeff(1);
T_0 = polycoeff(2);

end


x1 = 0.034925;
delta_x = 0.0127;
x_values = [x1 x1+delta_x x1+2*delta_x x1+3*delta_x x1+4*delta_x x1+5*delta_x x1+6*delta_x x1+7*delta_x];

x = linspace(0, max(x_values), 100);

AL_25V_data = 'Aluminum_25V_240mA';  % A
AL_30V_data = 'Aluminum_30V_290mA';  % B 
Brass_25V_data = 'Brass_25V_237mA';  % C
Brass_30V_data = 'Brass_30V_285mA';  % D
Steel_22V_data = 'Steel_22V_203mA';  % E

[T_0_A, Hexp_A, SS_Temps_A] = SS_slope_T_0(AL_25V_data);
f_x_A = Hexp_A*x + T_0_A;

[T_0_B, Hexp_B, SS_Temps_B] = SS_slope_T_0(AL_30V_data);
f_x_B = Hexp_B*x + T_0_B;

[T_0_C, Hexp_C, SS_Temps_C] = SS_slope_T_0(Brass_25V_data);
f_x_C = Hexp_C*x + T_0_C;

[T_0_D, Hexp_D, SS_Temps_D] = SS_slope_T_0(Brass_30V_data);
f_x_D = Hexp_D*x + T_0_D;

[T_0_E, Hexp_E, SS_Temps_E] = SS_slope_T_0(Steel_22V_data);
f_x_E = Hexp_E*x + T_0_E;

Hexp_values = [Hexp_A Hexp_B Hexp_C Hexp_D Hexp_E];
T_0_values = [T_0_A T_0_B T_0_C T_0_D T_0_E];


V = [25 30 25 30 22];
I = [.24 .29 .237 .285 .203];
Qdot = V.*I;
r = 0.0127;
A = pi*r^2;
k = [130 130 115 115 16.2];
H_an = Qdot./(k*A);



figure();
hold on;
scatter(x_values, SS_Temps_A, 40, 'filled', 'ro');
plot(x, f_x_A, 'r', 'LineWidth', 1.5);
plot(x, H_an(1)*x + T_0_A, 'b', 'LineWidth', 1.5)
title('Steady State Temperatures vs. Rod Location (Aluminum, 25V, 240 mA)');
xlabel('x (m)');
ylabel('Temperature (^\circC)');
legend('Experimental Data', 'Experimental Data Extrapolation', 'Analytical Data Extrapolation', 'Location','best');
hold off
exportgraphics(gcf, 'Aluminum_25V_240mA.png', 'Resolution', 300);

figure();
hold on;
scatter(x_values, SS_Temps_B, 40, 'filled', 'ro');
plot(x, f_x_B, 'r', 'LineWidth', 1.5);
plot(x, H_an(2)*x + T_0_B, 'b', 'LineWidth', 1.5)
title('Steady State Temperatures vs. Rod Location (Aluminum, 30V, 290 mA)');
xlabel('x (m)');
ylabel('Temperature (^\circC)');
legend('Experimental Data', 'Experimental Data Extrapolation', 'Analytical Data Extrapolation', 'Location','best');
hold off
exportgraphics(gcf, 'Aluminum_30V_290mA.png', 'Resolution', 300);

figure();
hold on;
scatter(x_values, SS_Temps_C, 40, 'filled', 'ro');
plot(x, f_x_C, 'r', 'LineWidth', 1.5);
plot(x, H_an(3)*x + T_0_C, 'b', 'LineWidth', 1.5)
title('Steady State Temperatures vs. Rod Location (Brass, 25V, 237 mA)');
xlabel('x (m)');
ylabel('Temperature (^\circC)');
legend('Experimental Data', 'Experimental Data Extrapolation', 'Analytical Data Extrapolation', 'Location','best');
hold off
exportgraphics(gcf, 'Brass_25V_237mA.png', 'Resolution', 300);

figure();
hold on;
scatter(x_values, SS_Temps_D, 40, 'filled', 'ro');
plot(x, f_x_D, 'r', 'LineWidth', 1.5);
plot(x, H_an(4)*x + T_0_D, 'b', 'LineWidth', 1.5)
title('Steady State Temperatures vs. Rod Location (Brass, 30V, 285 mA)');
xlabel('x (m)');
ylabel('Temperature (^\circC)');
legend('Experimental Data', 'Experimental Data Extrapolation', 'Analytical Data Extrapolation', 'Location','best');
hold off
exportgraphics(gcf, 'Brass_30V_285mA.png', 'Resolution', 300);

figure();
hold on;
scatter(x_values, SS_Temps_E, 40, 'filled', 'ro');
plot(x, f_x_E, 'r', 'LineWidth', 1.5);
plot(x, H_an(5)*x + T_0_E, 'b', 'LineWidth', 1.5)
title('Steady State Temperatures vs. Rod Location (Steel, 22V, 203 mA)');
xlabel('x (m)');
ylabel('Temperature (^\circC)');
legend('Experimental Data', 'Experimental Data Extrapolation', 'Analytical Data Extrapolation', 'Location','best');
hold off
exportgraphics(gcf, 'Steel_22V_203mA.png', 'Resolution', 300);