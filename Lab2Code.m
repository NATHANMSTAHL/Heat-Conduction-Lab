clear; clc; close all;

V = [25 30 25 30 22];
I = [.24 .29 .237 .285 .203];
Qdot = V.*I;
r = 0.0127;
A = pi*r^2;
k = [130 130 115 115 16.2];
H_an = Qdot./(k*A);