clc;
close all;
clear;

Ixx = 0;
Iyy = 0;
Izz = 0;

omega_roll = 0;
omega_pitch = 0;
omega_yawn = 0;

distancia_alerons = 0;
d_deflector_hortizontal = 0;
d_deflector_vertical = 0;


% ROLL equation
F_alerons = (Ixx * omega_roll)/(2 * distancia_alerons);

% PITCH equation
F_deflector_horitzontal = (Iyy * omega_pitch)/d_deflector_horitzontal;

% YAWN equation
F_deflector_vertical = (Iyy * omega_yawn)/d_deflector_vertical;

% Force equation
Cl = 1.2; %Cl de placa plana a 10-12º 
F = Cl * 0.5 * rho * v^2 * S;
