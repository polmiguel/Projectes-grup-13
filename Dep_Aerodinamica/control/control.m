clc;
close all;
clear;

%Conditions data
rho = 0.9091; % 3000 m
v = 138.8; % m/s


%Structural data plane
xcg = 7.89;
ycg=0;
zcg = 0.3;
Ixx = 3.09 *10^5;
Iyy = 5.38 *10^4;
Izz = 8.51 *10^4;

%Regulations ratios
omega_roll = 15*2*pi/360; %15deg in 11 seconds
omega_pitch = 20*2*pi/360; %10 deg in 1 second
omega_yawn = 15*2*pi/360; %10 deg in 1 second

%Plane geometry
pos_ales = 11;
distancia_alerons = 15;
d_deflector_horitzontal = 20;
d_deflector_vertical = 20;


% ROLL equation
F_alerons = (Ixx * omega_roll)/(2 * (distancia_alerons));

% PITCH equation
F_deflector_horitzontal = (Iyy * omega_pitch)/(d_deflector_horitzontal-xcg);

% YAWN equation
F_deflector_vertical = (Iyy * omega_yawn)/(d_deflector_vertical-xcg);

% Force equation
Cl_aleiron = 0.2;  
Cl_elevator = 0.2;
Cl_rudder = 0.2;

S_alerons = F_alerons/(Cl_aleiron * 0.5 * rho * v^2)
S_deflector_horitzontal = 1.5*F_deflector_horitzontal/(Cl_elevator * 0.5 * rho * (v*0.8)^2)
S_deflector_vertical = 1.5*F_deflector_vertical/(Cl_rudder * 0.5 * rho * (v*0.8)^2)
