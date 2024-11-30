%% Prop
C_D= 0.006; % En creuer
C_L = 0.2;% En creuer
S = 54.45;% Superficie alar (m^2)

T_cruise % Thrust en cruise
D_cruise % Drag en cruise

m_0 % Massa sense fuel ni bateríes

h_cruise = 5000; % m
v_cruise = 500; % km/h

%% Aero

CDG %Centre de gravetat
I %Matriu d'inercia
L_fus %Dimensions fuselatge

%% Estrus
Chord_root = 3; %Corda del perfil a la root
Chord_tip = 1.95; %Corda del perfil a la tip
b = 22; %Envergadura de l'ala
twist = -5; % (deg) rotació lineal de la root fins la tip de l'ala

vel_creuer = 138.8; %(m/s) velocitat de creuer
vel_takeoff = 55; %(m/s) velocitat de takeoff

