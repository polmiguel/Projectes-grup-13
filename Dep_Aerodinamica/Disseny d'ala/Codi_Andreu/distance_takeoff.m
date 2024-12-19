clc; clear all;

% Constants
m = 11654; 
rho = 1.225; 
Sdrag = pi * 1.28^2;
Slift = (3 + 1.95)/2 * 22;
C_L = 0.2;
mu_r = 0.03; 
C_D = 0.0622; 
V_despegue = 58.6466; 

% Functions
T = @(v) 3000000 * 0.7 ./ (1 + v); 
D = @(v) 0.5 * rho * C_D * Sdrag * v.^2;
L = @(v) 0.5 * rho * C_L * Slift * v.^2;
F_r = @(v) mu_r * (m * 9.81 - L(v));
F_net = @(v) T(v) - D(v) - F_r(v);
a = @(v) F_net(v) / m;

% Velocity and time calculation
v = linspace(0, V_despegue, 1000); 
t = cumtrapz(v, 1 ./ a(v)); % Time corresponding to each velocity
s = cumtrapz(v, v ./ a(v)); % Distance traveled

distancia_recorrida = s(end);

% Output
fprintf('Distance for V_despegue: %.2f m\n', distancia_recorrida);

% Graphs
figure;

% Subplot 1: Sum of forces vs velocity
subplot(3,1,1);
plot(v, F_net(v));
xlabel('Velocity (m/s)');
ylabel('Sum of forces (N)');
title('Sum of forces as a function of velocity');

% Subplot 2: Distance traveled vs velocity
subplot(3,1,2);
plot(v, s);
xlabel('Velocity (m/s)');
ylabel('Distance traveled (m)');
title('Distance traveled as a function of velocity');

% Subplot 3: Distance traveled vs time
subplot(3,1,3);
plot(t, s);
xlabel('Time (s)');
ylabel('Distance traveled (m)');
title('Distance traveled as a function of time');
