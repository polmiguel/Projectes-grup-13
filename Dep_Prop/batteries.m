%% Data
% Clear variables and screen
clear; clc;

% General parameters
C_D = 0.008; % Drag coefficient
C_L = 0.32;  % Lift coefficient
S = 36.3;    % Wing area in m^2

% Initial mass without batteries or fuel
m0 = 8000; % kg

% Combustion engines
eff_comb = 0.2684;

% Batteries
E_bat = 11.1; % kWh per battery
m_bat = 48;   % kg per battery
V_bat = 0.063; % m^3 per battery

% Fuel
P_c = 44.65e6; % J/kg (44.65 MJ/kg)
d_comb = 0.8;  % kg/L
SFC = 0.408; % kg/kWh

% Coolant
refrigerant_ratio = 0.10; % Coolant volume/battery volume ratio
densidad_refrigerante = 1070; % kg/m^3

% Generators
P_generator = 400e3; % Power of Safran generator in W (400 kW)
eff_generator = 0.85; % Efficiency of the generator (85%)
m_generator = 34; % Mass of each generator in kg
N_generators = 2; % Number of generators (one per turboprop)

% Function for air density as a function of altitude
rho = @(h) 1.225 * (1 - (0.0065 * h) / 288.15).^4.2561;

% Flight parameters
v_taxi = 7.5;           % Taxi speed in m/s
t_taxi_out = 12.3 * 60; % Taxi-out time in seconds
mu = 0.02;              % Rolling friction coefficient
v_cruise = 500 / 3.6;   % Cruise speed in m/s
h_cruise = 7600;        % Cruise altitude in m

% Takeoff and landing
theta_to = deg2rad(8);  % Climb angle in radians
t_to = 600;             % Takeoff time in seconds
theta_land = deg2rad(4); % Descent angle in radians
t_land = 20 * 60;       % Landing time in seconds

% Distances for cruise
delta_x_to = h_cruise / tan(theta_to);      % Horizontal distance during climb
delta_x_land = h_cruise / tan(theta_land);  % Horizontal distance during descent
x_elec = 2e5;   % Total distance for electric mode in m
x_hybrid = 8e5; % Total distance for hybrid mode in m

t_cruise_e = (x_elec - delta_x_to - delta_x_land) / v_cruise; % Electric cruise time
t_cruise_h = (x_hybrid - delta_x_to - delta_x_land) / v_cruise; % Hybrid cruise time

% Convergence
tol = 1; % kg
max_iter = 1000;
m = m0; % Initial mass
m_prev = m0 - 2 * tol; % Ensure at least one iteration
iter = 0;

%% Iteration for final mass
while abs(m - m_prev) > tol && iter < max_iter
    iter = iter + 1;
    m_prev = m;

    % Total weight
    W = m * 9.81; % Weight in N

    %% Energy and power
    % Taxi-out
    rho_taxi = rho(0); % Air density at sea level
    D_taxi = 0.5 * rho_taxi * v_taxi^2 * S * C_D; % Drag
    L_taxi = 0.5 * rho_taxi * v_taxi^2 * S * C_L; % Lift
    T_taxi_out = D_taxi + mu * (W - L_taxi);      % Thrust
    P_taxi = T_taxi_out * v_taxi;                % Power
    E_taxi_out = P_taxi * t_taxi_out;            % Taxi-out energy

    % Takeoff
    v_to_vals = linspace(v_taxi, v_cruise, 1000); % Speed during climb
    D_to_vals = 0.5 * rho(h_cruise) .* v_to_vals.^2 * S * C_D;
    T_to_vals = D_to_vals + W * sin(theta_to);
    P_to_vals = T_to_vals .* v_to_vals;
    E_to = trapz(linspace(0, t_to, 1000), P_to_vals); % Takeoff energy
    P_to_avg = mean(P_to_vals); % Average power during takeoff

    % Cruise
    D_cruise = 0.5 * rho(h_cruise) * v_cruise^2 * S * C_D;
    T_cruise = D_cruise; % Level flight
    P_cruise = T_cruise * v_cruise;

    % Hybrid cruise power (including generator power)
    P_hybrid = P_cruise/eff_comb; % Total power required during hybrid cruise
    E_cruise_h = P_hybrid * (t_cruise_h - t_cruise_e); % Hybrid cruise energy

    % Electric cruise energy
    E_cruise_e = P_cruise * t_cruise_e; % Electric cruise energy

    % Descent
    v_land_vals = linspace(v_cruise, 62, 1000); % Speed during descent
    D_land_vals = 0.5 * rho(h_cruise) .* v_land_vals.^2 * S * C_D;
    T_land_vals = D_land_vals;
    P_land_vals = T_land_vals .* v_land_vals;
    E_land = trapz(linspace(0, t_land, 1000), P_land_vals); % Descent energy
    P_land_avg = mean(P_land_vals); % Average power during landing

    % Taxi-in
    E_taxi_in = P_taxi * t_taxi_out;

    %% Total energies
    E_elec = E_taxi_out + E_to + E_cruise_e + E_land + E_taxi_in; % Total electric energy (J)
    E_total = E_taxi_out + E_to + E_cruise_h + E_land + E_taxi_in; % Total energy required (J)

    % Energy from generators during hybrid cruise
    E_gen = P_generator * N_generators * t_cruise_h * eff_generator; % Total energy from generators (J)

    %% Fuel consumption
    % Total power required from engines during hybrid cruise
    P_engines_hybrid = P_hybrid; % Already includes generator power input
    m_comb_cruise = SFC * (P_engines_hybrid / 1000) * (t_cruise_h / 3600); % Fuel mass during hybrid cruise (kg)

    % Fuel mass is calculated based on the power required and SFC
    m_comb = m_comb_cruise; % Total fuel mass (kg)
    V_comb = m_comb / d_comb; % Fuel volume (L)

    %% Batteries
    E_elec_kWh = E_elec / (1000 * 3600); % Convert energy to kWh
    N_bat = ceil(E_elec_kWh / E_bat);    % Number of batteries
    m_bat_tot = N_bat * m_bat;
    V_bat_tot = N_bat * V_bat;

    %% Mass and volume
    % Coolant
    V_refrigerante = V_bat_tot * refrigerant_ratio;
    m_refrigerante = V_refrigerante * densidad_refrigerante;

    % Total mass
    m = m0 + m_bat_tot + m_refrigerante + m_comb + (N_generators * m_generator);

    fprintf('Iteration %d: m = %.2f kg, m_bat_tot = %.2f kg, m_comb = %.2f kg\n', iter, m, m_bat_tot, m_comb);
end

% Final results
if iter == max_iter
    disp('No convergence achieved.');
else
    disp('Convergence achieved.');
end

fprintf('Final mass: %.2f kg\n', m);
fprintf('Number of batteries: %d\n', N_bat);
fprintf('Fuel mass: %.2f kg\n', m_comb);
fprintf('Fuel volume: %.2f L\n', V_comb);
fprintf('Energy from generators: %.2f kWh\n', E_gen / (1000 * 3600));
fprintf('Average power during takeoff: %.2f kW\n', P_to_avg / 1000); % Convert to kW
fprintf('Average power during landing: %.2f kW\n', P_land_avg / 1000); % Convert to kW
