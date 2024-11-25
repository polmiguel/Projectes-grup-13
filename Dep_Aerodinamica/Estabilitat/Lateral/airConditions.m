function [T, P, rho] = airConditions(altitude)
    % Constants
    R = 287.05;  % Specific gas constant for dry air (J/kg·K)
    g = 9.80665; % Acceleration due to gravity (m/s²)
    T0 = 288.15; % Sea level standard temperature (K)
    P0 = 101325; % Sea level standard pressure (Pa)
    rho0 = 1.225; % Sea level standard air density (kg/m³)
    L = 0.0065;  % Temperature lapse rate (K/m)
    
    % Determine which layer of the atmosphere you're in (troposphere model)
    if altitude <= 11000 % Troposphere (up to 11 km)
        T = T0 - L * altitude; % Temperature at given altitude (K)
        P = P0 * (1 - L * altitude / T0)^(g / (R * L)); % Pressure at given altitude (Pa)
        rho = P / (R * T); % Density at given altitude (kg/m³)
    else
        error('Altitude exceeds 11 km; this code is for the troposphere only.');
    end
    
    % Display results
    fprintf('At an altitude of %.2f meters:\n', altitude);
    fprintf('Temperature = %.2f K\n', T);
    fprintf('Pressure = %.2f Pa\n', P);
    fprintf('Density = %.4f kg/m³\n', rho);
end
