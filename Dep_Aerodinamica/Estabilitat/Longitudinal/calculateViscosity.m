function mu = calculateViscosity(T)
    % Constants
    mu_0 = 1.827e-5; % Reference dynamic viscosity in kg/(m·s)
    T_0 = 291.15;    % Reference temperature in K
    C = 120;         % Sutherland's constant for air in K

    % Sutherland's formula
    mu = mu_0 * (T / T_0)^(3/2) * (T_0 + C) / (T + C);
end