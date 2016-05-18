% load_wing_parameters.m
%
% Loads Cessna 177 Aircraft Parameters
%   - in MKS units
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load_conversions;

mu_sealvl  = 0.0000181206;   % Dynamic Viscosity at Sea Level   kg/ms
mu_altceil = 0.0000166411;   % Dynamic Viscosity at Ceiling     kg/ms

rho_sealvl  = 1.225000;       % Air Density at Sea Level     kg/m3
rho_altceil = 0.780933;       % Air Density at Ceiling       kg/m3

mass_max = 1100;                 % Max Gross Weight     kg
mass_emp = 680;                  % Empty Weight         kg
wgt_max  = mass_max*g;          % Max Gross Weight      N
wgt_emp  = mass_emp*g;          % Empty Weight          N
v_cruise = 230*km2m/hr2sec;     % Cruise Speed         m/s
v_maneuver = 250*km2m/hr2sec;   % Meaneuver Speed      m/s
alt_ceil = 4450;                % Service Ceiling      m
v_gust_cruise = 50;             % Gust Velocity Cruise ft/s
v_gust_dive   = 25;             % Gust Velocity Dive   ft/s

b = 10.82;          % Wing Span         m
c = 1.5;            % Chord Length      m
S = b*c;            % Wing Area         m^2
AR = b^2/S;          % Aspect Ratio      -
e = 0.79;           % Oswald Efficiency

n1      =  4.40;     % Positive Limit Maneuvering Load
n2      = -0.4*4.4;  % Negative Limit Maneuvering Load
n3      = -1;        % Negative Limit at Dive Speed
n_PHAA  = 4.40;      % Load at PHAA
n_NHAA  = -1.76;     % Load at NHAA

rho_material = 2780;          % Material Density kg/m^3 2024 T4 Aluminum
sigma_yield = 324;   % Yield Strength in MPa 2024 T4 Aluminum
