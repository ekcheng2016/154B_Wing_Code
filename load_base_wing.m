% load_base_wing.m
%
% Description:
%   Loads base wing that will pass buckling and von mises
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_wing.A_cap = 7e-4;   % m^2
base_wing.A_str = 3.5e-5;   % m^2
base_wing.t_spar = 0.0025;   % m
base_wing.t_skin = 0.001016;  % m
% % locations of spars, spar caps and stringers (nose at the origin of the coordinate)
base_wing.x_spar0 = 0.25*c;                 % front spar (2 cell beam)
base_wing.x_strU0 = [0 0.05 0.15 0.35 0.55 0.65]*c; % upper surface
base_wing.x_strL0 = [0 0.05 0.15 0.35 0.55 0.65]*c; % lower surface