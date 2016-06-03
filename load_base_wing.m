% load_base_wing.m
%
% Description:
%   Loads base wing that will pass buckling and von mises
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

airf_geo.A_cap = 7e-4;   % m^2
airf_geo.A_str = 5e-6;   % m^2
airf_geo.t_spar = 0.0025;   % m
airf_geo.t_skin = 0.0038;  % m
% % locations of spars, spar caps and stringers (nose at the origin of the coordinate)
airf_geo.x_spar0 = 0.25*c;                 % front spar (2 cell beam)
airf_geo.x_strU0 = [0 0.05 0.15 0.35 0.55 0.65]*c; % upper surface
airf_geo.x_strL0 = [0 0.05 0.15 0.35 0.55 0.65]*c; % lower surface
