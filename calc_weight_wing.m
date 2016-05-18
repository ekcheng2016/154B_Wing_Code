%calc_weight_wing.m
%   
% Description:
%   Calculate the wing weight based on material of the wing and geometric
%   properties of the wing.
% 
% Inputs:
%   airf_geo : geometric qualities of the wing
%   b        : span (full) m
%   rho      : density of the material kg/m3
%
% Outputs:
%   weight_wing : wing weight (in kg)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ weight_wing ] = calc_weight_wing( airf_geo, b, rho )

A_cap = airf_geo.A_cap;     % spar cap areas
num_str = length([airf_geo.x_strU airf_geo.x_strL]-1);  % -1 because nose stringer is shared by both surfaces
A_str = airf_geo.A_str;     % stringer areas
t_spar = airf_geo.t_spar;   % spar thickness
h_spar = airf_geo.h_spar;   % height of the spar
t_skin = airf_geo.t_skin;   % skin thickness
L_skin = sum(sqrt((airf_geo.x(1:end-1)-airf_geo.x(2:end)).^2+...
                (airf_geo.yU(1:end-1)-airf_geo.yU(2:end)).^2))+...
         sum(sqrt((airf_geo.x(1:end-1)-airf_geo.x(2:end)).^2+...
                (airf_geo.yL(1:end-1)-airf_geo.yL(2:end)).^2));  % MORE BOOMS

% number of caps = 6 (4 in the central spar, 2 in the rear spar)
weight_cap = 6*rho*b*A_cap;
weight_str = num_str*rho*b*A_str;
weight_spar = rho*b*t_spar*sum(h_spar);
weight_skin = rho*b*t_skin*L_skin;

weight_total = weight_cap+weight_str+weight_spar+weight_skin;

weight_wing.cap   = weight_cap;
weight_wing.str   = weight_str;
weight_wing.spar  = weight_spar;
weight_wing.skin  = weight_skin;
weight_wing.total = weight_total;


end

