% calc_random_wing.m
%
% Description:
%   This script calculates a random wing based off of the variation
%   parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load_variation_parameters;

airf_geo.A_cap = check_allowable(base_wing.A_cap,(rand-0.5)*variation.A_cap,0,2);
airf_geo.A_str = check_allowable(base_wing.A_str,(rand-0.5)*variation.A_str,0,2);
airf_geo.t_spar = check_allowable(base_wing.t_spar,(rand-0.5)*variation.t_spar,0,2);
airf_geo.t_skin = check_allowable(base_wing.t_skin,(rand-0.5)*variation.t_skin,0,2);
airf_geo.x_spar0 = check_allowable(base_wing.x_spar0,(rand-0.5)*variation.x_spar0,0.0000000001*c,0.75*c);

n_upperstr = floor(check_allowable(length(base_wing.x_strU0),(rand-0.5)*variation.n_upperstr,2,1000));
n_lowerstr = floor(check_allowable(length(base_wing.x_strL0),(rand-0.5)*variation.n_lowerstr,2,1000));

airf_geo.x_strU0 = [0 rand(1,n_upperstr-1)]*0.75*c; % randomly place upper stringers
airf_geo.x_strL0 = [0 rand(1,n_lowerstr-1)]*0.75*c; % randomly place lower stringers
