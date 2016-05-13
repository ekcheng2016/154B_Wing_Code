% calc_sigmazz.m
%
% Description:
%   This calculates the column buckling, skin buckling and yield and
%   fatigue of the wing structure.
%
% Inputs:
%       I        : Area moment of inertia of the stringers
%       E        : Young's Modulus
%       sigma_zz_max : Max direct stress of the wing section
%       A_str    : Stringer Area
%
% Outputs:
%       P_cr : Critical Load that causes column buckling
%       sig_cr : Critical Stress that causes skin buckling
%       l_e  : Effective length 
%       l    : Rib spacing
%       
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ buckling ] = calc_buckling( I_str, sigma_zz_max,A_str )

E = 1.5; % material property

%% Column Buckling
P_max = 1.5*sigma_zz_max*A_str; % max loading of the wing section
l_e = sqrt(pi^2*E*I_str/P_max); % effective length of the spacing 
l = 2*l_e; % rib spacing 

%% Skin Buckling
mu = 0.33; % possoin's ratio
k = 8.5; % from the chart in class 8-9
sigma_crit = k*pi^2*E/(12*(1-mu^2))*(t/b)^2; % critical stress for skin buckling

%% Yield & Fatigue

% Material Property
K_IC = 26; % material property of 2024 T3 Al
c = 1.6E-11; % material property of 2024 T3 Al
m = 3.59; % material property of 2024 T3 Al

% Wing Section
f = 1; % stress intensity factor
a_crit = (K_IC/(f*sigma_zz_max))^2/pi; % in m
a_i = 0; % in m
n = (2/m)*(a_crit^(m/2)-a_i^(m/2))/(c*(sigma_zz_max*sqrt(pi))^m); % number of cycles to fatigue




