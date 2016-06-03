% calc_sigmazz.m
%
% Description:
%   This calculates the column buckling, skin buckling and yield and
%   fatigue of the wing structure.
%
% Inputs:
%       I_str    : Area moment of inertia of the stringers
%       E        : Young's Modulus
%       sigma_zz_max : Max tensile stress of the wing section
%       sigma_zz_min : Max compressive stress of the wing section
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
function [ buckling ] = calc_buckling( I_str, sigma_zz_max,sigma_zz_min,A_str,t , airf_geo)

% dereference
i_strU = airf_geo.i_strU;
i_strL = airf_geo.i_strL;
i_spar = airf_geo.i_spar;

i_U = sort([i_strU i_spar]);
i_L = sort([i_strL i_spar]);

x = airf_geo.x;
y_U = airf_geo.yU;
y_L = airf_geo.yL;

b_upper = sqrt((y_U(i_U(2:end))-y_U(i_U(1:end-1))).^2 + (x(i_U(2:end))-x(i_U(1:end-1))).^2);
b_lower = sqrt((y_L(i_L(2:end))-y_L(i_L(1:end-1))).^2 + (x(i_L(2:end))-x(i_L(1:end-1))).^2);

b = sort([b_upper b_lower]);

E = 73.1*1000; % MPa

%% Column Buckling
P_max = 1.5*abs(sigma_zz_min)*A_str; % max loading of the wing section
l_e = sqrt(pi^2*E*I_str/P_max); % effective length of the spacing 
rib_l = 2*l_e; % rib spacing 

%% Skin Buckling
mu = 0.33; % possoin's ratio
k = 8.5; % from the chart in class 8-9
% b = 0.20; % TODO: ARBITRARY PICK. VALUE BETWEEN STRINGER
sigma_crit = (k*pi^2*E/(12*(1-mu^2))).*(t./b).^2; % critical stress for skin buckling

%% Yield & Fatigue
%
K_IC = 26; % material property of 2024 T4 Al
c = 1.6E-11; % material property of 2024 T4 Al
m = 3.59; % material property of 2024 T4 Al

% Wing Section
f = 1; % stress intensity factor
a_crit = (K_IC/(f*sigma_zz_max))^2/pi; % in m
a_i = 0.1*a_crit; % assumes to be 10% of a_crit in m
n = 2/(c*(m-2)*((sigma_zz_max-sigma_zz_min)*pi^0.5)^m)*...
    (1/(a_i^((m-2)/2))-1/(a_crit^((m-2)/2))); % number of cycles until Failure

%% Important Fail or Pass Values
buckling.rib_l = rib_l; % rib spacing in m
buckling.n = n; % cycles to fatigue in # of cycles
buckling.a_crit = a_crit; % critical crack length in m
buckling.sigma_crit = min(sigma_crit); % critical stress in MPa
