% calc_sigmazz.m
%
% Description:
%   This calculates the column buckling, skin buckling and yield and
%   fatigue of the wing structure.
%
% Inputs:
%       I : Moment Of Inertia
%       E : Young's Modulus
%       l_e: Effective Length
%       
%
% Outputs:
%       P_cr : Critical Load that causes buckling
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ buckling ] = calc_buckling( I,E, sigma_zz.max, A_stringer )

% Column Buckling
P_cr = 1.5*sigma_zz.max*A_stringer; 
l_e = sqrt(pi^2*E*I/P_cr); % represents the rib spacing

% Skin Buckling
k = (m*b/a + a/(mb))^2;
N_cr = k*pi^2*D/b^2;


