% calc_stress.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Description:
%   This calculates the stress of the airfoil and produces a contour plot of
%   the airfoil cross section.
%
% Inputs:
%        Mx
%        My
%        Ixx
%        Iyy
%        Ixy
%
% Outputs:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[stress] = calc_stress(Mx,My,Ixx,Ixx,Ixy);

sigma_z = Mx*(Iyy*y-Ixy*x)/(Ixx*Iyy-Ixy^2)+...
          My*(Ixx*x-Ixy*y)/(Ixx*Iyy-Ixy^2);

