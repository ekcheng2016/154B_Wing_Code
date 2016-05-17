% calc_sigmazz.m
%
% Description:
%   This calculates the stress in the z-direction at the root of the wing
%   z=0 in units of Pa (N/m2)
%
% Inputs:
%   Ixx : area moment of inertia in the x-direction (m^4)
%   Iyy : area moment of inertia in the y-direction (m^4)
%   Ixy : area momoent of inertia in the xy-direction (m^4)
%   Mx  : moment in the x-direction at root (N/m)
%   My  : moment in the y-direction at root (N/m)
%   xU  : x-coordinates of upper surface (m)
%   yU  : y-coordinates of upper surface (m)
%   xL  : x-coordinates of lower surface (m)
%   yL  : y-coordinates of lower surface (m)
%
% Outputs:
%   sigma_zz_U: sigma_zz at the root on the upper surface of wing (Pa)
%   sigma_zz_L: sigma_zz at the root on the lower surface of wing (Pa)
%   sigma_zz_max : max absolute value direct stress (Pa)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ sigma_zz ] = calc_sigmazz( Ixx,Iyy,Ixy,Mx,My,xU,yU,xL,yL )

for i = 1:length(xU)
    sigma_zz_U(i) = (Mx*(Iyy*yU(i) - Ixy*xU(i)))/(Ixx*Iyy - Ixy^2) + ...
                    (My*(Ixx*xU(i) - Ixy*yU(i)))/(Ixx*Iyy - Ixy^2);
end

for i = 1:length(xL)
    sigma_zz_L(i) = (Mx*(Iyy*yL(i) - Ixy*xL(i)))/(Ixx*Iyy - Ixy^2) + ...
                    (My*(Ixx*xL(i) - Ixy*yL(i)))/(Ixx*Iyy - Ixy^2);
end

sigma_zz.upper = sigma_zz_U;
sigma_zz.lower = sigma_zz_L;

sigma_zz.max   = max([sigma_zz_U sigma_zz_L]);
sigma_zz.min   = min([sigma_zz_U sigma_zz_L]);
end

