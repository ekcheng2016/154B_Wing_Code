% calc_shear_moments.m
%
% Description:
%   calculate shear force and moments along z-axis from root to tip
%
% Inputs:
%   b  :              span (full) (m)
%   nz :              number of z-increments
%   wx :              x-load distribution tip to root (N/m)
%   wy :              y-load distribution tip to root (N/m)
%   wx0:              x-load distribution root to tip (N/m)
%   wy0:              y-load distribution root to tip (N/m)
%
% Outputs:
%   shear_moment:     structure containing the shear and moment profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [shear,moment] = calc_shear_moments(b,nz,wx,wy,wx0,wy0)

L = b/2;                % half span         m
z = 0:L/nz:L;

Sx = zeros(1,nz+1);
Sx_sum = zeros(1,nz+1);
Sy = zeros(1,nz+1);
Sy_sum = zeros(1,nz+1);

Mx = zeros(1,nz+1);
Mx_sum = zeros(1,nz+1);
My = zeros(1,nz+1);
My_sum = zeros(1,nz+1);

for i = 1:nz;
    Sx(i+1) = -(z(i+1)-z(i))*(wx(i+1) + wx(i))/2;
    Sy(i+1) = -(z(i+1)-z(i))*(wy(i+1) + wy(i))/2;
    Sx_sum(i+1) = Sx_sum(i) + Sx(i+1);
    Sy_sum(i+1) = Sy_sum(i) + Sy(i+1);
    Mx(i+1) = (z(i+1)-z(i))*(Sy_sum(i+1) + Sy_sum(i))/2;
    My(i+1) = (z(i+1)-z(i))*(Sx_sum(i+1) + Sx_sum(i))/2;
    Mx_sum(i+1) = Mx_sum(i) + Mx(i+1);
    My_sum(i+1) = My_sum(i) + My(i+1);
end

Sx_sum = -Sx_sum;
Sy_sum = -Sy_sum;

% z-->L-z: from root to tip
Mx0 = zeros(1,nz+1);
My0 = zeros(1,nz+1);
Sx0 = zeros(1,nz+1);
Sy0 = zeros(1,nz+1);
for i = 1:nz+1
    Mx0(i) = Mx_sum(nz+2-i);
    My0(i) = My_sum(nz+2-i);
    Sx0(i) = Sx_sum(nz+2-i);
    Sy0(i) = Sy_sum(nz+2-i);
end

shear.z   = z;
moment.z  = z;
shear.Sx  = Sx_sum;
shear.Sy  = Sy_sum;
shear.Sx0 = Sx0;
shear.Sy0 = Sy0;
moment.Mx  = Mx_sum;
moment.My  = My_sum;
moment.Mx0 = Mx0;
moment.My0 = My0;