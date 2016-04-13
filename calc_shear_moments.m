% calc_shear_moments.m
%
% Description:
%   calculate shear force and moments along z-axis from root to tip
%
% Inputs:
%   alt:              altitude of the condition
%   name:             load condition
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
function [shear_moment] = calc_shear_moments(alt,name,b,nz,wx,wy,wx0,wy0)

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

fig = figure();
hold on; grid on;
subplot(2,3,1)
    plot(z,Mx0./1e3,'Linewidth',2)
    xlim([0 L])
    xlabel('z (m)')
    ylabel('Moments M_x (kN/m)')
subplot(2,3,2)
    plot(z,Sy0./1e3,'Linewidth',2)
    xlim([0 L])
    xlabel('z (m)')
    ylabel('Shear Force S_y (kN)')
    title([cellstr(alt) ' Shears, Moments and Load Distributions: ', cellstr(name)],'FontSize',14);
subplot(2,3,3)
    plot(z,wy0./1e3,'Linewidth',2)
    xlim([0 L])
    xlabel('z (m)')
    ylabel('Distributed Load w_y (kN/m)')
subplot(2,3,4)
    plot(z,My0./1e3,'Linewidth',2)
    xlim([0 L])
    xlabel('z (m)')
    ylabel('Moments M_y (kN/m)')
subplot(2,3,5)
    plot(z,Sx0./1e3,'Linewidth',2)
    xlim([0 L])
    xlabel('z (m)')
    ylabel('Shear Force S_x (kN)')
subplot(2,3,6)
    plot(z,wx0./1e3,'Linewidth',2)
    xlim([0 L])
    xlabel('z (m)')
    ylabel('Distributed Load w_x (kN/m)')
pos = get(fig, 'position');
set(fig,'position',[pos(1:2) pos(3:4)*2]);
print(fig,[pwd '/Load_Distribution_Figures/Shear_Moment_' alt '_' name{1}],'-djpeg');

shear_moment.z   = z;
shear_moment.Sx  = Sx_sum;
shear_moment.Sy  = Sy_sum;
shear_moment.Sx0 = Sx0;
shear_moment.Sy0 = Sy0;
shear_moment.Mx  = Mx_sum;
shear_moment.My  = My_sum;
shear_moment.Mx0 = Mx0;
shear_moment.My0 = My0;