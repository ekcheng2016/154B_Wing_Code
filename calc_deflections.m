% calc_deflections.m
%
% Description:
%   calculate deflections along z-axis from root to tip
%
% Inputs:
%   alt  :          altitude
%   name :          name of loading case
%   b    :          full-span (m)
%   Ixx  :          area moment of inertia about x-axis (m^4)
%   Iyy  :          area moment of inertia about y-axis (m^4)
%   Ixy  :          area moment of inertia about xy-axis (m^4)
%   Mx0  :          moment about x-axis (1xM) from root to tip
%   My0  :          moment about y-axis (1xM) from root to tip
%   nz   :          number of values in the z-direction
%   wx0  :          load distribution (x) - validation purposes
%   wy0  :          load distribution (y) - validation purposes
%
% Outputs:
%   deflections :   structure containing deflection profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [deflections] = calc_deflections(alt,name,b,Ixx,Iyy,Ixy,Mx0,My0,nz,wx0,wy0)

    L = b/2;                % half span         m
    z = 0:L/nz:L;
    
    E = 70e9;   % modulus of material TODO: CHANGE LATER (don't hard code this here)
    K = 1/(E*(Ixx*Iyy - Ixy^2));
    ddu = -K.*(-Ixy.*Mx0 + Ixx.*My0);
    ddv = -K.*(Iyy.*Mx0 - Ixy.*My0);
    
    du = zeros(1,nz+1);
    du_sum = zeros(1,nz+1);
    dv = zeros(1,nz+1);
    dv_sum = zeros(1,nz+1);
    u = zeros(1,nz+1);
    u_sum = zeros(1,nz+1);
    v = zeros(1,nz+1);
    v_sum = zeros(1,nz+1);
    
    for i = 1:nz
        du(i+1) = (z(i+1)-z(i))*(ddu(i+1)+ddu(i))/2;
        dv(i+1) = (z(i+1)-z(i))*(ddv(i+1)+ddv(i))/2;
        du_sum(i+1) = du_sum(i) + du(i+1);
        dv_sum(i+1) = dv_sum(i) + dv(i+1);
        u(i+1) = (z(i+1)-z(i))*(du_sum(i+1)+du_sum(i))/2;
        v(i+1) = (z(i+1)-z(i))*(dv_sum(i+1)+dv_sum(i))/2;
        u_sum(i+1) = u_sum(i) + u(i+1);
        v_sum(i+1) = v_sum(i) + v(i+1);
    end
    
    fig = figure();
    hold on; grid on;
    subplot(1,2,1)
        plot(z,u_sum*1000,'LineWidth',2);
        xlim([0 L])
        xlabel('z (m)')
        ylabel('x-deflection (u) (mm)')
        title('Wing Deflections for : ','FontSize',14);
    subplot(1,2,2)
        plot(z,v_sum*1000,'LineWidth',2);
        xlim([0 L])
        xlabel('z (m)')
        ylabel('y-deflection (v) (mm)')
        title([alt ' ' name{1}],'FontSize',14);
    pos = get(fig,'position');
    set(fig,'position',[pos(1:2) pos(3)*1.5 pos(4)]);
    print(fig,[pwd '/Deflection_Figures/Deflection_' alt '_' name{1}],'-djpeg');
    
    deflections.z = z;
    deflections.u = u_sum;
    deflections.v = v_sum;
end