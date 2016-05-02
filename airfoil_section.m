% airfoil_section.m
%
% Description:
%   This calculates the airfoil profile and the thin-walled moment of
%   inertia approximations
%
% Inputs:
%   c:          chord length                      (m)
%   A_cap:      cap area                          (m^2)
%   A_str:      stringer area                     (m^2)
%   t_spar:     spar thickness                    (m)
%   t_skin:     skin thickness                    (m)
%   x_spar:     spar locations(1xM)               (m)
%   x_strU:     upper stringer locations(1xN)     (m)
%   x_strL:     lower stringer locations(1xP)     (m)
%
% Outputs:
%   Cx:         x-centroid of wing section        (m)
%   Cy:         y-centroid of wing section        (m)
%   Ixx:        xx-area moment of inertia         (m^4)
%   Iyy:        yy-area moment of inertia         (m^4)
%   Ixy:        xy-area moment of inertia         (m^4)
%   xU:         x-coordinate of upper surface from centroid    (m)
%   yU:         y-coordinate of upper surface from centroid    (m)
%   xL:         x-coordinate of lower surface from centroid    (m)
%   yL:         y-coordinate of lower surface from centroid    (m)
%   x_strU:     x-coordinate of upper stringers from centroid  (m)
%   x_strL:     x-coordinate of lower stringers from centroid  (m)
%   x_spar:     x-coordinate of spars from centroid            (m)
%   h_spar:     spar heights                                   (m)
%   i_spar:     spar indices                                   (-)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Cx,Cy,Ixx,Iyy,Ixy,xU,yU,xL,yL,x_strU,x_strL,x_spar,h_spar,i_spar] = ...
    airfoil_section(c, A_cap,A_str,t_spar,t_skin,x_spar,x_strU,x_strL,PLOT_AIRFOIL)
%% airfoil section profile
% NACA 2415
m = 0.02;
p = 0.4;
t = 0.15;

n = 500;            % number of increments
dx = c/n;
x = 0:dx:c;         % even spacing

yc = zeros(1,n+1);
yt = zeros(1,n+1);
xU0 = zeros(1,n+1);
yU0 = zeros(1,n+1);
xL0 = zeros(1,n+1);
yL0 = zeros(1,n+1);
theta = zeros(1,n+1);
xb = p*c;
i_xb = xb/dx + 1;  

for i = 1:n+1
    % MEAN THICKNESS LINE yt
    yt(i) = 5*t*c*(0.2969*sqrt(x(i)/c) - 0.1260*(x(i)/c) - 0.3516*(x(i)/c)^2 + 0.2843*(x(i)/c)^3 - 0.1015*(x(i)/c)^4);
    
    % MEAN CAMBER LINE yc
    if i <= i_xb
        yc(i) = m*x(i)/p^2*(2*p - x(i)/c);
        theta(i) = atan(2*m/p^2*(p - x(i)/c));
    else
        yc(i) = m*(c - x(i))/(1-p)^2*(1 + x(i)/c - 2*p);
        theta(i) = atan(2*m/(1-p)^2*(p-x(i)/c));
    end
    
    % U == upper surface x,y coordinates
    % L == lower surface x,y, coordinates
    xU0(i) = x(i) - yt(i)*sin(theta(i));
    yU0(i) = yc(i) + yt(i)*cos(theta(i));
    xL0(i) = x(i) + yt(i)*sin(theta(i));
    yL0(i) = yc(i) - yt(i)*cos(theta(i));
end

% The last 25% of the chord length of the airfoil was neglected under the assumption that
% this section contained flaps and ailerons, and would therefore not support aerodynamic loads.

i_xend = round(0.75*c/dx)+1;
x = x(1:i_xend);
xU = xU0(1:i_xend);
yU = yU0(1:i_xend);
xL = xL0(1:i_xend);
yL = yL0(1:i_xend);

% Here the airfoil profile is approximated by assuming xU0=x & xL0=x.

%% adding stringers, spar caps and spars
% spars & spar caps
x_spar = [x_spar,x(end)];         % add the 3rd spar
i_spar = round(x_spar./dx)+1;        
h_spar = yU(i_spar) - yL(i_spar);
Cy_spar = (yU(i_spar) + yL(i_spar))/2;
A_spar = t_spar.*h_spar;

% stringers
n_strU = length(x_strU);
i_strU = round(x_strU./dx)+1;
n_strL = length(x_strL);
i_strL = round(x_strL./dx)+1;

% skins
x_skinU = [x_spar,x_strU];
x_skinU = sort(x_skinU);
i_skinU = round(x_skinU/dx)+1;
n_skinU = length(x_skinU)-1;
L_skinU = zeros(1,n_skinU);
A_skinU = zeros(1,n_skinU);
Cx_skinU = zeros(1,n_skinU);
Cy_skinU = zeros(1,n_skinU);
for i = 1:n_skinU
    L_skinU(i) = sqrt((yU(i_skinU(i+1)) - yU(i_skinU(i)))^2 + (x_skinU(i+1) - x_skinU(i))^2);
    A_skinU(i) = t_skin*L_skinU(i);
    Cx_skinU(i) = (x_skinU(i+1) + x_skinU(i))/2;
    Cy_skinU(i) = (yU(i_skinU(i+1)) + yU(i_skinU(i)))/2;
end

x_skinL = [x_spar,x_strL];
x_skinL = sort(x_skinL);
i_skinL = round(x_skinL/dx)+1;
n_skinL = length(x_skinL)-1;
L_skinL = zeros(1,n_skinL);
A_skinL = zeros(1,n_skinL);
Cx_skinL = zeros(1,n_skinL);
Cy_skinL = zeros(1,n_skinL);
for i = 1:n_skinL
    L_skinL(i) = sqrt((yL(i_skinL(i+1)) - yL(i_skinL(i)))^2 + (x_skinL(i+1) - x_skinL(i))^2);
    A_skinL(i) = t_skin*L_skinL(i);
    Cx_skinL(i) = (x_skinL(i+1) + x_skinL(i))/2;
    Cy_skinL(i) = (yL(i_skinL(i+1)) + yL(i_skinL(i)))/2;
end

%% centroid of the wing section
Cx_sum = 0;
Cy_sum = 0;
A_sum = 0;

% 2 spars
for i = 1:2
    Cx_sum = Cx_sum + x_spar(i)*A_spar(i);
    Cy_sum = Cy_sum + Cy_spar(i)*A_spar(i);
    A_sum = A_sum + A_spar(i);
end

% spar caps
for i = 1:2
    Cx_sum = Cx_sum + x_spar(i)*A_cap;
    Cy_sum = Cy_sum + yU(i_spar(i))*A_cap;
    Cy_sum = Cy_sum + yL(i_spar(i))*A_cap;
    A_sum  = A_sum  + 2*A_cap;
end

% upper stringers
for i = 1:length(n_strU)
    Cx_sum = Cx_sum + x_strU(i)*A_str(i);
    Cy_sum = Cy_sum + Cy_skinU(i)*A_str(i);
    A_sum = A_sum + A_str(i);
end

% lower stringers
for i = 1:length(n_strL)
    Cx_sum = Cx_sum + x_strL(i)*A_str(i);
    Cy_sum = Cy_sum + Cy_skinL(i)*A_str(i);
    A_sum = A_sum + A_str(i);
end

% upper skin
for i = 1:n_skinU
   Cx_sum = Cx_sum + x_skinU(i)*A_skinU(i);
   Cy_sum = Cy_sum + Cy_skinU(i)*A_skinU(i);
   A_sum = A_sum + A_skinU(i);
end

% lower skin
for i = 1:n_skinL
   Cx_sum = Cx_sum + x_skinL(i)*A_skinL(i);
   Cy_sum = Cy_sum + Cy_skinL(i)*A_skinL(i);
   A_sum = A_sum + A_skinL(i);
end

Cx = Cx_sum/A_sum;
Cy = Cy_sum/A_sum;

if PLOT_AIRFOIL
    afc = figure();
    plot(x,yU,'k',x,yL,'k','Linewidth',2);
    ylim([-0.3 0.3])
    hold on; grid on;
    plot(x_strU,yU(i_strU),'or',x_strL,yL(i_strL),'or','markersize',5);
    plot([x_spar(1),x_spar(1)],[yU(i_spar(1)),yL(i_spar(1))],'b',[x_spar(2),x_spar(2)],[yU(i_spar(2)),yL(i_spar(2))],'b',[x(end),x(end)],[yU(end),yL(end)],'b','Linewidth',3);
    plot(x_spar,yU(i_spar),'sg',x_spar,yL(i_spar),'sg','markersize',6);
    scatter(Cx,Cy,'m*')
    title('NACA 2415 with Origin at Nose','fontsize',14);
    ylabel('y (m)')
    xlabel('x (m)')
    print(afc,[pwd '/Airfoil_Section/REAL_Airfoil_Origin'],'-djpeg','-r300');
end

%% Area moments of inertia
Ixx = 0;
Iyy = 0;
Ixy = 0;

% 3 spars
for i = 1:2
    Ixx = Ixx + t_spar*h_spar(i)^3/12 + A_spar(i)*(Cy_spar(i)-Cy)^2;
    Iyy = Iyy + t_spar^3*h_spar(i)/12 + A_spar(i)*(x_spar(i)-Cx)^2;
    Ixy = Ixy + A_spar(i)*(Cy_spar(i)-Cy)*(x_spar(i)-Cx);
end

% spar caps (2 in the rear, 4 in the front)
Ixx = Ixx + 2*((A_cap*(yU(i_spar(1))-Cy)^2)+(A_cap*(yL(i_spar(1))-Cy)^2));
Ixx = Ixx + ((A_cap*(yU(i_spar(2))-Cy)^2)+(A_cap*(yL(i_spar(2))-Cy)^2));
Iyy = Iyy + 4*A_cap*(x_spar(1)-Cx)^2;    % 4 spar caps - front location
Iyy = Iyy + 2*A_cap*(x_spar(2)-Cx)^2;    % 2 spar caps - rear location
Ixy = Ixy + 2*A_cap*(yU(i_spar(1))-Cy)*(x_spar(1)-Cx); % upper spar cap - front
Ixy = Ixy + 2*A_cap*(yL(i_spar(1))-Cy)*(x_spar(1)-Cx); % lower spar cap - front
Ixy = Ixy + A_cap*(yU(i_spar(2))-Cy)*(x_spar(2)-Cx); % upper spar cap - rear
Ixy = Ixy + A_cap*(yL(i_spar(2))-Cy)*(x_spar(2)-Cx); % lower spar cap - rear
   
% upper stringers
for i = 1:length(n_strU)
    Ixx = Ixx + A_str*(Cy_skinU(i)-Cy)^2;
    Iyy = Iyy + A_str*(x_strU(i)-Cx)^2;
    Ixy = Ixy + A_str*(Cy_skinU(i)-Cy)*(x_strU(i)-Cx);
end

% lower stringers
for i = 1:length(n_strU)
    Ixx = Ixx + A_str*(Cy_skinL(i)-Cy)^2;
    Iyy = Iyy + A_str*(x_strL(i)-Cx)^2;
    Ixy = Ixy + A_str*(Cy_skinL(i)-Cy)*(x_strL(i)-Cx);
end

% upper skin
for i = 1:n_skinU
   Ixx = Ixx + A_skinU(i)*(Cy_skinU(i)-Cy)^2;
   Iyy = Iyy + A_skinU(i)*(x_skinU(i)-Cx)^2;
   Ixy = Ixy + A_skinU(i)*(Cy_skinU(i)-Cy)*(x_skinU(i)-Cx);
end

% lower skin
for i = 1:n_skinL
   Ixx = Ixx + A_skinL(i)*(Cy_skinL(i)-Cy)^2;
   Iyy = Iyy + A_skinL(i)*(x_skinL(i)-Cx)^2;
   Ixy = Ixy + A_skinL(i)*(Cy_skinL(i)-Cy)*(x_skinL(i)-Cx);
end

%% convert coordinates to centroid as origin
xU = xU - Cx;
yU = yU - Cy;
xL = xL - Cx;
yL = yL - Cy;
x_strU = x_strU - Cx;
x_strL = x_strL - Cx;
x_spar = x_spar - Cx;

if PLOT_AIRFOIL
    afc = figure();
    plot(xU,yU,'k',xL,yL,'k','Linewidth',2);
    ylim([-0.3 0.3])
    hold on; grid on;
    plot(x_strU,yU(i_strU),'or',x_strL,yL(i_strL),'or','markersize',5);
    plot([x_spar(1),x_spar(1)],[yU(i_spar(1)),yL(i_spar(1))],'b',[x_spar(2),x_spar(2)],[yU(i_spar(2)),yL(i_spar(2))],'b',[x(end),x(end)]-Cx,[yU(end),yL(end)],'b','Linewidth',3);
    plot(x_spar,yU(i_spar),'sg',x_spar,yL(i_spar),'sg','markersize',6);
    scatter(Cx-Cx,Cy-Cy,'m*')
    title('NACA 2415 Normalized to Centroid','fontsize',14);
    ylabel('y (m)')
    xlabel('x (m)')
    print(afc,[pwd '/Airfoil_Section/REAL_Airfoil_Centroid'],'-djpeg','-r300');
end