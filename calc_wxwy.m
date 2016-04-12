% calc_wxwy.m
%
% Description:
%   This function calculates the force distribution in the x and y
%   directions.  The notation is defined from ROOT to TIP
%
% Inputs:
%   n       : load factor
%   rho     : air density               (kg/m^3)
%   V       : velocity                  (m/s)
%   AoA     : angle of attack           (rad)
%   Cd      : airfoil drag coefficient  (-)
%   nz      : number of z increments    (-)
%
% Outputs:
%   wx : force distribution in the x-direction  (N/m)
%   wy : force distribution in the y-direction  (N/m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [wx,wy] = calc_wxwy(alt,name,n,rho,V,AoA,Cd,nz)
load_aircraft_parameters;
load_conversions;

W = wgt_max;                            % N          weight

L = n*W;                                 % N          lift
CL = 2*L/rho/V^2/S;                      %            lift coefficient
D = 0.5*rho*V^2*S*(Cd + CL^2/pi/AR/e);   % N          drag 

z = 0:b/2/nz:b/2;                        % root to tip (half span)

% lift distribution
l_rect = L/b.*ones(1,nz+1);
l_ellip = (4*L/pi/b).*sqrt(1-(2*z./b).^2);
l = (l_rect + l_ellip)./2;               % N/m

fig = figure();
hold on;
subplot(2,2,3)
    hold on;
    plot(z,l_ellip,'r',z,l_rect,'g','Linewidth',1)
    plot(z,l,'Linewidth',2)
    xlim([0 b/2]);
    xlabel('z (m)','FontSize',9)
    ylabel('Lift Distribution (N/m)','FontSize',9)
    title('Lift Distribution Half-Span Profile:','FontSize',9);
    if l(1) < 0
        legend({'lift elliptic distribution','lift rectangular distribution','combined lift distribution'},...
        'FontSize',8,'location','northwest')
    else
        legend({'lift elliptic distribution','lift rectangular distribution','combined lift distribution'},...
        'FontSize',8,'location','southeast')
    end
    
% drag distribution
%  Assume drag is constant, when it is 20% from the tip, drag increases by
%  10%
tmp = abs(z - (0.8*0.5*b));
[val idx] = min(tmp);
d = [D*ones(1,idx) 1.1*D*ones(1,length(z)-idx)]/b;

subplot(2,2,4)
    plot(z,d,'LineWidth',2)
    xlim([0 b/2]);
    xlabel('z (m)','FontSize',9)
    ylabel('Drag Distribution (N/m)','FontSize',9)
    title('Drag Distribution Half-Span Profile:','FontSize',9);
    legend({'drag distribution'}, 'FontSize',8,'location','west');

% rotate into x-y coordinate
wy = cos(AoA).*l + sin(AoA).*d;
wx = -sin(AoA).*l + cos(AoA).*d;

% Note: wx and wy are defined from root to tip

subplot(2,2,1:2)
    plot(z,wy,z,wx,'linewidth',2)
    xlim([0 b/2]);
    xlabel('z (m)','FontSize',9)
    ylabel('Distributed Load (N/m)','FontSize',9)
    title([' Load Distribution Half-Span Profile: ',cellstr([alt ' - ' name{1}])],'FontSize',9);
    legend({'w_y','w_x'},'FontSize',8)
pos = get(fig, 'position');
set(fig,'position',[pos(1:2) pos(3:4)*2]);
print(fig,[pwd '/Load_Distribution_Figures/' alt '_Load Distribution Half Span Profile_' name{1}],'-djpeg','-r300');
