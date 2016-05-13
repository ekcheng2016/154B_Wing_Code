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
%   wx0 : force distribution in the x-direction from root to tip  (N/m)
%   wy0 : force distribution in the y-direction from root to tip  (N/m)
%   wx  : force distribution in the x-direction from tip to root  (N/m)
%   wy  : force distribution in the y-direction from tip to root  (N/m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [load] = calc_wxwy(n,rho,V,AoA,Cd,CM,nz)
load_aircraft_parameters;
load_conversions;

W = wgt_max;                            % N          weight

L = n*W;                                 % N          lift
CL = 2*L/rho/V^2/S;                      %            lift coefficient
D = 0.5*rho*V^2*S*(Cd + CL^2/pi/AR/e);   % N          drag 

M0 = 0.5*rho*V^2*S*CM*c;

z = 0:b/2/nz:b/2;                        % root to tip (half span)

% lift distribution
l_rect = L/b.*ones(1,nz+1);
l_ellip = (4*L/pi/b).*sqrt(1-(2*z./b).^2);
l = (l_rect + l_ellip)./2;               % N/m

% drag distribution
%  Assume drag is constant, when it is 20% from the tip, drag increases by
%  10%
tmp = abs(z - (0.8*0.5*b));
[val idx] = min(tmp);
d = [D*ones(1,idx) 1.1*D*ones(1,length(z)-idx)]/b;

% rotate into x-y coordinate
wy0 = cos(AoA).*l + sin(AoA).*d;
wx0 = -sin(AoA).*l + cos(AoA).*d;

% Note: wx0 and wy0 are defined from root to tip
% Note: wx and wy are defined from tip to root

wy = zeros(1,nz+1);
wx = zeros(1,nz+1);
for i = 1:nz+1
    wy(i) = wy0(nz+2-i);
    wx(i) = wx0(nz+2-i);
end

load.z = z;
load.wx0 = wx0;
load.wy0 = wy0;
load.wx = wx;
load.wy = wy;
load.l_ellip = l_ellip;
load.l_rect  = l_rect;
load.l       = l;
load.d       = d;
load.M0      = M0;