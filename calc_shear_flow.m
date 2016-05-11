% calc_shear_flow.m
%
% Description:
%   This calculates the shear flow of the wing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ output_args ] = calc_shear_flow(Ixx,Iyy,Ixy,airf_geo,Mx0,My0,...
                                Sx,Sy,M0,c,Cx,Cy)

% dereference airfoil geometry structure
x       = airf_geo.x;
dx      = airf_geo.dx;
yU      = airf_geo.yU;
yL      = airf_geo.yL;
x_boomU = airf_geo.xU;
x_boomL = airf_geo.xL;
y_boomU = airf_geo.y_skinU;
y_boomL = airf_geo.y_skinL;
%y_boomU = airf_geo.yU;
%y_boomL = airf_geo.yL; 
L_boomU = airf_geo.L_boomU;
L_boomL = airf_geo.L_boomL;
%L_boomU = sqrt((airf_geo.xU(1:end-1)-airf_geo.xU(2:end)).^2+...
%                (airf_geo.yU(1:end-1)-airf_geo.yU(2:end)).^2);
%L_boomL = sqrt((airf_geo.xL(1:end-1)-airf_geo.xL(2:end)).^2+...
%                (airf_geo.yL(1:end-1)-airf_geo.yL(2:end)).^2);
x_strU  = airf_geo.x_strU;
x_strL  = airf_geo.x_strL;
x_spar  = airf_geo.x_spar;
h_spar  = airf_geo.h_spar;
A_str   = airf_geo.A_str;
A_cap   = airf_geo.A_cap;
t_skin  = airf_geo.t_skin;
t_spar  = airf_geo.t_spar;

% boom area upper part
nBU = length(x_boomU);
i_BU = 1:nBU;

% stress at each boom on the top surface at the root of the wing
sz_RBU = zeros(1,nBU);
% sz_RBU(:) = Mx0(1)*(Iyy*y_boomU(i_BU(:))-Ixy*x_boomU(i_BU(:)))/(Ixx*Iyy-Ixy^2) +...
%             My0(1)*(Ixx*x_boomU(i_BU(:))-Ixy*y_boomU(i_BU(:)))/(Ixx*Iyy-Ixy^2);
% simplified version below
sz_RBU(:) = Mx0(1)*(Iyy*y_boomU-Ixy*x_boomU)/(Ixx*Iyy-Ixy^2) +...
            My0(1)*(Ixx*x_boomU-Ixy*y_boomU)/(Ixx*Iyy-Ixy^2);

% find boom node that is closest to stringers
for ii = 1:length(x_strU)
   [Y,ind] = min(abs(x_boomU - x_strU(ii)));
   x_strU_align(ii) = x_boomU(ind); % stringer x-position algined to closest boom
   ind_strU(ii) = ind;
end

% find boom node that is closest to spars
for ii = 1:length(x_spar)
   [Y,ind] = min(abs(x_boomU - x_spar(ii)));
   x_sparU_align(ii) = x_boomU(ind);
   ind_capU(ii) = ind; 
end

BU = zeros(1,nBU);

BU(ind_strU) = A_str; % add stringer area
% BU(ind_capU(1:end-1)) = BU(ind_capU(1:end-1))+2*A_cap; % add 2 spar caps area
BU(ind_capU(1:end-1)) = 2*A_cap;
% BU(ind_capU(end))     = BU(ind_capU(end))+A_cap;   % add 1 spar caps for last spar
BU(ind_capU(end))     = A_cap;

BU(2:nBU-1) = BU(2:nBU-1) ...
            + t_skin*L_boomU(1:nBU-2)/6.*(2+sz_RBU(1:nBU-2)./sz_RBU(2:nBU-1)) ...
            + t_skin*L_boomU(2:nBU-1)/6.*(2+sz_RBU(3:nBU)./sz_RBU(2:nBU-1));
        
% ^above replace below
% for i = 2:nBU-1
%     BU(i) = BU(i) + t_skin*L_boomU(i-1)/6*(2+sz_RBU(i-1)/sz_RBU(i)) + t_skin*L_boomU(i)/6*(2+sz_RBU(i+1)/sz_RBU(i));
% end
% BU(end) = BU(end) + t_skin*L_boomU(end)/6*(2+sz_RBU(end-1)/sz_RBU(end));
BU(end) = BU(end) + t_skin*L_boomU(end-1)/6*(2+sz_RBU(end-1)/sz_RBU(end));


% boom area lower part
nBL = length(x_boomL);
i_BL = 1:nBL;

% stress at each boom on the bottom surface at the root of the wing
sz_RBL = zeros(1,nBL);
% sz_RBL(:) = Mx0(1)*(Iyy*y_boomL(i_BL(:))-Ixy*x_boomL(i_BL(:)))/(Ixx*Iyy-Ixy^2)... 
%           + My0(1)*(Ixx*x_boomL(i_BL(:))-Ixy*y_boomL(i_BL(:)))/(Ixx*Iyy-Ixy^2);
sz_RBL(:) = Mx0(1)*(Iyy*y_boomL-Ixy*x_boomL)/(Ixx*Iyy-Ixy^2)... 
          + My0(1)*(Ixx*x_boomL-Ixy*y_boomL)/(Ixx*Iyy-Ixy^2);


% find boom node that is closest to stringers
for ii = 1:length(x_strL)
   [Y,ind] = min(abs(x_boomL - x_strL(ii)));
   x_strL_align(ii) = x_boomL(ind); % stringer x-position algined to closest boom
   ind_strL(ii) = ind;
end

for ii = 1:length(x_spar)
   [Y,ind] = min(abs(x_boomL - x_spar(ii)));
   x_sparL_align(ii) = x_boomL(ind);
   ind_capL(ii) = ind; 
end

% if_strL = ismember(x_boomL,x_strL);
% if_capL = ismember(x_boomL,x_spar);

BL = zeros(1,nBL);
BL(ind_strL) = A_str;
%BL(ind_capL(1:end-1)) = BL(ind_capL(1:end-1))+2*A_cap;
BL(ind_capL(1:end-1)) = 2*A_cap;
%BL(ind_capL(end))     = BL(ind_capL(end))+A_cap;   % add 1 spar caps for last spar
BL(ind_capL(end))     = A_cap;

BL(2:nBL-1) = BL(2:nBL-1)...
              + t_skin*L_boomL(1:nBL-2)/6.*(2+sz_RBL(1:nBL-2)./sz_RBL(2:nBL-1))...
              + t_skin*L_boomL(2:nBL-1)/6.*(2+sz_RBL(3:nBL)./sz_RBL(2:nBL-1));

% for i = 2:nBL-1
%     BL(i) = BL(i) + t_skin*L_boomL(i-1)/6*(2+sz_RBL(i-1)/sz_RBL(i)) + t_skin*L_boomL(i)/6*(2+sz_RBL(i+1)/sz_RBL(i));
% end
BL(end) = BL(end) + t_skin*L_boomL(end-1)/6*(2+sz_RBL(end-1)/sz_RBL(end));

% the first stringer shared by upper and lower surface
BU(1) = BU(1) ...
       + t_skin*L_boomU(1)/6*(2+sz_RBU(2)/sz_RBU(1)) ...
       + t_skin*L_boomL(1)/6*(2+sz_RBL(2)/sz_RBL(1));
BL(1) = BU(1);

% add the contribution of spars
% index for spar location
BU(ind_capU) = BU(ind_capU) + t_spar*h_spar/6.*(2+sz_RBL(ind_capL)./sz_RBU(ind_capU));
BL(ind_capU) = BL(ind_capU) + t_spar*h_spar/6.*(2+sz_RBU(ind_capU)./sz_RBL(ind_capL));

% combine top and bottom (start from top right, go counter clockwise)
B_total = [fliplr(BU), BL(2:end)];
xB_total = [fliplr(x_boomU), x_boomL(2:end)];
yB_total = [fliplr(y_boomU), y_boomL(2:end)];
LB_total = [fliplr(L_boomU), L_boomL, h_spar(2)];

figure()
plot(xB_total,yB_total,'ro','markersize',6);
hold on
plot(x,yU,'k',x,yL,'k','linewidth',1.5);
plot([x(end), x(end)],[y_boomU(end),y_boomL(end)],'b',[x_sparU_align(1),x_sparL_align(1)],[y_boomU(ind_capU(1)),y_boomL(ind_capL(1))],'b','linewidth',2)
ylim([-0.3 0.3]) 
xlabel('x (m)')
ylabel('y (m)')
title('Boom Distribution')
grid on

% for i = 1:2
%     BU(i_BsparU(i)) = BU(i_BsparU(i)) + t_spar*h_spar(i)/6*(2+sz_RBL(i_BsparL(i))/sz_RBU(i_BsparU(i)));
%     BL(i_BsparL(i)) = BL(i_BsparL(i)) + t_spar*h_spar(i)/6*(2+sz_RBU(i_BsparU(i))/sz_RBL(i_BsparL(i)));
% end


% calculate area of triangle formed by two nodes on the airfoil profile and point(x_spar(1),0)
nq = length(LB_total);
A = zeros(1,nq);
for i = 1:nq-1
    A(i) = abs(xB_total(i)*(yB_total(i+1)-0) + xB_total(i+1)*(0-yB_total(i)) + x_spar(1)*(yB_total(i)-yB_total(i+1)))/2;
end
A(end) = abs(xB_total(end)*(yB_total(1)-0) + xB_total(1)*(0-yB_total(end)) + x_spar(1)*(yB_total(end)-yB_total(1)))/2;
Asum = sum(A);

% calculate the area of each cell
[Y, i_A1(1)] = min(abs(xB_total(1:length(x_boomU))-x_spar(1)));
[Y, i_A1(2)] = min(abs(xB_total(length(x_boomU)+1:end)-x_spar(1)));
i_A1(2) = i_A1(2) + length(x_boomU);
A1 = A(i_A1(1):i_A1(2)-1);
A1sum = sum(A1);
A2sum = Asum - A1sum;

% calcualte qb at the root of the wing 
% Calculate shear flow (LECTURE 4C pg 11)
K1 = (Sx(1)*Ixx-Sy(1)*Ixy)/(Ixx*Iyy-Ixy^2);
K2 = (Sy(1)*Iyy-Sx(1)*Ixy)/(Ixx*Iyy-Ixy^2);

qb_temp = 0;

for ii = 1:length(B_total)
    qb_temp = qb_temp - K1*B_total(ii)*xB_total(ii) - K2*B_total(ii)*yB_total(ii);
    qb(ii) = qb_temp;
end

% separate cell 1 and cell 2 
qb1 = qb(i_A1(1):i_A1(2)-1);
L_boom1 = LB_total(i_A1(1):i_A1(2)-1);
L1sum = sum(L_boom1);

qb2 = [qb(1:i_A1(1)-1),qb(i_A1(2):end)];
L_boom2 = [LB_total(1:i_A1(1)-1),LB_total(i_A1(2):end)];
% modify distance between booms at the rear spar to compensate the change of the thickness
L_boom2(end) = L_boom2(end)*t_skin/t_spar; 
L2sum = sum(L_boom2);

% equations
syms q01 q02
% eq1: equating moments of applied shear and pitch moment to moments of internal shear flow
% find aerodynamic center
x_AC = (0.25*c)-Cx;
y_AC = 0-Cy;
[Y, ind_AC(1)] = min(abs(xB_total-x_AC));
[Y, ind_AC(2)] = min(abs(yB_total-y_AC));
x_s = xB_total(ind_AC(1));
y_s = yB_total(ind_AC(2));
eq1 = (2*A1sum*q01)+(2*A2sum*q02)+sum(2*qb(:).*A(:))-M0-(Sy(1)*x_s)-(Sx(1)*y_s);

% eq2: set the angle of twist of each cell to be the same
% Shear Flow Analysis PDF pg. 16
C1 = -K1*(sum(B_total(i_A1(1):i_A1(2)-1).*xB_total(i_A1(1):i_A1(2)-1))*L1sum/t_skin);
C2 = -K2*(sum(B_total(i_A1(1):i_A1(2)-1).*yB_total(i_A1(1):i_A1(2)-1))*L1sum/t_skin);
C3 = -K1*(sum(B_total(1:i_A1(1)-1).*xB_total(1:i_A1(1)-1).*LB_total(1:i_A1(1)-1))+...
          sum(B_total(i_A1(2):end).*xB_total(i_A1(2):end).*LB_total(i_A1(2):end)))/t_skin;
C4 = -K2*(sum(B_total(1:i_A1(1)-1).*yB_total(1:i_A1(1)-1).*LB_total(1:i_A1(1)-1))+...
          sum(B_total(i_A1(2):end).*yB_total(i_A1(2):end).*LB_total(i_A1(2):end)))/t_skin;

eq2 = (((q01*(L1sum/t_skin))-((q01-q02)*h_spar(1)/t_spar)+C1+C2)/A1sum)-...
      (((q02*(L2sum/t_skin))+((q02-q01)*h_spar(1)/t_spar)-(q02*h_spar(2)/t_spar)+C3+C4)/A2sum);

[q01,q02] = solve(eq1==0,eq2==0);

% shear flow along airfoil contour from top right corner to top right corner CCW 
q = zeros(1,nq);   
% shear flow from top right corner to boom index right before central spar
for i = 1:i_A1(1)-1
    q(i) = qb2(i) + q02;
end

% shear flow from top central spar around airfoil to bottom of central spar
for i = 1 : length(qb1)
    j = i + i_A1(1) - 1;
    q(j) = qb1(i) + q01;
end

% shear flow from bottom of central spar to bottom right corner
for i = 1 : nq - i_A1(2) + 1
    j = i + i_A1(2) - 1;
    k = i + i_A1(1) - 1;
    q(j) = qb2(k) + q02;
end

% TODO: calculate the shear flow in the central spar here

% verification:
% 1: Check whether last element of the array qb is zero or close to zero
if abs(qb(end)) <= 10
    output_args(1) = 1; % TRUE IF 1
else
    output_args(1) = 0;
end

% 2: Multiply stress at each boom by the corresponding boom areas and add
% them up, the sum should come to zero
chk = sum(sz_RBU(:).*BU(:))+sum(sz_RBL(2:end).*BL(2:end));
if chk == 0
   output_args(2) = 1;
else
   output_args(2) = 0;
end

% 3: Sx: please see pdf ?shear flow analysis? page 5 for details
% TODO:

% 4: Sy: please see pdf ?shear flow analysis? page 5 for details
% TODO:

% shear stress tau
%
% please calculate the shear stress from the shear flow
%


figure()
plot(qb)
title('Shear Flow')
xlabel('boom index (counterclockwise'), ylabel('shear flow');


end

