% calc_shear_flow.m
%
% Description:
%   This calculates the shear flow of the wing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ output_args ] = calc_shear_flow(Ixx,Iyy,Ixy,airf_geo,Mx0,My0,...
                                Sx,Sy)

% dereference airfoil geometry structure
x       = airf_geo.x;
dx      = airf_geo.dx;
x_boomU = airf_geo.xU;
x_boomL = airf_geo.xL;
y_boomU = airf_geo.yU;
y_boomL = airf_geo.yL; 
L_boomU = sqrt((airf_geo.xU(1:end-1)-airf_geo.xU(2:end)).^2+...
                (airf_geo.yU(1:end-1)-airf_geo.yU(2:end)).^2);
L_boomL = sqrt((airf_geo.xL(1:end-1)-airf_geo.xL(2:end)).^2+...
                (airf_geo.yL(1:end-1)-airf_geo.yL(2:end)).^2);
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
BU(ind_capU(1:end-1)) = BU(ind_capU(1:end-1))+2*A_cap; % add 2 spar caps area
BU(ind_capU(end))     = BU(ind_capU(end))+A_cap;   % add 1 spar caps for last spar

BU(2:end-1) = BU(2:end-1) ...
            + t_skin*L_boomU(1:end-1)/6.*(2+sz_RBU(1:end-2)./sz_RBU(2:end-1)) ...
            + t_skin*L_boomU(2:end)/6.*(2+sz_RBU(3:end)./sz_RBU(2:end-1));
        
% ^above replace below
% for i = 2:nBU-1
%     BU(i) = BU(i) + t_skin*L_boomU(i-1)/6*(2+sz_RBU(i-1)/sz_RBU(i)) + t_skin*L_boomU(i)/6*(2+sz_RBU(i+1)/sz_RBU(i));
% end
BU(end) = BU(end) + t_skin*L_boomU(end)/6*(2+sz_RBU(end-1)/sz_RBU(end));


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
BL(ind_capL(1:end-1)) = BL(ind_capL(1:end-1))+2*A_cap;

BL(2:end-1) = BL(2:end-1)...
              + t_skin*L_boomL(1:end-1)/6.*(2+sz_RBL(1:end-2)./sz_RBL(2:end-1))...
              + t_skin*L_boomL(2:end)/6.*(2+sz_RBL(3:end)./sz_RBL(2:end-1));

% for i = 2:nBL-1
%     BL(i) = BL(i) + t_skin*L_boomL(i-1)/6*(2+sz_RBL(i-1)/sz_RBL(i)) + t_skin*L_boomL(i)/6*(2+sz_RBL(i+1)/sz_RBL(i));
% end
BL(end) = BL(end) + t_skin*L_boomL(end)/6*(2+sz_RBL(end-1)/sz_RBL(end));

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

% for i = 1:2
%     BU(i_BsparU(i)) = BU(i_BsparU(i)) + t_spar*h_spar(i)/6*(2+sz_RBL(i_BsparL(i))/sz_RBU(i_BsparU(i)));
%     BL(i_BsparL(i)) = BL(i_BsparL(i)) + t_spar*h_spar(i)/6*(2+sz_RBU(i_BsparU(i))/sz_RBL(i_BsparL(i)));
% end

% Calculate shear flow (LECTURE 4C pg 11)
K1 = (Sx(1)*Ixx-Sy(1)*Ixy)/(Ixx*Iyy-Ixy^2);
K2 = (Sy(1)*Iyy-Sx(1)*Ixy)/(Ixx*Iyy-Ixy^2);

qs_temp = 0;

for ii = 1:length(B_total)
    qs_temp = qs_temp - K1*B_total(ii)*xB_total(ii) - K2*B_total(ii)*yB_total(ii);
    qs(ii) = qs_temp;
end

figure()
plot(qs)
title('Shear Flow')
xlabel('boom index (counterclockwise'), ylabel('shear flow');

output_args = 1;

end

