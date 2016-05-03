% calc_shear_flow.m
%
% Description:
%   This calculates the shear flow of the wing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ output_args ] = calc_shear_flow(Ixx,Iyy,Ixy,airf_geo,Mx0,My0)

% dereference airfoil geometry structure
x       = airf_geo.x;
dx      = airf_geo.dx;
x_boomU = airf_geo.xU;      % this is ok, double check this
x_boomL = airf_geo.xL;      % this is ok, double check this
yU      = airf_geo.yU;
yL      = airf_geo.yL;
x_strU  = airf_geo.x_strU;
x_strL  = airf_geo.x_strL;
x_spar  = airf_geo.x_spar;
L_boomU = airf_geo.L_boomU; % this is wrong
L_boomL = airf_geo.L_boomL; % this is wrong
h_spar  = airf_geo.h_spar;
A_str   = airf_geo.A_str;
A_cap   = airf_geo.A_cap;
t_skin  = airf_geo.t_skin;
t_spar  = airf_geo.t_spar;
                            
% boom area upper part
nBU = length(x_boomU);
i_BU = zeros(1,nBU);
i_BU(:) = floor((x_boomU(:)-x_boomU(1))/dx)+1

% stress at each boom on the top surface at the root of the wing
sz_RBU = zeros(1,nBU);
sz_RBU(:) = Mx0(1)*(Iyy*yU(i_BU(:))-Ixy*x(i_BU(:)))/(Ixx*Iyy-Ixy^2) +...
            My0(1)*(Ixx*x(i_BU(:))-Ixy*yU(i_BU(:)))/(Ixx*Iyy-Ixy^2);

if_strU = ismember(x_boomU,x_strU);
if_capU = ismember(x_boomU,x_spar);

BU = zeros(1,nBU);
for i = 1:nBU
    if if_strU(i) == 1
        BU(i) = A_str;
    elseif if_capU(i) == 1
        BU(i) = A_cap;
    end
end

for i = 2:nBU-1
    BU(i) = BU(i) + t_skin*L_boomU(i-1)/6*(2+sz_RBU(i-1)/sz_RBU(i)) + t_skin*L_boomU(i)/6*(2+sz_RBU(i+1)/sz_RBU(i));
end
BU(end) = BU(end) + t_skin*L_boomU(end-1)/6*(2+sz_RBU(end-1)/sz_RBU(end));


% boom area lower part
nBL = length(x_boomL);
i_BL = zeros(1,nBL);
i_BL(:) = floor((x_boomL(:)-x_boomL(1))/dx)+1;

% stress at each boom on the bottom surface at the root of the wing
sz_RBL = zeros(1,nBL);
sz_RBL(:) = Mx0(1)*(Iyy*yL(i_BL(:))-Ixy*x(i_BL(:)))/(Ixx*Iyy-Ixy^2) + My0(1)*(Ixx*x(i_BL(:))-Ixy*yL(i_BL(:)))/(Ixx*Iyy-Ixy^2);

if_strL = ismember(x_boomL,x_strL);
if_capL = ismember(x_boomL,x_spar);

BL = zeros(1,nBL);
for i = 1:nBL
    if if_strL(i) == 1
        BL(i) = A_str;
    elseif if_capL(i) == 1
        BL(i) = A_cap;
    end
end

for i = 2:nBL-1
    BL(i) = BL(i) + t_skin*L_boomL(i-1)/6*(2+sz_RBL(i-1)/sz_RBL(i)) + t_skin*L_boomL(i)/6*(2+sz_RBL(i+1)/sz_RBL(i));
end
BL(end) = BL(end) + t_skin*L_boomL(end-1)/6*(2+sz_RBL(end-1)/sz_RBL(end));

% the first stringer shared by upper and lower surface
BU(1) = BU(1) + t_skin*L_boomU(1)/6*(2+sz_RBU(2)/sz_RBU(1)) + t_skin*L_boomL(1)/6*(2+sz_RBL(2)/sz_RBL(1));
BL(1) = BU(1);

% add the contribution of spars
% index for spar location
i_BsparU = find(if_capU);
i_BsparL = find(if_capL);
for i = 1:2
    BU(i_BsparU(i)) = BU(i_BsparU(i)) + t_spar*h_spar(i)/6*(2+sz_RBL(i_BsparL(i))/sz_RBU(i_BsparU(i)));
    BL(i_BsparL(i)) = BL(i_BsparL(i)) + t_spar*h_spar(i)/6*(2+sz_RBU(i_BsparU(i))/sz_RBL(i_BsparL(i)));
end

end

