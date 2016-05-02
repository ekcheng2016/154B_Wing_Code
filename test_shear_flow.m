% test_shear_flow.m 
%   This file is currently under development.
% 

clear all
close all
clc

load_aircraft_parameters
load('Airfoil_Data/NACA2415_coordinates');
xy_scaled = c*xy_NACA_2415; %[m]
idx_zero = find(xy_scaled(:,1) == 0);

% design parameters
t = 0.0015; %[m] airfoil thickness (assumed the same number as EXCEL)
A_cap = 0.0001; %[m^2] spar cap area
N_str = 10; % number of stringers
x_spar = [0.3 0.75]*c; %[m] spar locations. Note: last element must be the
                       % spar next to control surface

twice_deltaA = xy_scaled(1:end-1,1).*xy_scaled(2:end,2)-...
               xy_scaled(2:end,1).*xy_scaled(1:end-1,2); %cross product of two consecutive points
twice_A = sum(twice_deltaA); %[m^2] twice the cross sectional area
d_panel = sqrt((xy_scaled(1:end-1,1)-xy_scaled(2:end,1)).^2+...
            (xy_scaled(1:end-1,2)-xy_scaled(2:end,2)).^2); %[m] distances between each nodes on panel
delta_A_skin = d_panel*t;

% for each spar, find corresponding node
for ii = 1:length(x_spar)
    
    % upper end of spars
    [val, idx] = min(abs(xy_scaled(1:idx_zero,1)-x_spar(ii)));
    x_spar_upper(ii) = xy_scaled(idx,1); %[m] x-location of upper panel node closest to spar1 
    y_spar_upper(ii) = xy_scaled(idx,2);
    idx_spar_upper(ii) = idx;
    
    %   find index immediately RIGHT of where spar would go on UPPER surface
    [idx] = find((xy_scaled(1:idx_zero,1)-x_spar(ii))>=0, 1,'last');
    idx_spar_Panel_upper(ii) = idx;
    
    % lower end of spars
    [val, idx] = min(abs(xy_scaled(idx_zero:end,1)-x_spar(ii)));
    x_spar_lower(ii) = xy_scaled(idx+idx_zero-1,1); %[m]xy_ x-location of lower panel node closest to spar1 
    y_spar_lower(ii) = xy_scaled(idx+idx_zero-1,2);
    idx_spar_lower(ii) = idx+idx_zero-1;
    
    % find index immediately LEFT of where spar would go on LOWER surface
    % (compared to RIGHT for UPPER surface)
    [idx] = find((xy_scaled(idx_zero:end,1)-x_spar(ii))<=0, 1,'last');
    idx_spar_Panel_lower(ii) = idx+idx_zero-1;
    
end

% Spar cap area vector
delta_A_cap = zeros(size(delta_A_skin));
delta_A_cap(idx_spar_Panel_upper(1:end-1)) = 2*A_cap; %[m^2] 
delta_A_cap(idx_spar_Panel_upper(end)) = A_cap; %[m^2] 
delta_A_cap(idx_spar_Panel_lower(1:end-1)) = 2*A_cap; %[m^2] 
delta_A_cap(idx_spar_Panel_lower(end)) = A_cap; %[m^2] 

% Stringer area vector
delta_A_str = zeros(size(delta_A_skin));
%TODO: implement stringer area


%


figure() 
plot(xy_scaled(:,1), xy_scaled(:,2)), hold on, grid on
for ii = 1:length(x_spar)
    plot([x_spar_upper(ii) x_spar_lower(ii)],[y_spar_upper(ii) y_spar_lower(ii)],'r')
end
title('NACA2415 Airfoil Contour');
xlabel('x(m)'), ylabel('y(m)'),axis('equal')

