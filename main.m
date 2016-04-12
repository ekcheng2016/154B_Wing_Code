clear all; 
close all; 
clc;

load_aircraft_parameters;
load_conversions;

Re_sealvl = calc_Re(rho_sealvl,c,v_maneuver,mu_sealvl);
Re_alceil = calc_Re(rho_altceil,c,v_maneuver,mu_altceil);

airfoil_to_wing;
n_allow_slvl = calc_flgt_envel(naca2415(1),rho_sealvl,'Sea Level')
n_allow_ceil = calc_flgt_envel(naca2415(2),rho_altceil,'Ceiling Altitude (14600 feet)')

calc_centroid_momentinertia;

% SEA LEVEL Load Distributions
% AT ALL CRITICAL CONDITIONS
% NOTE:
%  wx0 and wy0 here are defined from root to tip.
%  wy and wx below are defined from tip to root.
nz = 500;
for ii = 1:length(n_allow_slvl.n)
    if ~isnan(n_allow_slvl.n(ii))
        [load_slvl(ii).wx0 load_slvl(ii).wy0] = ...
                          calc_wxwy('Sea Level',...
                                    n_allow_slvl.name(ii),...
                                    n_allow_slvl.n(ii),...
                                    rho_sealvl,...
                                    n_allow_slvl.V(ii),...
                                    n_allow_slvl.AoA(ii),...
                                    n_allow_slvl.Cd(ii),nz);
        load_slvl(ii).wy = zeros(1,nz+1);
        load_slvl(ii).wx = zeros(1,nz+1);
        for i = 1:nz+1
            load_slvl(ii).wy(i) = load_slvl(ii).wx0(nz+2-i);
            load_slvl(ii).wx(i) = load_slvl(ii).wy0(nz+2-i);
        end
        
        shear_moment_slvl(ii) = calc_shear_moments('Sea Level',...
                                    n_allow_slvl.name(ii),...
                                    b, nz, load_slvl(ii).wx,...
                                    load_slvl(ii).wy,...
                                    load_slvl(ii).wx0,...
                                    load_slvl(ii).wy0);
        
    end
end

% CEILING ALTITUDE Load Distributions
% AT ALL CRITICAL CONDITIONS
for ii = 1:length(n_allow_ceil.n)
    if ~isnan(n_allow_ceil.n(ii))
        [load_ceil(ii).wx0 load_ceil(ii).wy0] = ...
                          calc_wxwy('Ceiling Altitude',...
                                    n_allow_ceil.name(ii),...
                                    n_allow_ceil.n(ii),...
                                    rho_altceil,...
                                    n_allow_ceil.V(ii),...
                                    n_allow_ceil.AoA(ii),...
                                    n_allow_ceil.Cd(ii),nz);
        load_ceil(ii).wy = zeros(1,nz+1);
        load_ceil(ii).wx = zeros(1,nz+1);
        for i = 1:nz+1
            load_ceil(ii).wy(i) = load_ceil(ii).wx0(nz+2-i);
            load_ceil(ii).wx(i) = load_ceil(ii).wy0(nz+2-i);
        end
        
        shear_moment_ceil(ii) = calc_shear_moments('Ceiling Altitude',...
                                    n_allow_ceil.name(ii),...
                                    b, nz, load_ceil(ii).wx,...
                                    load_ceil(ii).wy,...
                                    load_ceil(ii).wx0,...
                                    load_ceil(ii).wy0);   
    end
end

% TODO: airfoil section properties
% A_cap = ;
% A_str = ;
% t_spar = ;
% t_skin = ;
% % locations of spars, spar caps and stringers (nose at the origin of the coordinate)
% x_spar0 = ;                       % front spar (2 cell beam)
% x_strU0 = [];                     % upper surface
% x_strL0 = [];                     % lower surface
% % new coordinate with origin at the centroid is used for the output below
% [c,Ixx,Iyy,Ixy,x,yU,yL,x_strU,x_strL,x_boomU,x_boomL,L_boomU,L_boomL,x_spar,h_spar,i_spar,dx] = airfoil_section(A_cap,A_str,t_spar,t_skin,x_spar0,x_strU0,x_strL0);

