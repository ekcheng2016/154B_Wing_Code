clear all; close all; clc;

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
for ii = 1:length(n_allow_slvl.n)
    if ~isnan(n_allow_slvl.n(ii))
        [wx_slvl wy_slvl] = calc_wxwy(n_allow_slvl.name(ii),...
                                    n_allow_slvl.n(ii),...
                                    rho_sealvl,...
                                    n_allow_slvl.V(ii),...
                                    n_allow_slvl.AoA(ii),...
                                    n_allow_slvl.Cd(ii),500);
    end
end