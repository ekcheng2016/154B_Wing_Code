% optimization.m

clear all; 
close all; 
clc;

% NUMBER OF MONTE CARLO ITERATIONS
N_ITERATIONS = 1;

PLOT_FLTENVELOPE = 0; % set 1 to plot and save flight envelope plots
PLOT_AIRFOIL = 1;     % set 1 to plot and save airfoil section plots
PLOT_LIFTCURVE = 0;   % set 1 to plot and save lift curve slope plot
PLOT_LOADS = 0;
PLOT_DISTRIBUTION = 1;
if PLOT_DISTRIBUTION
    fid = fopen([pwd '/Optimized/Shear_Flow_Figure/shear_flow_numbers.txt'],'w');
end

load_aircraft_parameters;
load_conversions;

Re_sealvl = calc_Re(rho_sealvl,c,v_maneuver,mu_sealvl);
Re_alceil = calc_Re(rho_altceil,c,v_maneuver,mu_altceil);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 CALCULATE FLIGHT ENVELOPE & LOADS                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
airfoil_to_wing;
n_allow_slvl = calc_flgt_envel(naca2415(1),rho_sealvl,'Sea Level',PLOT_FLTENVELOPE);
n_allow_ceil = calc_flgt_envel(naca2415(2),rho_altceil,'Ceiling Altitude (14600 feet)',PLOT_FLTENVELOPE);
% LOADS ARE CALCULATED AT BOTH SEA LEVEL AND CEILING ALTITUDE
nz = 500;
calc_loads;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% SIMPLIFIED MODEL
% calc_centroid_momentinertia;

% ACTUAL AIRFOIL SECTION
load_base_wing;

airf_geo = base_wing;
% % new coordinate with origin at the centroid is used for the output below
[Cx,Cy,Ixx,Iyy,Ixy,I_str,airf_geo] = airfoil_section(c,airf_geo);

disp('Begin Monte Carlo Optimization Process.');
BASE_WING_WGT = 420;
NUM_SUCCESS = 0;
for kk = 1:N_ITERATIONS

if mod(kk,2) == 0
    fprintf('.')
end
if mod(kk,10) == 1 && kk > 1
    fprintf(' %d iterations. Number of Successes: %d \n',kk-1,NUM_SUCCESS);
end

[Cx,Cy,Ixx,Iyy,Ixy,I_str,airf_geo] = airfoil_section(c,airf_geo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           CALCULATE FORCES, MOMENTS, STRESSES @ SEA LEVEL           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE:
%  wx0 and wy0 here are defined from root to tip.
%  wy and wx below are defined from tip to root.
for ii = 1:length(n_allow_slvl.n)
    if ~isnan(n_allow_slvl.n(ii))
        % DETERMINE SHEARS AND MOMENTS          
        [shear_slvl(ii) moment_slvl(ii)] = calc_shear_moments(b, nz,...
                                    load_slvl(ii).wx,load_slvl(ii).wy,...
                                    load_slvl(ii).wx0,load_slvl(ii).wy0);
        
        % calculate shear flow
        tau_sz_slvl(ii) = calc_shear_flow(Ixx,Iyy,Ixy,airf_geo,...
                                moment_slvl(ii).Mx0, moment_slvl(ii).My0,...
                                shear_slvl(ii).Sx0, shear_slvl(ii).Sy0,...
                                load_slvl(ii).M0,c,Cx,Cy,ii);           
                  
        % DETERMINE DEFLECTIONS
        deflection_slvl(ii) = calc_deflections(b, Ixx, Iyy, Ixy,nz,...
                                    moment_slvl(ii).Mx0,moment_slvl(ii).My0,...
                                    load_slvl(ii).wx0,load_slvl(ii).wy0);
                              
        
        % SIGMA_ZZ AT THE ROOT
        [sigma_zz_slvl(ii)] = calc_sigmazz(Ixx,Iyy,Ixy,...
                                moment_slvl(ii).Mx0(1),moment_slvl(ii).My0(1),...
                                airf_geo.x,airf_geo.yU,airf_geo.x,airf_geo.yL);

    end
end
    
sigma_zz_MAX_slvl_val = max([sigma_zz_slvl(1:end).max])/1e6;
sigma_zz_MIN_slvl_val = min([sigma_zz_slvl(1:end).min])/1e6;
sigma_zz_MAX_slvl_ind = find([sigma_zz_slvl(1:end).max]/1e6 == sigma_zz_MAX_slvl_val);
sigma_zz_MIN_slvl_ind = find([sigma_zz_slvl(1:end).min]/1e6 == sigma_zz_MIN_slvl_val);
tau_sz_MAX_slvl_val = max([tau_sz_slvl(1:end).max])/1e6;
tau_sz_MAX_slvl_ind = find([tau_sz_slvl(1:end).max]/1e6 == tau_sz_MAX_slvl_val);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           LOAD DISTRIBUTIONS @ CEILING   (all critical pts)         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:length(n_allow_ceil.n)
    if ~isnan(n_allow_ceil.n(ii))
        % DETERMINE SHEARS AND MOMENTS          
        [shear_ceil(ii) moment_ceil(ii)] = calc_shear_moments(b, nz,...
                                    load_ceil(ii).wx,load_ceil(ii).wy,...
                                    load_ceil(ii).wx0,load_ceil(ii).wy0);
        
        % calculate shear flow
        tau_sz_ceil(ii) = calc_shear_flow(Ixx,Iyy,Ixy,airf_geo,...
                                moment_ceil(ii).Mx0, moment_ceil(ii).My0,...
                                shear_ceil(ii).Sx0, shear_ceil(ii).Sy0,...
                                load_ceil(ii).M0,c,Cx,Cy,0);                     % PLOT
        
        % DETERMINE DEFLECTIONS
        deflection_ceil(ii) = calc_deflections(b, Ixx, Iyy, Ixy,nz,...
                                    moment_ceil(ii).Mx0,moment_ceil(ii).My0,...
                                    load_ceil(ii).wx0,load_ceil(ii).wy0);
        
        % SIGMA_ZZ AT THE ROOT
        [sigma_zz_ceil(ii)] = calc_sigmazz(Ixx,Iyy,Ixy,...
                            moment_ceil(ii).Mx0(1),moment_ceil(ii).My0(1),...
                            airf_geo.x,airf_geo.yU,airf_geo.x,airf_geo.yL);
                  
    end
end

sigma_zz_MAX_ceil_val = max([sigma_zz_ceil(1:end).max])/1e6;
sigma_zz_MIN_ceil_val = min([sigma_zz_ceil(1:end).min])/1e6;
sigma_zz_MAX_ceil_ind = find([sigma_zz_ceil(1:end).max]/1e6 == sigma_zz_MAX_ceil_val);
sigma_zz_MIN_ceil_ind = find([sigma_zz_ceil(1:end).min]/1e6 == sigma_zz_MIN_ceil_val);
tau_sz_MAX_ceil_val = max([tau_sz_ceil(1:end).max])/1e6;
tau_sz_MAX_ceil_ind = find([tau_sz_ceil(1:end).max]/1e6 == tau_sz_MAX_ceil_val);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                        Buckling & Fatigue                         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

buckling = calc_buckling(I_str,max([sigma_zz_MAX_ceil_val sigma_zz_MAX_slvl_val]),...
                        min([sigma_zz_MIN_ceil_val sigma_zz_MIN_slvl_val]),airf_geo.A_str,airf_geo.t_skin,airf_geo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                        Von Mises Failure                          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigma_eq = von_mises([sigma_zz_slvl(:).max],[sigma_zz_ceil(:).max],...
                    [tau_sz_slvl(:).max],[tau_sz_ceil(:).max]);

if (((sigma_eq.val/1e6)*1.5 >= sigma_yield) && ...
   ((1.5*max(abs([sigma_zz_MIN_ceil_val sigma_zz_MIN_slvl_val]))) >= buckling.sigma_crit))
    calc_random_wing;
    fields = {'x','xU','y_skinU','yU','xL','y_skinL','yL','L_skinU',...
          'L_skinL','x_strU','x_strL','x_spar','L_boomU','L_boomL','h_spar',...
          'i_spar','dx','i_strU','i_strL'};
    airf_geo = rmfield(airf_geo,fields);
    continue;
end

weight_wing = calc_weight_wing(airf_geo,b,rho_material);

if weight_wing.total < BASE_WING_WGT
   NUM_SUCCESS = NUM_SUCCESS + 1;
   WING_SUCCESS(NUM_SUCCESS).weight_wing = weight_wing;
   WING_SUCCESS(NUM_SUCCESS).airf_geo = airf_geo;
   WING_SUCCESS(NUM_SUCCESS).Cx       = Cx;
   WING_SUCCESS(NUM_SUCCESS).Cy       = Cy;
   WING_SUCCESS(NUM_SUCCESS).Ixx      = Ixx;
   WING_SUCCESS(NUM_SUCCESS).Iyy      = Iyy;
   WING_SUCCESS(NUM_SUCCESS).Ixy      = Ixy;
   WING_SUCCESS(NUM_SUCCESS).shear_slvl = shear_slvl;
   WING_SUCCESS(NUM_SUCCESS).shear_ceil = shear_ceil;
   WING_SUCCESS(NUM_SUCCESS).moment_slvl = moment_slvl;
   WING_SUCCESS(NUM_SUCCESS).moment_ceil = moment_ceil;
   WING_SUCCESS(NUM_SUCCESS).deflection_slvl = deflection_slvl;
   WING_SUCCESS(NUM_SUCCESS).deflection_ceil = deflection_ceil;
   WING_SUCCESS(NUM_SUCCESS).sigma_zz_slvl = sigma_zz_slvl;
   WING_SUCCESS(NUM_SUCCESS).sigma_zz_ceil = sigma_zz_ceil;
   WING_SUCCESS(NUM_SUCCESS).tau_sz_slvl = tau_sz_slvl;
   WING_SUCCESS(NUM_SUCCESS).tau_sz_ceil = tau_sz_ceil;
   WING_SUCCESS(NUM_SUCCESS).sigma_zz_MAX_slvl_ind = sigma_zz_MAX_slvl_ind;
   WING_SUCCESS(NUM_SUCCESS).sigma_zz_MAX_slvl_val = sigma_zz_MAX_slvl_val;
   WING_SUCCESS(NUM_SUCCESS).sigma_zz_MAX_ceil_ind = sigma_zz_MAX_ceil_ind;
   WING_SUCCESS(NUM_SUCCESS).sigma_zz_MAX_ceil_val = sigma_zz_MAX_ceil_val;
   WING_SUCCESS(NUM_SUCCESS).sigma_zz_MIN_slvl_ind = sigma_zz_MIN_slvl_ind;
   WING_SUCCESS(NUM_SUCCESS).sigma_zz_MIN_slvl_val = sigma_zz_MIN_slvl_val;
   WING_SUCCESS(NUM_SUCCESS).sigma_zz_MIN_ceil_ind = sigma_zz_MIN_ceil_ind;
   WING_SUCCESS(NUM_SUCCESS).sigma_zz_MIN_ceil_val = sigma_zz_MIN_ceil_val;
   WING_SUCCESS(NUM_SUCCESS).tau_sz_MAX_slvl_ind   = tau_sz_MAX_slvl_ind;
   WING_SUCCESS(NUM_SUCCESS).tau_sz_MAX_slvl_val   = tau_sz_MAX_slvl_val;
   WING_SUCCESS(NUM_SUCCESS).tau_sz_MAX_ceil_ind   = tau_sz_MAX_ceil_ind;
   WING_SUCCESS(NUM_SUCCESS).tau_sz_MAX_ceil_val   = tau_sz_MAX_ceil_val;
   WING_SUCCESS(NUM_SUCCESS).sigma_eq = sigma_eq;
   WING_SUCCESS(NUM_SUCCESS).buckling = buckling;
end

fields = {'x','xU','y_skinU','yU','xL','y_skinL','yL','L_skinU',...
          'L_skinL','x_strU','x_strL','x_spar','L_boomU','L_boomL','h_spar',...
          'i_spar','dx','i_strU','i_strL'};
airf_geo = rmfield(airf_geo,fields);
calc_random_wing;

end
display('End Monte Carlo Optimization Process');

[min_wgt ind] = min(arrayfun(@(x) x.weight_wing.total, WING_SUCCESS));
disp(strjoin(['Max sigma_zz at Sea Level : ' num2str(WING_SUCCESS(ind).sigma_zz_MAX_slvl_val) ...
    'MPa, occurs at : ' n_allow_ceil.name(WING_SUCCESS(ind).sigma_zz_MAX_slvl_ind)]))
disp(strjoin(['Min sigma_zz at Sea Level : ' num2str(WING_SUCCESS(ind).sigma_zz_MIN_slvl_val) ...
    'MPa, occurs at : ' n_allow_ceil.name(WING_SUCCESS(ind).sigma_zz_MIN_slvl_ind)]))
disp(sprintf(strjoin(['Max tau_sz at Sea Level : ' num2str(WING_SUCCESS(ind).tau_sz_MAX_slvl_val) ...
    'MPa, occurs at : ' n_allow_slvl.name(WING_SUCCESS(ind).tau_sz_MAX_slvl_ind) '\n'])))

disp(strjoin(['Max sigma_zz at Ceiling : ' num2str(WING_SUCCESS(ind).sigma_zz_MAX_ceil_val) ...
    'MPa, occurs at : ' n_allow_ceil.name(WING_SUCCESS(ind).sigma_zz_MAX_ceil_ind)]))
disp(strjoin(['Min sigma_zz at Ceiling : ' num2str(WING_SUCCESS(ind).sigma_zz_MIN_ceil_val) ...
    'MPa, occurs at : ' n_allow_ceil.name(WING_SUCCESS(ind).sigma_zz_MIN_ceil_ind)]))
disp(sprintf(strjoin(['Max tau_sz at Ceiling : ' num2str(WING_SUCCESS(ind).tau_sz_MAX_ceil_val) ...
    'MPa, occurs at : ' n_allow_slvl.name(WING_SUCCESS(ind).tau_sz_MAX_ceil_ind) '\n'])))

disp(sprintf(strjoin(['Buckling Critical Stress: ' num2str(WING_SUCCESS(ind).buckling.sigma_crit) 'MPa \n'])));

disp(sprintf(strjoin(['Von Mises Equivalent Stress : ' num2str(WING_SUCCESS(ind).sigma_eq.val/1e6) ...
    'MPa, occurs at : ' n_allow_slvl.name(WING_SUCCESS(ind).sigma_eq.ind) ', Flight Condition: ' sigma_eq.fgt_cond '\n'])));
                
disp(strjoin(['Based on the current configuration, the wing weighs : ' ...
             num2str(WING_SUCCESS(ind).weight_wing.total) 'kg']));
disp(strjoin(['The wing weight ' num2str(100*WING_SUCCESS(ind).weight_wing.total/(mass_emp)) ...
              '% of the entire aircraft empty weight.']));

          
plot_airfoil(WING_SUCCESS(ind).airf_geo,WING_SUCCESS(ind).Cx,WING_SUCCESS(ind).Cy);

if PLOT_DISTRIBUTION
plot_distributions('Sea Level',fid,b,c,WING_SUCCESS(ind).Cx,WING_SUCCESS(ind).airf_geo,n_allow_slvl,WING_SUCCESS(ind).shear_slvl,WING_SUCCESS(ind).moment_slvl,...
    WING_SUCCESS(ind).deflection_slvl,WING_SUCCESS(ind).sigma_zz_slvl(sigma_zz_MAX_slvl_ind),WING_SUCCESS(ind).tau_sz_slvl);

plot_distributions('Ceiling Altitude',fid,b,c,WING_SUCCESS(ind).Cx,WING_SUCCESS(ind).airf_geo,n_allow_ceil,WING_SUCCESS(ind).shear_ceil,WING_SUCCESS(ind).moment_ceil,...
    WING_SUCCESS(ind).deflection_ceil,WING_SUCCESS(ind).sigma_zz_ceil(sigma_zz_MAX_ceil_ind),WING_SUCCESS(ind).tau_sz_ceil);
end

fprintf(fid,'\n');
fprintf(fid,strjoin(['Max sigma_zz at Sea Level : ' num2str(WING_SUCCESS(ind).sigma_zz_MAX_slvl_val) ...
    'MPa, occurs at : ' n_allow_ceil.name(WING_SUCCESS(ind).sigma_zz_MAX_slvl_ind)]));
fprintf(fid,'\n');
fprintf(fid,strjoin(['Min sigma_zz at Sea Level : ' num2str(WING_SUCCESS(ind).sigma_zz_MIN_slvl_val) ...
    'MPa, occurs at : ' n_allow_ceil.name(WING_SUCCESS(ind).sigma_zz_MIN_slvl_ind)]));
fprintf(fid,'\n');
fprintf(fid,sprintf(strjoin(['Max tau_sz at Sea Level : ' num2str(WING_SUCCESS(ind).tau_sz_MAX_slvl_val) ...
    'MPa, occurs at : ' n_allow_slvl.name(WING_SUCCESS(ind).tau_sz_MAX_slvl_ind) '\n'])));

fprintf(fid,'\n');
fprintf(fid,strjoin(['Max sigma_zz at Ceiling : ' num2str(WING_SUCCESS(ind).sigma_zz_MAX_ceil_val) ...
    'MPa, occurs at : ' n_allow_ceil.name(WING_SUCCESS(ind).sigma_zz_MAX_ceil_ind)]));
fprintf(fid,'\n');
fprintf(fid,strjoin(['Min sigma_zz at Ceiling : ' num2str(WING_SUCCESS(ind).sigma_zz_MIN_ceil_val) ...
    'MPa, occurs at : ' n_allow_ceil.name(WING_SUCCESS(ind).sigma_zz_MIN_ceil_ind)]));
fprintf(fid,'\n');
fprintf(fid,sprintf(strjoin(['Max tau_sz at Ceiling : ' num2str(WING_SUCCESS(ind).tau_sz_MAX_ceil_val) ...
    'MPa, occurs at : ' n_allow_slvl.name(WING_SUCCESS(ind).tau_sz_MAX_ceil_ind) '\n'])));

fprintf(fid,'\n');
fprintf(fid,sprintf(strjoin(['Buckling Critical Stress: ' num2str(WING_SUCCESS(ind).buckling.sigma_crit) 'MPa \n'])));

fprintf(fid,'\n');
fprintf(fid,sprintf(strjoin(['Von Mises Equivalent Stress : ' num2str(WING_SUCCESS(ind).sigma_eq.val/1e6) ...
    'MPa, occurs at : ' n_allow_slvl.name(WING_SUCCESS(ind).sigma_eq.ind) ', Flight Condition: ' sigma_eq.fgt_cond '\n'])));
                
fprintf(fid,'\n');
fprintf(fid,strjoin(['Based on the current configuration, the wing weighs : ' ...
             num2str(WING_SUCCESS(ind).weight_wing.total) 'kg']));
fprintf(fid,'\n');
fprintf(fid,strjoin(['The wing weight ' num2str(100*WING_SUCCESS(ind).weight_wing.total/(mass_emp)) ...
              '% of the entire aircraft empty weight.']));
