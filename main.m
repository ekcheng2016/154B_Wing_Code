clear all; 
close all; 
clc;

PLOT_FLTENVELOPE = 0; % set 1 to plot and save flight envelope plots
PLOT_AIRFOIL = 0;     % set 1 to plot and save airfoil section plots
PLOT_LIFTCURVE = 0;   % set 1 to plot and save lift curve slope plot
PLOT_SHEAR_FLOW = 0;
PLOT_DEFLECTION = 0;    % set 1 to plot and save all load/shear plots
PLOT_SHEAR = 0;
PLOT_MOMENT = 0;
PLOT_LOADS = 0;
PLOT_SIGMAZZ = 0;
if PLOT_SHEAR_FLOW
    fid = fopen([pwd '/Shear_Flow_Figure/shear_flow_numbers.txt'],'w');
end

clrstring = 'bgkrc';

load_aircraft_parameters;
load_conversions;

Re_sealvl = calc_Re(rho_sealvl,c,v_maneuver,mu_sealvl);
Re_alceil = calc_Re(rho_altceil,c,v_maneuver,mu_altceil);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 CALCULATE FLIGHT ENVELOPE                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
airfoil_to_wing;
n_allow_slvl = calc_flgt_envel(naca2415(1),rho_sealvl,'Sea Level',PLOT_FLTENVELOPE);
n_allow_ceil = calc_flgt_envel(naca2415(2),rho_altceil,'Ceiling Altitude (14600 feet)',PLOT_FLTENVELOPE);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% SIMPLIFIED MODEL
% calc_centroid_momentinertia;

% ACTUAL AIRFOIL SECTION
airf_geo.A_cap = 7e-4;   % m^2
airf_geo.A_str = 3.5e-5;   % m^2
airf_geo.t_spar = 0.0025;   % m
airf_geo.t_skin = 0.001016;  % m
% % locations of spars, spar caps and stringers (nose at the origin of the coordinate)
airf_geo.x_spar0 = 0.25*c;                 % front spar (2 cell beam)
airf_geo.x_strU0 = [0 0.05 0.15 0.35 0.55 0.65]*c; % upper surface
airf_geo.x_strL0 = [0 0.05 0.15 0.35 0.55 0.65]*c; % lower surface
% % new coordinate with origin at the centroid is used for the output below
[Cx,Cy,Ixx,Iyy,Ixy,I_str,airf_geo] = ...
    airfoil_section(c,airf_geo);
                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           LOAD DISTRIBUTIONS @ SEA LEVEL (all critical pts)         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE:
%  wx0 and wy0 here are defined from root to tip.
%  wy and wx below are defined from tip to root.
nz = 500;
for ii = 1:length(n_allow_slvl.n)
    if ~isnan(n_allow_slvl.n(ii))
        % DETERMINE LOAD DISTRIBUTION
        [load_slvl(ii)] = calc_wxwy(n_allow_slvl.n(ii),rho_sealvl,...
                                    n_allow_slvl.V(ii),n_allow_slvl.AoA(ii),...
                                    n_allow_slvl.Cd(ii),n_allow_slvl.CM(ii),nz);
        
        %PLOT DISTRIBUTIONS
        if PLOT_LOADS
        lift_ellip_fig = figure(100);
        hold on; box on; grid on;
        lef(ii) = plot(load_slvl(ii).z,load_slvl(ii).l_ellip,'Color',clrstring(ii),'linewidth',2);
                  plot(-load_slvl(ii).z,load_slvl(ii).l_ellip,'Color',clrstring(ii),'linewidth',2);
        
        lift_rect_fig = figure(101);
        hold on; box on; grid on;
        lrf(ii) = plot(load_slvl(ii).z,load_slvl(ii).l_rect,'Color',clrstring(ii),'linewidth',2);
                  plot(-load_slvl(ii).z,load_slvl(ii).l_rect,'Color',clrstring(ii),'linewidth',2);
        
        lift_fig = figure(102);
        hold on; box on; grid on;
        lf(ii) = plot(load_slvl(ii).z,load_slvl(ii).l,'Color',clrstring(ii),'linewidth',2);
                 plot(-load_slvl(ii).z,load_slvl(ii).l,'Color',clrstring(ii),'linewidth',2);
        
        drag_fig = figure(103);
        hold on; box on; grid on;
        df(ii) = plot(load_slvl(ii).z,load_slvl(ii).d,'Color',clrstring(ii),'linewidth',2);
                 plot(-load_slvl(ii).z,load_slvl(ii).d,'Color',clrstring(ii),'linewidth',2);
        
        wx_fig = figure(104);
        hold on; box on; grid on;
        wxf(ii) = plot(load_slvl(ii).z,load_slvl(ii).wx0,'Color',clrstring(ii),'linewidth',2);
                  plot(-load_slvl(ii).z,load_slvl(ii).wx0,'Color',clrstring(ii),'linewidth',2);
        
        wy_fig = figure(105);
        hold on; box on; grid on;
        wyf(ii) = plot(load_slvl(ii).z,load_slvl(ii).wy0,'Color',clrstring(ii),'linewidth',2);
                  plot(-load_slvl(ii).z,load_slvl(ii).wy0,'Color',clrstring(ii),'linewidth',2);
        end
                  
        % DETERMINE SHEARS AND MOMENTS          
        [shear_slvl(ii) moment_slvl(ii)] = calc_shear_moments(b, nz,...
                                    load_slvl(ii).wx,load_slvl(ii).wy,...
                                    load_slvl(ii).wx0,load_slvl(ii).wy0);
        
        % calculate shear flow
        tau_sz_slvl(ii) = calc_shear_flow(Ixx,Iyy,Ixy,airf_geo,...
                                moment_slvl(ii).Mx0, moment_slvl(ii).My0,...
                                shear_slvl(ii).Sx0, shear_slvl(ii).Sy0,...
                                load_slvl(ii).M0,c,Cx,Cy,ii);             
        if PLOT_SHEAR_FLOW
            sf_fig = figure(113);
            hold on; box on; grid on;
            sff(ii) = plot(tau_sz_slvl(ii).qb/1e3,'color',clrstring(ii),'linewidth',2);
            tau_fig = figure(114);
            hold on; box on; grid on;
            tauf(ii) = plot(tau_sz_slvl(ii).skin/1e6,'color',clrstring(ii),'linewidth',2);
            fprintf(fid,'Flight Condition: Sea Level\n');
            fprintf(fid,'Loading Condition: %20s\n',char(n_allow_slvl.name(ii)));
            fprintf(fid,'q01 = %5.4f N/m\n',tau_sz_slvl(ii).q01);
            fprintf(fid,'q02 = %5.4f N/m\n',tau_sz_slvl(ii).q02);
            fprintf(fid,'tau_spars = %5.4f MPa, %5.4f MPa\n',tau_sz_slvl(ii).spar/1e6);
            fprintf(fid,'Max Shear Stress = %5.4f MPa\n',tau_sz_slvl(ii).max/1e6);
            fprintf(fid,'Last Open Shear Flow Value = %5.4f N/m\n',tau_sz_slvl(ii).qb(end));
            fprintf(fid,'Average Open Shear Flow Value = %5.4f N/m\n',mean(tau_sz_slvl(ii).qb));            
            fprintf(fid,'Last Open Shear Flow/Avg Shear Flow = %5.4f %% \n\n',(tau_sz_slvl(ii).qb(end)/mean(tau_sz_slvl(ii).qb))*100);
        end
        
        % PLOT
        if PLOT_SHEAR
        sx_fig = figure(106);
        hold on; box on; grid on;
        sxf(ii) = plot(shear_slvl(ii).z,shear_slvl(ii).Sx0/1e3,'Color',clrstring(ii),'linewidth',2);
                  plot(-shear_slvl(ii).z,shear_slvl(ii).Sx0/1e3,'Color',clrstring(ii),'linewidth',2);
        
        sy_fig = figure(107);
        hold on; box on; grid on;
        syf(ii) = plot(shear_slvl(ii).z,shear_slvl(ii).Sy0/1e3,'Color',clrstring(ii),'linewidth',2);
                  plot(-shear_slvl(ii).z,shear_slvl(ii).Sy0/1e3,'Color',clrstring(ii),'linewidth',2);
        end
        
        if PLOT_MOMENT
        mx_fig = figure(108);
        hold on; box on; grid on;
        mxf(ii) = plot(moment_slvl(ii).z,moment_slvl(ii).Mx0/1e3,'Color',clrstring(ii),'linewidth',2);
                  plot(-moment_slvl(ii).z,moment_slvl(ii).Mx0/1e3,'Color',clrstring(ii),'linewidth',2);
                  
        my_fig = figure(109);
        hold on; box on; grid on;
        myf(ii) = plot(moment_slvl(ii).z,moment_slvl(ii).My0/1e3,'Color',clrstring(ii),'linewidth',2);
                  plot(-moment_slvl(ii).z,moment_slvl(ii).My0/1e3,'Color',clrstring(ii),'linewidth',2);        
        end
                  
        % DETERMINE DEFLECTIONS
        deflection_slvl(ii) = calc_deflections(b, Ixx, Iyy, Ixy,nz,...
                                    moment_slvl(ii).Mx0,moment_slvl(ii).My0,...
                                    load_slvl(ii).wx0,load_slvl(ii).wy0);
                                
        if PLOT_DEFLECTION
        u_fig = figure(110);
        hold on; box on; grid on;
        uf(ii) = plot(deflection_slvl(ii).z,deflection_slvl(ii).u*1000,'Color',clrstring(ii),'linewidth',2);
                 plot(-deflection_slvl(ii).z,deflection_slvl(ii).u*1000,'Color',clrstring(ii),'linewidth',2);
        
        v_fig = figure(111);
        hold on; box on; grid on;
        vf(ii) = plot(deflection_slvl(ii).z,deflection_slvl(ii).v*1000,'Color',clrstring(ii),'linewidth',2);
                 plot(-deflection_slvl(ii).z,deflection_slvl(ii).v*1000,'Color',clrstring(ii),'linewidth',2);
        end
        
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
disp(strjoin(['Max sigma_zz at Sea Level : ' num2str(sigma_zz_MAX_slvl_val) ...
    'MPa, occurs at : ' n_allow_ceil.name(sigma_zz_MAX_slvl_ind)]))
disp(strjoin(['Min sigma_zz at Sea Level : ' num2str(sigma_zz_MIN_slvl_val) ...
    'MPa, occurs at : ' n_allow_ceil.name(sigma_zz_MIN_slvl_ind)]))
tau_sz_MAX_slvl_val = max([tau_sz_slvl(1:end).max])/1e6;
tau_sz_MAX_slvl_ind = find([tau_sz_slvl(1:end).max]/1e6 == tau_sz_MAX_slvl_val);
disp(strjoin(['Max tau_sz at Sea Level : ' num2str(tau_sz_MAX_slvl_val) ...
    'MPa, occurs at : ' n_allow_slvl.name(tau_sz_MAX_slvl_ind)]))


if PLOT_LOADS
figure(100);    xlabel('Span (m)');     ylabel('Elliptical Lift Distribution (N/m)');
                legend(lef(:),n_allow_slvl(:).name);    xlim([-b/2 b/2]);
                title('Sea Level'); pos = get(lift_ellip_fig,'position');
                set(lift_ellip_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(lift_ellip_fig,[pwd '/Load_Distribution_Figures/SLVL_Lift_Elliptical'],'-djpeg','-r300');
                
figure(101);    xlabel('Span (m)');     ylabel('Rectangular Lift Distribution (N/m)');
                legend(lrf(:),n_allow_slvl(:).name);    xlim([-b/2 b/2]);
                title('Sea Level'); pos = get(lift_rect_fig,'position');
                set(lift_rect_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(lift_rect_fig,[pwd '/Load_Distribution_Figures/SLVL_Lift_Rectangular'],'-djpeg','-r300');
                
figure(102);    xlabel('Span (m)');     ylabel('Lift Distribution (N/m)');
                legend(lf(:),n_allow_slvl(:).name);    xlim([-b/2 b/2]);
                title('Sea Level'); pos = get(lift_fig,'position');
                set(lift_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(lift_fig,[pwd '/Load_Distribution_Figures/SLVL_Lift_TTL'],'-djpeg','-r300');
                
figure(103);    xlabel('Span (m)');     ylabel('Drag Distribution (N/m)');
                legend(df(:),n_allow_slvl(:).name);    xlim([-b/2 b/2]);
                title('Sea Level'); pos = get(drag_fig,'position');
                set(drag_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(drag_fig,[pwd '/Load_Distribution_Figures/SLVL_Drag'],'-djpeg','-r300');
                
figure(104);    xlabel('Span (m)');     ylabel('w_x Distribution (N/m)');
                legend(wxf(:),n_allow_slvl(:).name);    xlim([-b/2 b/2]);
                title('Sea Level'); pos = get(wx_fig,'position');
                set(wx_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(wx_fig,[pwd '/Load_Distribution_Figures/SLVL_wx'],'-djpeg','-r300');
                
figure(105);    xlabel('Span (m)');     ylabel('w_y Distribution (N/m)');
                legend(wyf(:),n_allow_slvl(:).name);    xlim([-b/2 b/2]);
                title('Sea Level'); pos = get(wy_fig,'position');
                set(wy_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(wy_fig,[pwd '/Load_Distribution_Figures/SLVL_wy'],'-djpeg','-r300');
end

if PLOT_SHEAR
figure(106);    xlabel('Span (m)');     ylabel('S_x Distribution (kN)');
                legend(sxf(:),n_allow_slvl(:).name);    xlim([-b/2 b/2]);
                title('Sea Level'); pos = get(sx_fig,'position');
                set(sx_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(sx_fig,[pwd '/Load_Distribution_Figures/SLVL_Sx'],'-djpeg','-r300');
                
figure(107);    xlabel('Span (m)');     ylabel('S_y Distribution (kN)');
                legend(syf(:),n_allow_slvl(:).name);    xlim([-b/2 b/2]);
                title('Sea Level'); pos = get(sy_fig,'position');
                set(sy_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(sy_fig,[pwd '/Load_Distribution_Figures/SLVL_Sy'],'-djpeg','-r300');
end

if PLOT_MOMENT
figure(108);    xlabel('Span (m)');     ylabel('M_x Distribution (kNm)');
                legend(mxf(:),n_allow_slvl(:).name);    xlim([-b/2 b/2]);
                title('Sea Level'); pos = get(mx_fig,'position');
                set(mx_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(mx_fig,[pwd '/Load_Distribution_Figures/SLVL_Mx'],'-djpeg','-r300');
                
figure(109);    xlabel('Span (m)');     ylabel('M_y Distribution (kNm)');
                legend(myf(:),n_allow_slvl(:).name);    xlim([-b/2 b/2]);
                title('Sea Level'); pos = get(my_fig,'position');
                set(my_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(my_fig,[pwd '/Load_Distribution_Figures/SLVL_My'],'-djpeg','-r300');
end

if PLOT_DEFLECTION
figure(110);    xlabel('Span (m)');     ylabel('u deflection (mm)');
                legend(uf(:),n_allow_slvl(:).name);    xlim([-b/2 b/2]);
                title('Sea Level'); pos = get(u_fig,'position');
                set(u_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(u_fig,[pwd '/Deflection_Figures/SLVL_u'],'-djpeg','-r300');
                
figure(111);    xlabel('Span (m)');     ylabel('v deflection (mm)');
                legend(vf(:),n_allow_slvl(:).name);    xlim([-b/2 b/2]);
                title('Sea Level'); pos = get(v_fig,'position');
                set(v_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(v_fig,[pwd '/Deflection_Figures/SLVL_v'],'-djpeg','-r300');                
end

if PLOT_SHEAR_FLOW
figure(113);    sff(end+1) = plot([0 length(tau_sz_slvl(1).qb)],[0 0],'r','linewidth',2);
                xlabel('Node (counterclockwise)');  ylabel('Shear Flow (kN/m)');
                legend(sff(:),[n_allow_slvl(:).name 'Zero']);
                title('Sea Level'); pos = get(sf_fig,'position');
                set(sf_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(sf_fig,[pwd '/Shear_Flow_Figure/SLVL_sf_open'],'-djpeg','-r300');
figure(114);    xlabel('Node (counterclockwise)');  ylabel('Shear Force (\tau_{sz}) (MPa)');
                legend(tauf(:),n_allow_slvl(:).name);
                title('Sea Level'); pos = get(tau_fig,'position');
                set(tau_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(tau_fig,[pwd '/Shear_Flow_Figure/SLVL_tau_sz'],'-djpeg','-r300');
end

if PLOT_SIGMAZZ
sigmazz_fig = figure(112);
hold on; box on; grid on;
plot(airf_geo.x,sigma_zz_slvl(sigma_zz_MAX_slvl_ind).upper/1e6,'b','linewidth',2);
plot(airf_geo.x,sigma_zz_slvl(sigma_zz_MAX_slvl_ind).lower/1e6,'g','linewidth',2);
xlabel('Chord (m)');    ylabel('Direct Stress (\sigma_{zz}) (MPa)');
legend('Upper Surface','Lower Surface');    xlim([0-Cx c-Cx]);
title('Sea Level'); pos = get(sigmazz_fig,'position');
set(sigmazz_fig,'position',[pos(1:2) pos(3:4)*1.5]);
print(sigmazz_fig,[pwd '/Stress_Figure/SLVL_DirectStress'],'-djpeg','-r300');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           LOAD DISTRIBUTIONS @ CEILING   (all critical pts)         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:length(n_allow_ceil.n)
    if ~isnan(n_allow_ceil.n(ii))
        % DETERMINE LOAD DISTRIBUTION
        [load_ceil(ii)] = calc_wxwy(n_allow_ceil.n(ii),rho_altceil,...
                                    n_allow_ceil.V(ii),n_allow_ceil.AoA(ii),...
                                    n_allow_ceil.Cd(ii),n_allow_ceil.CM(ii),nz);
        
        %PLOT DISTRIBUTIONS
        if PLOT_LOADS
        lift_ellip_fig = figure(200);
        hold on; box on; grid on;
        lef(ii) = plot(load_ceil(ii).z,load_ceil(ii).l_ellip,'Color',clrstring(ii),'linewidth',2);
                  plot(-load_ceil(ii).z,load_ceil(ii).l_ellip,'Color',clrstring(ii),'linewidth',2);
        
        lift_rect_fig = figure(201);
        hold on; box on; grid on;
        lrf(ii) = plot(load_ceil(ii).z,load_ceil(ii).l_rect,'Color',clrstring(ii),'linewidth',2);
                  plot(-load_ceil(ii).z,load_ceil(ii).l_rect,'Color',clrstring(ii),'linewidth',2);
        
        lift_fig = figure(202);
        hold on; box on; grid on;
        lf(ii) = plot(load_ceil(ii).z,load_ceil(ii).l,'Color',clrstring(ii),'linewidth',2);
                 plot(-load_ceil(ii).z,load_ceil(ii).l,'Color',clrstring(ii),'linewidth',2);
        
        drag_fig = figure(203);
        hold on; box on; grid on;
        df(ii) = plot(load_ceil(ii).z,load_ceil(ii).d,'Color',clrstring(ii),'linewidth',2);
                 plot(-load_ceil(ii).z,load_ceil(ii).d,'Color',clrstring(ii),'linewidth',2);
        
        wx_fig = figure(204);
        hold on; box on; grid on;
        wxf(ii) = plot(load_ceil(ii).z,load_ceil(ii).wx0,'Color',clrstring(ii),'linewidth',2);
                  plot(-load_ceil(ii).z,load_ceil(ii).wx0,'Color',clrstring(ii),'linewidth',2);
        
        wy_fig = figure(205);
        hold on; box on; grid on;
        wyf(ii) = plot(load_ceil(ii).z,load_ceil(ii).wy0,'Color',clrstring(ii),'linewidth',2);
                  plot(-load_ceil(ii).z,load_ceil(ii).wy0,'Color',clrstring(ii),'linewidth',2);
        end
                  
        % DETERMINE SHEARS AND MOMENTS          
        [shear_ceil(ii) moment_ceil(ii)] = calc_shear_moments(b, nz,...
                                    load_ceil(ii).wx,load_ceil(ii).wy,...
                                    load_ceil(ii).wx0,load_ceil(ii).wy0);
        
        % calculate shear flow
        tau_sz_ceil(ii) = calc_shear_flow(Ixx,Iyy,Ixy,airf_geo,...
                                moment_ceil(ii).Mx0, moment_ceil(ii).My0,...
                                shear_ceil(ii).Sx0, shear_ceil(ii).Sy0,...
                                load_ceil(ii).M0,c,Cx,Cy,0);                     % PLOT
        
        if PLOT_SHEAR_FLOW
            sf_fig = figure(213);
            hold on; box on; grid on;
            sff(ii) = plot(tau_sz_ceil(ii).qb/1e3,'color',clrstring(ii),'linewidth',2);
            tau_fig = figure(214);
            hold on; box on; grid on;
            tauf(ii) = plot(tau_sz_ceil(ii).skin/1e6,'color',clrstring(ii),'linewidth',2);        
            fprintf(fid,'Flight Condition: Ceiling\n');
            fprintf(fid,'Loading Condition: %20s\n',char(n_allow_ceil.name(ii)));
            fprintf(fid,'q01 = %5.4f N/m\n',tau_sz_ceil(ii).q01);
            fprintf(fid,'q02 = %5.4f N/m\n',tau_sz_ceil(ii).q02);
            fprintf(fid,'tau_spars = %5.4f MPa, %5.4f MPa\n',tau_sz_ceil(ii).spar/1e6);
            fprintf(fid,'Max Shear Stress = %5.4f MPa\n',tau_sz_ceil(ii).max/1e6);
            fprintf(fid,'Last Open Shear Flow Value = %5.4f N/m\n',tau_sz_ceil(ii).qb(end));
            fprintf(fid,'Average Open Shear Flow Value = %5.4f N/m\n',mean(tau_sz_ceil(ii).qb));
            fprintf(fid,'Last Open Shear Flow/Avg Shear Flow = %5.4f %% \n\n',(tau_sz_ceil(ii).qb(end)/mean(tau_sz_ceil(ii).qb))*100);
        end
        
        if PLOT_SHEAR
        sx_fig = figure(206);
        hold on; box on; grid on;
        sxf(ii) = plot(shear_ceil(ii).z,shear_ceil(ii).Sx0/1e3,'Color',clrstring(ii),'linewidth',2);
                  plot(-shear_ceil(ii).z,shear_ceil(ii).Sx0/1e3,'Color',clrstring(ii),'linewidth',2);
        
        sy_fig = figure(207);
        hold on; box on; grid on;
        syf(ii) = plot(shear_ceil(ii).z,shear_ceil(ii).Sy0/1e3,'Color',clrstring(ii),'linewidth',2);
                  plot(-shear_ceil(ii).z,shear_ceil(ii).Sy0/1e3,'Color',clrstring(ii),'linewidth',2);
        end
        
        if PLOT_MOMENT
        mx_fig = figure(208);
        hold on; box on; grid on;
        mxf(ii) = plot(moment_ceil(ii).z,moment_ceil(ii).Mx0/1e3,'Color',clrstring(ii),'linewidth',2);
                  plot(-moment_ceil(ii).z,moment_ceil(ii).Mx0/1e3,'Color',clrstring(ii),'linewidth',2);
                  
        my_fig = figure(209);
        hold on; box on; grid on;
        myf(ii) = plot(moment_ceil(ii).z,moment_ceil(ii).My0/1e3,'Color',clrstring(ii),'linewidth',2);
                  plot(-moment_ceil(ii).z,moment_ceil(ii).My0/1e3,'Color',clrstring(ii),'linewidth',2);        
        end
        
        % DETERMINE DEFLECTIONS
        deflection_ceil(ii) = calc_deflections(b, Ixx, Iyy, Ixy,nz,...
                                    moment_ceil(ii).Mx0,moment_ceil(ii).My0,...
                                    load_ceil(ii).wx0,load_ceil(ii).wy0);
        
        if PLOT_DEFLECTION
        u_fig = figure(210);
        hold on; box on; grid on;
        uf(ii) = plot(deflection_ceil(ii).z,deflection_ceil(ii).u*1000,'Color',clrstring(ii),'linewidth',2);
                 plot(-deflection_ceil(ii).z,deflection_ceil(ii).u*1000,'Color',clrstring(ii),'linewidth',2);
        
        v_fig = figure(211);
        hold on; box on; grid on;
        vf(ii) = plot(deflection_ceil(ii).z,deflection_ceil(ii).v*1000,'Color',clrstring(ii),'linewidth',2);
                 plot(-deflection_ceil(ii).z,deflection_ceil(ii).v*1000,'Color',clrstring(ii),'linewidth',2);
        end
        
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
disp(strjoin(['Max sigma_zz at Ceiling : ' num2str(sigma_zz_MAX_ceil_val) ...
    'MPa, occurs at : ' n_allow_ceil.name(sigma_zz_MAX_ceil_ind)]))
disp(strjoin(['Min sigma_zz at Ceiling : ' num2str(sigma_zz_MIN_ceil_val) ...
    'MPa, occurs at : ' n_allow_ceil.name(sigma_zz_MIN_ceil_ind)]))
tau_sz_MAX_ceil_val = max([tau_sz_ceil(1:end).max])/1e6;
tau_sz_MAX_ceil_ind = find([tau_sz_ceil(1:end).max]/1e6 == tau_sz_MAX_ceil_val);
disp(strjoin(['Max tau_sz at Ceiling : ' num2str(tau_sz_MAX_ceil_val) ...
    'MPa, occurs at : ' n_allow_slvl.name(tau_sz_MAX_ceil_ind)]))

if PLOT_LOADS
figure(200);    xlabel('Span (m)');     ylabel('Elliptical Lift Distribution (N/m)');
                legend(lef(:),n_allow_ceil(:).name);    xlim([-b/2 b/2]);
                title('Ceiling Altitude'); pos = get(lift_ellip_fig,'position');
                set(lift_ellip_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(lift_ellip_fig,[pwd '/Load_Distribution_Figures/CEIL_Lift_Elliptical'],'-djpeg','-r300');
                
figure(201);    xlabel('Span (m)');     ylabel('Rectangular Lift Distribution (N/m)');
                legend(lrf(:),n_allow_ceil(:).name);    xlim([-b/2 b/2]);
                title('Ceiling Altitude'); pos = get(lift_rect_fig,'position');
                set(lift_rect_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(lift_rect_fig,[pwd '/Load_Distribution_Figures/CEIL_Lift_Rectangular'],'-djpeg','-r300');
                
figure(202);    xlabel('Span (m)');     ylabel('Lift Distribution (N/m)');
                legend(lf(:),n_allow_ceil(:).name);    xlim([-b/2 b/2]);
                title('Ceiling Altitude'); pos = get(lift_fig,'position');
                set(lift_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(lift_fig,[pwd '/Load_Distribution_Figures/CEIL_Lift_TTL'],'-djpeg','-r300');
                
figure(203);    xlabel('Span (m)');     ylabel('Drag Distribution (N/m)');
                legend(df(:),n_allow_ceil(:).name);    xlim([-b/2 b/2]);
                title('Ceiling Altitude'); pos = get(drag_fig,'position');
                set(drag_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(drag_fig,[pwd '/Load_Distribution_Figures/CEIL_Drag'],'-djpeg','-r300');
                
figure(204);    xlabel('Span (m)');     ylabel('w_x Distribution (N/m)');
                legend(wxf(:),n_allow_ceil(:).name);    xlim([-b/2 b/2]);
                title('Ceiling Altitude'); pos = get(wx_fig,'position');
                set(wx_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(wx_fig,[pwd '/Load_Distribution_Figures/CEIL_wx'],'-djpeg','-r300');
                
figure(205);    xlabel('Span (m)');     ylabel('w_y Distribution (N/m)');
                legend(wyf(:),n_allow_ceil(:).name);    xlim([-b/2 b/2]);
                title('Ceiling Altitude'); pos = get(wy_fig,'position');
                set(wy_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(wy_fig,[pwd '/Load_Distribution_Figures/CEIL_wy'],'-djpeg','-r300');
end

if PLOT_SHEAR
figure(206);    xlabel('Span (m)');     ylabel('S_x Distribution (kN)');
                legend(sxf(:),n_allow_ceil(:).name);    xlim([-b/2 b/2]);
                title('Ceiling Altitude'); pos = get(sx_fig,'position');
                set(sx_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(sx_fig,[pwd '/Load_Distribution_Figures/CEIL_Sx'],'-djpeg','-r300');
                
figure(207);    xlabel('Span (m)');     ylabel('S_y Distribution (kN)');
                legend(syf(:),n_allow_ceil(:).name);    xlim([-b/2 b/2]);
                title('Ceiling Altitude'); pos = get(sy_fig,'position');
                set(sy_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(sy_fig,[pwd '/Load_Distribution_Figures/CEIL_Sy'],'-djpeg','-r300');
end

if PLOT_MOMENT
figure(208);    xlabel('Span (m)');     ylabel('M_x Distribution (kNm)');
                legend(mxf(:),n_allow_ceil(:).name);    xlim([-b/2 b/2]);
                title('Ceiling Altitude'); pos = get(mx_fig,'position');
                set(mx_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(mx_fig,[pwd '/Load_Distribution_Figures/CEIL_Mx'],'-djpeg','-r300');
                
figure(209);    xlabel('Span (m)');     ylabel('M_y Distribution (kNm)');
                legend(myf(:),n_allow_ceil(:).name);    xlim([-b/2 b/2]);
                title('Ceiling Altitude'); pos = get(my_fig,'position');
                set(my_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(my_fig,[pwd '/Load_Distribution_Figures/CEIL_My'],'-djpeg','-r300');
end

if PLOT_DEFLECTION
figure(210);    xlabel('Span (m)');     ylabel('u deflection (mm)');
                legend(uf(:),n_allow_ceil(:).name);    xlim([-b/2 b/2]);
                title('Ceiling Altitude'); pos = get(u_fig,'position');
                set(u_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(u_fig,[pwd '/Deflection_Figures/CEIL_u'],'-djpeg','-r300');
                
figure(211);    xlabel('Span (m)');     ylabel('v deflection (mm)');
                legend(vf(:),n_allow_ceil(:).name);    xlim([-b/2 b/2]);
                title('Ceiling Altitude'); pos = get(v_fig,'position');
                set(v_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(v_fig,[pwd '/Deflection_Figures/CEIL_v'],'-djpeg','-r300');
end

if PLOT_SHEAR_FLOW
figure(213);    sff(end+1) = plot([0 length(tau_sz_ceil(1).qb)],[0 0],'r','linewidth',2);
                xlabel('Node (counterclockwise)');  ylabel('Shear Flow (kN/m)');
                legend(sff(:),[n_allow_slvl(:).name 'Zero']);
                title('Ceiling Altitude'); pos = get(sf_fig,'position');
                set(sf_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(sf_fig,[pwd '/Shear_Flow_Figure/CEIL_sf_open'],'-djpeg','-r300');
figure(214);    xlabel('Node (counterclockwise)');  ylabel('Shear Force (\tau_{sz}) (MPa)');
                legend(tauf(:),n_allow_slvl(:).name);
                title('Ceiling Altitude'); pos = get(tau_fig,'position');
                set(tau_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(tau_fig,[pwd '/Shear_Flow_Figure/CEIL_tau_sz'],'-djpeg','-r300');
end

if PLOT_SIGMAZZ
sigmazz_fig = figure(212);
hold on; box on; grid on;
plot(airf_geo.x,sigma_zz_slvl(sigma_zz_MAX_ceil_ind).upper/1e6,'b','linewidth',2);
plot(airf_geo.x,sigma_zz_slvl(sigma_zz_MAX_ceil_ind).lower/1e6,'g','linewidth',2);
xlabel('Chord (m)');    ylabel('Direct Stress (\sigma_{zz}) (MPa)');
legend('Upper Surface','Lower Surface');    xlim([0-Cx c-Cx]);
title('Ceiling Altitude'); pos = get(sigmazz_fig,'position');
set(sigmazz_fig,'position',[pos(1:2) pos(3:4)*1.5]);
print(sigmazz_fig,[pwd '/Stress_Figure/CEIL_DirectStress'],'-djpeg','-r300');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                        Buckling & Fatigue                         %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

buckling = calc_buckling(I_str,max([sigma_zz_MAX_ceil_val sigma_zz_MAX_slvl_val]),...
                        min([sigma_zz_MIN_ceil_val sigma_zz_MIN_slvl_val]),airf_geo.A_str,airf_geo.t_skin, airf_geo);

disp(strjoin(['Buckling Critical Stress: ' num2str(buckling.sigma_crit) 'MPa']));

if (1.5*max(abs([sigma_zz_MIN_ceil_val sigma_zz_MIN_slvl_val]))) >= buckling.sigma_crit
   disp('WING WILL BUCKLE');
else
    disp('WING WILL NOT BUCKLE!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                        Von Mises Failure                          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigma_eq = von_mises([sigma_zz_slvl(:).max],[sigma_zz_ceil(:).max],...
                    [tau_sz_slvl(:).max],[tau_sz_ceil(:).max]);

disp(strjoin(['Von Mises Equivalent Stress : ' num2str(sigma_eq.val/1e6) ...
    'MPa, occurs at : ' n_allow_slvl.name(sigma_eq.ind) ', Flight Condition: ' sigma_eq.fgt_cond]));
                
if (sigma_eq.val/1e6)*1.5 >= sigma_yield
   disp('WING WILL FAIL UNDER VON MISES STRESS CRITERIA');
else
    disp('WING WILL NOT FAIL UNDER VON MISES STRESS CRITERIA!');
end

%%% ADD IF ELSE STATEMENT HERE FOR YIELD AND BUCKLING FIRST (BEFORE
%%% CALCULATING WEIGHT)
weight_wing = calc_weight_wing(airf_geo,b,rho_material);

disp(strjoin(['Based on the current configuration, the wing weighs : ' ...
             num2str(weight_wing.total) 'kg']));
disp(strjoin(['The wing weight ' num2str(100*weight_wing.total/(mass_emp)) ...
              '% of the entire aircraft empty weight.']));

