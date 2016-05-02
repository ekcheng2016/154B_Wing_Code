clear all; 
close all; 
clc;

PLOT_PREVIOUS = 0; % set 1 to plot and save all previous plots

clrstring = 'bgkrc';

load_aircraft_parameters;
load_conversions;

Re_sealvl = calc_Re(rho_sealvl,c,v_maneuver,mu_sealvl);
Re_alceil = calc_Re(rho_altceil,c,v_maneuver,mu_altceil);

airfoil_to_wing;
n_allow_slvl = calc_flgt_envel(naca2415(1),rho_sealvl,'Sea Level',PLOT_PREVIOUS);
n_allow_ceil = calc_flgt_envel(naca2415(2),rho_altceil,'Ceiling Altitude (14600 feet)',PLOT_PREVIOUS);

calc_centroid_momentinertia;

% SEA LEVEL Load Distributions
% AT ALL CRITICAL CONDITIONS
% NOTE:
%  wx0 and wy0 here are defined from root to tip.
%  wy and wx below are defined from tip to root.
nz = 500;
for ii = 1:length(n_allow_slvl.n)
    if ~isnan(n_allow_slvl.n(ii))
        % DETERMINE LOAD DISTRIBUTION
        [load_slvl(ii).z, ...
         load_slvl(ii).wx0     load_slvl(ii).wy0,...
         load_slvl(ii).l_ellip load_slvl(ii).l_rect,...
         load_slvl(ii).l       load_slvl(ii).d] = ...
                          calc_wxwy(n_allow_slvl.n(ii),...
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
        
        %PLOT DISTRIBUTIONS
        if PLOT_PREVIOUS
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
        shear_moment_slvl(ii) = calc_shear_moments(b, nz, load_slvl(ii).wx,...
                                    load_slvl(ii).wy,...
                                    load_slvl(ii).wx0,...
                                    load_slvl(ii).wy0);
        
        % PLOT
        if PLOT_PREVIOUS
        sx_fig = figure(106);
        hold on; box on; grid on;
        sxf(ii) = plot(shear_moment_slvl(ii).z,shear_moment_slvl(ii).Sx0,'Color',clrstring(ii),'linewidth',2);
                  plot(-shear_moment_slvl(ii).z,shear_moment_slvl(ii).Sx0,'Color',clrstring(ii),'linewidth',2);
        
        sy_fig = figure(107);
        hold on; box on; grid on;
        syf(ii) = plot(shear_moment_slvl(ii).z,shear_moment_slvl(ii).Sy0,'Color',clrstring(ii),'linewidth',2);
                  plot(-shear_moment_slvl(ii).z,shear_moment_slvl(ii).Sy0,'Color',clrstring(ii),'linewidth',2);
                  
        mx_fig = figure(108);
        hold on; box on; grid on;
        mxf(ii) = plot(shear_moment_slvl(ii).z,shear_moment_slvl(ii).Mx0,'Color',clrstring(ii),'linewidth',2);
                  plot(-shear_moment_slvl(ii).z,shear_moment_slvl(ii).Mx0,'Color',clrstring(ii),'linewidth',2);
                  
        my_fig = figure(109);
        hold on; box on; grid on;
        myf(ii) = plot(shear_moment_slvl(ii).z,shear_moment_slvl(ii).My0,'Color',clrstring(ii),'linewidth',2);
                  plot(-shear_moment_slvl(ii).z,shear_moment_slvl(ii).My0,'Color',clrstring(ii),'linewidth',2);        
        end
                  
        % DETERMINE DEFLECTIONS
        deflection_slvl(ii) = calc_deflections(b, Ixx, Iyy, Ixy,...
                                    shear_moment_slvl(ii).Mx0,...
                                    shear_moment_slvl(ii).My0,...
                                    nz,...
                                    load_slvl(ii).wx0,...
                                    load_slvl(ii).wy0);
                                
        if PLOT_PREVIOUS
        u_fig = figure(110);
        hold on; box on; grid on;
        uf(ii) = plot(deflection_slvl(ii).z,deflection_slvl(ii).u*1000,'Color',clrstring(ii),'linewidth',2);
                 plot(-deflection_slvl(ii).z,deflection_slvl(ii).u*1000,'Color',clrstring(ii),'linewidth',2);
        
        v_fig = figure(111);
        hold on; box on; grid on;
        vf(ii) = plot(deflection_slvl(ii).z,deflection_slvl(ii).v*1000,'Color',clrstring(ii),'linewidth',2);
                 plot(-deflection_slvl(ii).z,deflection_slvl(ii).v*1000,'Color',clrstring(ii),'linewidth',2);
        end
    end
end


if PLOT_PREVIOUS
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

figure(106);    xlabel('Span (m)');     ylabel('S_x Distribution (N)');
                legend(sxf(:),n_allow_slvl(:).name);    xlim([-b/2 b/2]);
                title('Sea Level'); pos = get(sx_fig,'position');
                set(sx_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(sx_fig,[pwd '/Load_Distribution_Figures/SLVL_Sx'],'-djpeg','-r300');
                
figure(107);    xlabel('Span (m)');     ylabel('S_y Distribution (N)');
                legend(syf(:),n_allow_slvl(:).name);    xlim([-b/2 b/2]);
                title('Sea Level'); pos = get(sy_fig,'position');
                set(sy_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(sy_fig,[pwd '/Load_Distribution_Figures/SLVL_Sy'],'-djpeg','-r300');
                
figure(108);    xlabel('Span (m)');     ylabel('M_x Distribution (Nm)');
                legend(mxf(:),n_allow_slvl(:).name);    xlim([-b/2 b/2]);
                title('Sea Level'); pos = get(mx_fig,'position');
                set(mx_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(mx_fig,[pwd '/Load_Distribution_Figures/SLVL_Mx'],'-djpeg','-r300');
                
figure(109);    xlabel('Span (m)');     ylabel('M_y Distribution (Nm)');
                legend(myf(:),n_allow_slvl(:).name);    xlim([-b/2 b/2]);
                title('Sea Level'); pos = get(my_fig,'position');
                set(my_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(my_fig,[pwd '/Load_Distribution_Figures/SLVL_My'],'-djpeg','-r300');

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
                
                
% CEILING ALTITUDE Load Distributions
% AT ALL CRITICAL CONDITIONS
for ii = 1:length(n_allow_ceil.n)
    if ~isnan(n_allow_ceil.n(ii))
        % DETERMINE LOAD DISTRIBUTION
        [load_ceil(ii).z, ...
         load_ceil(ii).wx0     load_ceil(ii).wy0,...
         load_ceil(ii).l_ellip load_ceil(ii).l_rect,...
         load_ceil(ii).l       load_ceil(ii).d] = ...
                          calc_wxwy(n_allow_ceil.n(ii),...
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
        
        %PLOT DISTRIBUTIONS
        if PLOT_PREVIOUS
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
        shear_moment_ceil(ii) = calc_shear_moments(b, nz, load_ceil(ii).wx,...
                                    load_ceil(ii).wy,...
                                    load_ceil(ii).wx0,...
                                    load_ceil(ii).wy0);
        
        % PLOT
        if PLOT_PREVIOUS
        sx_fig = figure(206);
        hold on; box on; grid on;
        sxf(ii) = plot(shear_moment_ceil(ii).z,shear_moment_ceil(ii).Sx0,'Color',clrstring(ii),'linewidth',2);
                  plot(-shear_moment_ceil(ii).z,shear_moment_ceil(ii).Sx0,'Color',clrstring(ii),'linewidth',2);
        
        sy_fig = figure(207);
        hold on; box on; grid on;
        syf(ii) = plot(shear_moment_ceil(ii).z,shear_moment_ceil(ii).Sy0,'Color',clrstring(ii),'linewidth',2);
                  plot(-shear_moment_ceil(ii).z,shear_moment_ceil(ii).Sy0,'Color',clrstring(ii),'linewidth',2);
                  
        mx_fig = figure(208);
        hold on; box on; grid on;
        mxf(ii) = plot(shear_moment_ceil(ii).z,shear_moment_ceil(ii).Mx0,'Color',clrstring(ii),'linewidth',2);
                  plot(-shear_moment_ceil(ii).z,shear_moment_ceil(ii).Mx0,'Color',clrstring(ii),'linewidth',2);
                  
        my_fig = figure(209);
        hold on; box on; grid on;
        myf(ii) = plot(shear_moment_ceil(ii).z,shear_moment_ceil(ii).My0,'Color',clrstring(ii),'linewidth',2);
                  plot(-shear_moment_ceil(ii).z,shear_moment_ceil(ii).My0,'Color',clrstring(ii),'linewidth',2);        
      
        end
        
        % DETERMINE DEFLECTIONS
        deflection_ceil(ii) = calc_deflections(b, Ixx, Iyy, Ixy,...
                                    shear_moment_ceil(ii).Mx0,...
                                    shear_moment_ceil(ii).My0,...
                                    nz,...
                                    load_ceil(ii).wx0,...
                                    load_ceil(ii).wy0);
        
        if PLOT_PREVIOUS
        u_fig = figure(210);
        hold on; box on; grid on;
        uf(ii) = plot(deflection_ceil(ii).z,deflection_ceil(ii).u*1000,'Color',clrstring(ii),'linewidth',2);
                 plot(-deflection_ceil(ii).z,deflection_ceil(ii).u*1000,'Color',clrstring(ii),'linewidth',2);
        
        v_fig = figure(211);
        hold on; box on; grid on;
        vf(ii) = plot(deflection_ceil(ii).z,deflection_ceil(ii).v*1000,'Color',clrstring(ii),'linewidth',2);
                 plot(-deflection_ceil(ii).z,deflection_ceil(ii).v*1000,'Color',clrstring(ii),'linewidth',2);
        end
    end
end

if PLOT_PREVIOUS
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

figure(206);    xlabel('Span (m)');     ylabel('S_x Distribution (N)');
                legend(sxf(:),n_allow_ceil(:).name);    xlim([-b/2 b/2]);
                title('Ceiling Altitude'); pos = get(sx_fig,'position');
                set(sx_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(sx_fig,[pwd '/Load_Distribution_Figures/CEIL_Sx'],'-djpeg','-r300');
                
figure(207);    xlabel('Span (m)');     ylabel('S_y Distribution (N)');
                legend(syf(:),n_allow_ceil(:).name);    xlim([-b/2 b/2]);
                title('Ceiling Altitude'); pos = get(sy_fig,'position');
                set(sy_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(sy_fig,[pwd '/Load_Distribution_Figures/CEIL_Sy'],'-djpeg','-r300');
                
figure(208);    xlabel('Span (m)');     ylabel('M_x Distribution (Nm)');
                legend(mxf(:),n_allow_ceil(:).name);    xlim([-b/2 b/2]);
                title('Ceiling Altitude'); pos = get(mx_fig,'position');
                set(mx_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(mx_fig,[pwd '/Load_Distribution_Figures/CEIL_Mx'],'-djpeg','-r300');
                
figure(209);    xlabel('Span (m)');     ylabel('M_y Distribution (Nm)');
                legend(myf(:),n_allow_ceil(:).name);    xlim([-b/2 b/2]);
                title('Ceiling Altitude'); pos = get(my_fig,'position');
                set(my_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                print(my_fig,[pwd '/Load_Distribution_Figures/CEIL_My'],'-djpeg','-r300');

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
                
% TODO: airfoil section properties
A_cap = 5/1000;   % m^2
A_str = 3/1000;   % m^2
t_spar = 0.005;   % m
t_skin = 0.0025;  % m
% % locations of spars, spar caps and stringers (nose at the origin of the coordinate)
x_spar0 = 0.25*c;                 % front spar (2 cell beam)
x_strU0 = [0.05 0.15 0.35 0.55 0.65]*c; % upper surface
x_strL0 = [0.05 0.15 0.35 0.55 0.65]*c; % lower surface
% % new coordinate with origin at the centroid is used for the output below
[Cx_test Cy_test Ixx_test Iyy_test Ixy_test] = airfoil_section(c,A_cap,A_str,...
                                                t_spar,t_skin,x_spar0,x_strU0,x_strL0);
