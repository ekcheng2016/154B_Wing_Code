% calc_loads.m
%
% Description:
%   Determine the loads at sea level and ceiling altitude and plot if
%   plotting functionality is turned on
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clrstring = 'bgkrc';

% AT SEA LEVEL
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
    end
end

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

% AT CEILING ALTITUDE
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
    end
end

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