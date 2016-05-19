% plot_distributions.m
%
% Description:
%   Plots shears, moments, deflections, direct stress, shear stress
%   distributions
%
% Inputs:
%   Flight Condition
%   fid
%   b 
%   c
%   Cx
%   airf_geo
%   n_allow
%   shear
%   moment
%   deflection
%   sigmazz
%   tausz
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_distributions(fgt_cond,fid,b,c,Cx,airf_geo,n_allow,shear,moment,deflection,sigmazz,tausz)

clrstring = 'bgkrc';
if fgt_cond == 'Sea Level'
    fignum = 100;
    savekey = 'SLVL';
else
    fignum = 200;
    savekey = 'CEIL';
end

for ii = 1:length(n_allow.n)
    if ~isnan(n_allow.n(ii))
        % Open shear flow plot
        sf_fig = figure(fignum+13);
        hold on; box on; grid on;
        sff(ii) = plot(tausz(ii).qb/1e3,'color',clrstring(ii),'linewidth',2);
        % Shear Stress of entire skin
        tau_fig = figure(fignum+14);
        hold on; box on; grid on;
        tauf(ii) = plot(tausz(ii).skin/1e6,'color',clrstring(ii),'linewidth',2);
        % Shear characteristics
        fprintf(fid,['Flight Condition:' fgt_cond '\n']);
        fprintf(fid,'Loading Condition: %20s\n',char(n_allow.name(ii)));
        fprintf(fid,'q01 = %5.4f N/m\n',tausz(ii).q01);
        fprintf(fid,'q02 = %5.4f N/m\n',tausz(ii).q02);
        fprintf(fid,'tau_spars = %5.4f MPa, %5.4f MPa\n',tausz(ii).spar/1e6);
        fprintf(fid,'Max Shear Stress = %5.4f MPa\n',tausz(ii).max/1e6);
        fprintf(fid,'Last Open Shear Flow Value = %5.4f N/m\n',tausz(ii).qb(end));
        fprintf(fid,'Average Open Shear Flow Value = %5.4f N/m\n',mean(tausz(ii).qb));            
        fprintf(fid,'Last Open Shear Flow/Avg Shear Flow = %5.4f %% \n\n',(tausz(ii).qb(end)/mean(tausz(ii).qb))*100);
        
        % Shear in the x-direction S_x
        sx_fig = figure(fignum+6);
        hold on; box on; grid on;
        sxf(ii) = plot(shear(ii).z,shear(ii).Sx0/1e3,'Color',clrstring(ii),'linewidth',2);
                  plot(-shear(ii).z,shear(ii).Sx0/1e3,'Color',clrstring(ii),'linewidth',2);
        
        % Shear in the y-direction S_y
        sy_fig = figure(fignum+7);
        hold on; box on; grid on;
        syf(ii) = plot(shear(ii).z,shear(ii).Sy0/1e3,'Color',clrstring(ii),'linewidth',2);
                  plot(-shear(ii).z,shear(ii).Sy0/1e3,'Color',clrstring(ii),'linewidth',2);
        
        % Moment in the x-direction M_x
        mx_fig = figure(fignum+8);
        hold on; box on; grid on;
        mxf(ii) = plot(moment(ii).z,moment(ii).Mx0/1e3,'Color',clrstring(ii),'linewidth',2);
                  plot(-moment(ii).z,moment(ii).Mx0/1e3,'Color',clrstring(ii),'linewidth',2);
        
        % Moment in the y-direction M_y          
        my_fig = figure(fignum+9);
        hold on; box on; grid on;
        myf(ii) = plot(moment(ii).z,moment(ii).My0/1e3,'Color',clrstring(ii),'linewidth',2);
                  plot(-moment(ii).z,moment(ii).My0/1e3,'Color',clrstring(ii),'linewidth',2);        
        
        % Deflection in the x-direction u
        u_fig = figure(fignum+10);
        hold on; box on; grid on;
        uf(ii) = plot(deflection(ii).z,deflection(ii).u*1000,'Color',clrstring(ii),'linewidth',2);
                 plot(-deflection(ii).z,deflection(ii).u*1000,'Color',clrstring(ii),'linewidth',2);
        
        % Deflection in the y-direction v
        v_fig = figure(fignum+11);
        hold on; box on; grid on;
        vf(ii) = plot(deflection(ii).z,deflection(ii).v*1000,'Color',clrstring(ii),'linewidth',2);
                 plot(-deflection(ii).z,deflection(ii).v*1000,'Color',clrstring(ii),'linewidth',2);
     
    end 
end

figure(fignum+6);    xlabel('Span (m)');     ylabel('S_x Distribution (kN)');
                     legend(sxf(:),n_allow(:).name);    xlim([-b/2 b/2]);
                     title(fgt_cond); pos = get(sx_fig,'position');
                     set(sx_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                     print(sx_fig,[pwd '/Load_Distribution_Figures/' savekey '_Sx'],'-djpeg','-r300');
                
figure(fignum+7);    xlabel('Span (m)');     ylabel('S_y Distribution (kN)');
                     legend(syf(:),n_allow(:).name);    xlim([-b/2 b/2]);
                     title(fgt_cond); pos = get(sy_fig,'position');
                     set(sy_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                     print(sy_fig,[pwd '/Load_Distribution_Figures/' savekey '_Sy'],'-djpeg','-r300');

figure(fignum+8);    xlabel('Span (m)');     ylabel('M_x Distribution (kNm)');
                     legend(mxf(:),n_allow(:).name);    xlim([-b/2 b/2]);
                     title(fgt_cond); pos = get(mx_fig,'position');
                     set(mx_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                     print(mx_fig,[pwd '/Load_Distribution_Figures/' savekey '_Mx'],'-djpeg','-r300');
                
figure(fignum+9);    xlabel('Span (m)');     ylabel('M_y Distribution (kNm)');
                     legend(myf(:),n_allow(:).name);    xlim([-b/2 b/2]);
                     title(fgt_cond); pos = get(my_fig,'position');
                     set(my_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                     print(my_fig,[pwd '/Load_Distribution_Figures/' savekey '_My'],'-djpeg','-r300');

figure(fignum+10);    xlabel('Span (m)');     ylabel('u deflection (mm)');
                      legend(uf(:),n_allow(:).name);    xlim([-b/2 b/2]);
                      title(fgt_cond); pos = get(u_fig,'position');
                      set(u_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                      print(u_fig,[pwd '/Deflection_Figures/' savekey '_u'],'-djpeg','-r300');
                
figure(fignum+11);    xlabel('Span (m)');     ylabel('v deflection (mm)');
                      legend(vf(:),n_allow(:).name);    xlim([-b/2 b/2]);
                      title(fgt_cond); pos = get(v_fig,'position');
                      set(v_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                      print(v_fig,[pwd '/Deflection_Figures/' savekey '_v'],'-djpeg','-r300');                

figure(fignum+13);    sff(end+1) = plot([0 length(tausz(1).qb)],[0 0],'r','linewidth',2);
                      xlabel('Node (counterclockwise)');  ylabel('Shear Flow (kN/m)');
                      legend(sff(:),[n_allow(:).name 'Zero']);
                      title(fgt_cond); pos = get(sf_fig,'position');
                      set(sf_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                      print(sf_fig,[pwd '/Shear_Flow_Figure/' savekey '_sf_open'],'-djpeg','-r300');
figure(fignum+14);    xlabel('Node (counterclockwise)');  ylabel('Shear Force (\tau_{sz}) (MPa)');
                      legend(tauf(:),n_allow(:).name);
                      title(fgt_cond); pos = get(tau_fig,'position');
                      set(tau_fig,'position',[pos(1:2) pos(3:4)*1.5]);
                      print(tau_fig,[pwd '/Shear_Flow_Figure/' savekey '_tau_sz'],'-djpeg','-r300');

sigmazz_fig = figure(fignum+12);
hold on; box on; grid on;
plot(airf_geo.x,sigmazz.upper/1e6,'b','linewidth',2);
plot(airf_geo.x,sigmazz.lower/1e6,'g','linewidth',2);
xlabel('Chord (m)');    ylabel('Direct Stress (\sigma_{zz}) (MPa)');
legend('Upper Surface','Lower Surface');    xlim([0-Cx c-Cx]);
title(fgt_cond); pos = get(sigmazz_fig,'position');
set(sigmazz_fig,'position',[pos(1:2) pos(3:4)*1.5]);
print(sigmazz_fig,[pwd '/Stress_Figure/' savekey '_DirectStress'],'-djpeg','-r300');
end
                