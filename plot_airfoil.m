% plot_airfoil.m
% 
% Description:
%   This function plots the airfoil plot (nose at origin and centered at
%   centroid)
%
% Inputs:
%   airf_geo, Cx, Cy
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = plot_airfoil(airf_geo,Cx,Cy)

% DEREFERENCE
x = airf_geo.x+Cx;
yU = airf_geo.yU+Cy;
yL = airf_geo.yL+Cy;
x_strU = airf_geo.x_strU+Cx;
x_strL = airf_geo.x_strL+Cx;
i_strU = airf_geo.i_strU;
i_strL = airf_geo.i_strL;
x_spar = airf_geo.x_spar+Cx;
i_spar = airf_geo.i_spar;


afc = figure();
plot(x,yU,'k',x,yL,'k','Linewidth',2);
ylim([-0.3 0.3])
hold on; grid on;
plot(x_strU,yU(i_strU),'or',x_strL,yL(i_strL),'or','markersize',5);
for i = 1:length(x_spar)
    plot([x_spar(i),x_spar(i)],[yU(i_spar(i)),yL(i_spar(i))],'b','Linewidth',3)
end
plot(x_spar,yU(i_spar),'sg',x_spar,yL(i_spar),'sg','markersize',6);
scatter(Cx,Cy,'m*')
title('NACA 2415 with Origin at Nose','fontsize',14);
ylabel('y (m)')
xlabel('x (m)')
print(afc,[pwd '/Airfoil_Section/REAL_Airfoil_Origin'],'-djpeg','-r300');

afc = figure();
plot(x-Cx,yU-Cy,'k',x-Cx,yL-Cy,'k','Linewidth',2);
hold on; grid on;
%  plot(x,yU,'xb',x,yL,'xb');
ylim([-0.3 0.3])
plot(x_strU-Cx,yU(i_strU)-Cy,'or',x_strL-Cx,yL(i_strL)-Cy,'or','markersize',5);
for i = 1:length(x_spar)
    plot([x_spar(i)-Cx,x_spar(i)-Cx],[yU(i_spar(i))-Cy,yL(i_spar(i))-Cy],'b','Linewidth',3)
end
plot(x_spar-Cx,yU(i_spar)-Cy,'sg',x_spar-Cx,yL(i_spar)-Cy,'sg','markersize',6);
scatter(Cx-Cx,Cy-Cy,'m*')
title('NACA 2415 Normalized to Centroid','fontsize',14);
ylabel('y (m)')
xlabel('x (m)')
print(afc,[pwd '/Airfoil_Section/REAL_Airfoil_Centroid'],'-djpeg','-r300');

end