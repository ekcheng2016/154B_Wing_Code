% airfoil_to_wing.m
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
condition = cellstr(['slvl';'ceil']);
file_in = './Airfoil_Data/';
airfoil_data = strcat(file_in,'naca2415_',condition,'.txt');

% Airfoil(ii) is airfoil listed in airfoils vector
% Text file is formatted in columns:
% Alpha (deg) CL CD CDp CM Top_Xtr Bot_Xtr
for ii = 1:length(airfoil_data)
    raw = dlmread(char(airfoil_data(ii)));
   
    naca2415(ii).condition = condition(ii);
    naca2415(ii).alpha     = raw(:,1)*deg2rad;
    naca2415(ii).Cl        = raw(:,2);
    naca2415(ii).Cd        = raw(:,3);
    naca2415(ii).Cdp       = raw(:,4);
    naca2415(ii).CM        = raw(:,5);

    % Calculate the zero lift angle of attack --> alpha_0 (radians)
    tmp = abs(naca2415(ii).Cl-0);                   % find when Cl = 0
    [idx idx] = min(tmp);                           % locate the index
    naca2415(ii).alpha0 = naca2415(ii).alpha(idx);  % in radians

    % calculate max/min Cl, max/min alpha
    [maxCl maxInd] = max(naca2415(ii).Cl);
    naca2415(ii).Clmax    = maxCl;
    naca2415(ii).maxAlpha = naca2415(ii).alpha(maxInd);
    [minCl minInd] = min(naca2415(ii).Cl);
    naca2415(ii).Clmin    = minCl;
    naca2415(ii).minAlpha = naca2415(ii).alpha(minInd);
    
    naca2415(ii).minInd = minInd;
    naca2415(ii).maxInd = maxInd;

    % Determine indices where AoA is -5 and 5 degrees; where all airfoil
    % have straight line data
    start_ind = find(naca2415(ii).alpha == -5*deg2rad);
    end_ind = find(naca2415(ii).alpha == 5*deg2rad);

    % Calculate the slope both /rad --> Cl_alpha
    naca2415(ii).Cl_alpha = (naca2415(ii).Cl(end_ind)-naca2415(ii).Cl(start_ind))/...
                            (naca2415(ii).alpha(end_ind)-naca2415(ii).alpha(start_ind));
end

% Plot Cl vs. Alpha graph (airfoil)
fig = figure(1);
hold on; grid on;
plot(naca2415(1).alpha,naca2415(1).Cl,'--r');
plot(naca2415(2).alpha,naca2415(2).Cl,'--k');
xlabel('\alpha (rad)','FontSize',12); ylabel('Coefficient of Lift','FontSize',12);

for ii = 1:length(airfoil_data)
    % Calculate the slope /rad --> CL_alpha
    naca2415(ii).CL_alpha = naca2415(ii).Cl_alpha/(1+(naca2415(ii).Cl_alpha/(pi*AR*e)));

    %Calculate the 3D CL values
    naca2415(ii).CL = (naca2415(ii).CL_alpha/naca2415(ii).Cl_alpha)*naca2415(ii).Cl;
    
    % Calculate the max/min CL value (stall CL)
%    naca2415(ii).CLmax = naca2415(ii).CL(end);
%    naca2415(ii).CLmin = naca2415(ii).CL(1);
    naca2415(ii).CLmax = (naca2415(ii).CL_alpha/naca2415(ii).Cl_alpha)*naca2415(ii).Clmax;
    naca2415(ii).CLmin = (naca2415(ii).CL_alpha/naca2415(ii).Cl_alpha)*naca2415(ii).Clmin;

    if ii == 1
        plot(naca2415(ii).alpha,naca2415(ii).CL,'r','LineWidth',2);
    else
        plot(naca2415(ii).alpha,naca2415(ii).CL,'k','LineWidth',2);
    end
end

legend({'Sea Level C_l','Ceiling C_l','Sea Level C_L','Ceiling C_L'},'FontSize',12,'Location','southeast');
title('C_L vs. \alpha','FontSize',14);

print(fig,[pwd '/Lift_Curve_Slope_Figure/NACA2415_LiftCurve_Figure'],'-djpeg');
