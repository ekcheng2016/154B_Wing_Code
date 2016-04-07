% find_n_conditions.m
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AoA,Cd] = find_n_conditions(n,V,rho,S,naca2415,wgt_max)
    CL = (2*n*wgt_max)/(rho*V^2*S);
    
    AoA = interp1(naca2415.CL,naca2415.alpha(naca2415.minInd:naca2415.maxInd),CL);
    Cd  = interp1(naca2415.CL,naca2415.Cd(naca2415.minInd:naca2415.maxInd),CL);
end