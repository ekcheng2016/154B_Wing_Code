% find_n_conditions.m
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AoA,Cd] = find_n_conditions(n,V,rho,S,naca2415,wgt_max)
    CL = (2*n*wgt_max)/(rho*V^2*S)
    
    tmp = abs(naca2415.CL-CL);
    [idx idx] = min(tmp);
    
%    AoA = interp1(naca2415.CL,naca2415.alpha,CL,'linear');
%    Cd  = interp1(naca2415.CL,naca2415.Cd,CL,'linear');

    AoA = naca2415.alpha(idx);
    Cd  = naca2415.Cd(idx);
end