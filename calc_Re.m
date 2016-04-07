% calc_Re.m
%
% Calculates Reynolds Number
%   
% Inputs:
%   alt : altitude (m)
%   l   : effective length (m)
%   v   : velocity (m/s)
%   mu  : dynamic viscosity (kg/ms)
%
% Outputs:
%   Re  : Reynolds Number
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Re ] = calc_Re(rho,l,v,mu)

    Re = (rho*l*v)/mu;
    
end