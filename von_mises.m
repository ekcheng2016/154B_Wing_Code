% von_mises.m
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ sigma_eq ] = von_mises( sigma_zz_slvl,sigma_zz_ceil,tau_sz_slvl,tau_sz_ceil )

%  AT SEA LEVEL
sigma_eq_temp = sqrt((sigma_zz_slvl.^2) + (3*tau_sz_slvl.^2));
sigma_eq_slvl_val = max(sigma_eq_temp);
sigma_eq_slvl_ind = find(sigma_eq_temp == sigma_eq_slvl_val);

%  AT CEILING
sigma_eq_temp = sqrt((2*sigma_zz_slvl.^2) + (6*tau_sz_ceil.^2));
sigma_eq_ceil_val = max(sigma_eq_temp);
sigma_eq_ceil_ind = find(sigma_eq_temp == sigma_eq_ceil_val);

if sigma_eq_ceil_val > sigma_eq_slvl_val
    sigma_eq.fgt_cond = 'Ceil';
    sigma_eq.val = sigma_eq_ceil_val;
    sigma_eq.ind = sigma_eq_ceil_ind;
else
    sigma_eq.fgt_cond = 'Slvl';
    sigma_eq.val = sigma_eq_slvl_val;
    sigma_eq.ind = sigma_eq_slvl_ind;
   
end

end

