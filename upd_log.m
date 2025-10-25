function [mu_nat, sd_nat, cov_nat] = upd_log(muPhi, SigmaPhi)
% upd_log — exact moments for exp-transformed Gaussian params.
% Inputs:
%   muPhi     : n×1 mean in log-space
%   SigmaPhi  : n×n covariance in log-space
% Outputs:
%   mu_nat    : n×1 E[exp(phi)]
%   sd_nat    : n×1 SD[exp(phi)]
%   cov_nat   : n×n Cov[exp(phi)]
%
% Example:
%   [m,s,c] = upd_log(post.muPhi, post.SigmaPhi);
%
% AUTHOR / DATE
%   Matthew D. Greaves, University of Melbourne. 
%   Last updated: 25/10/2025.

muPhi = muPhi(:);
n = numel(muPhi);
v = diag(SigmaPhi);
mu_nat = exp(muPhi + 0.5*v);

cov_nat = zeros(n,n);
for i = 1:n
    for j = 1:n
        cov_nat(i,j) = exp(muPhi(i)+muPhi(j) +...
            0.5*(SigmaPhi(i,i)+SigmaPhi(j,j)+2*SigmaPhi(i,j))) ...
            - mu_nat(i)*mu_nat(j);
    end
end
sd_nat = sqrt(diag(cov_nat));
end
