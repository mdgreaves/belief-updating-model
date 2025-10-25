function p_bpa = upd_ffx_bpa(posts, field)
% Fixed-effects Bayesian Parameter Averaging (precision pooling)
%
% p_bpa = upd_ffx_bpa(posts)
% p_bpa = upd_ffx_bpa(posts, field)
%
% Inputs
% -------
% posts : 1xS or Sx1 cell array
%     Per-subject posterior structs
%     each containing 'muPhi' and 'SigmaPhi' (or the field you choose).
%
% field : char/string (optional, default = 'Phi')
%     Which parameter block to pool: 'Phi' (default), 'Theta', or 'X0'.
%     The function reads 'mu[field]' and 'Sigma[field]' from each model.
%
% Output
% ------
% p_bpa : struct with fields
%     .mu[field]    : d x 1 group mean (log space)
%     .Sigma[field] : d x d group covariance (log space)
%     .muNat[field] : d x 1 group mean on natural scale (exp), via upd_log
%     .sdNat[field] : d x 1 group SD on natural scale
%     .covNat[field]: d x d group covariance on natural scale
%
% Notes
% -----
% * This is pure FFX BPA:  
%   Σ_g^{-1} = Σ_k Σ_k^{-1},  
%   μ_g = Σ_g (Σ_k Σ_k^{-1} μ_k)
%
% AUTHOR / DATE
%   Matthew D. Greaves, University of Melbourne. 
%   Last updated: 25/10/2025.

    if nargin < 2 || isempty(field), field = 'Phi'; end

    % Pull field names, e.g., 'muPhi', 'SigmaPhi'
    muName    = ['mu'    field];
    SigmaName = ['Sigma' field];

    S = numel(posts);
    assert(S >= 1, ['posts must be a non-empty cell array of subject',...
        ' posteriors.']);
    d = numel(posts{1}.(muName));
    Psum = zeros(d);        % Sum of precisions
    Wsum = zeros(d,1);      % Sum of precision-weighted means

    % Pool across subjects: Σ^{-1} and Σ^{-1}μ
    for s = 1:S
        mu_s = posts{s}.(muName);
        S_s  = posts{s}.(SigmaName);

        % Accumulate precision and precision-weighted mean
        Psum = Psum + (S_s \ eye(d));   % Σ_s^{-1}
        Wsum = Wsum + (S_s \ mu_s);     % Σ_s^{-1} μ_s
    end

    % Group covariance and mean (log space)
    Sigma_g = Psum \ eye(d);
    mu_g    = Sigma_g * Wsum;

    % Pack outputs using the same naming convention as VBA
    p_bpa = struct();
    p_bpa.(muName)    = mu_g;
    p_bpa.(SigmaName) = Sigma_g;

    % Also provide natural-scale moments (exact lognormal mapping)
    [muNat, sdNat, covNat] = upd_log(mu_g, Sigma_g);
    p_bpa.(['muNat'    field]) = muNat;
    p_bpa.(['sdNat'    field]) = sdNat;
    p_bpa.(['covNat'   field]) = covNat;
end
