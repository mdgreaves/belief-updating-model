%==========================================================================
% EMPIRICAL ANALYSIS: STATIC LEARNING MODEL
%  - Fixed-effects Bayesian Parameter Averaging (BPA)
%  - One-tailed posterior inference for positivity-constrained parameters
%
% PURPOSE
%   Invert the static belief-updating model per participant, pool parameter
%   posteriors with precision-weighted FFX BPA, and test whether the
%   group-level learning parameters are credibly > 0.
%
% DATA LAYOUT (Why cell arrays?)
%   We assume S subjects with subject-varying trial counts, missing trials,
%   or pre-cleaned exclusions. Cell arrays let each subject keep their own
%   [T_s × 1] vectors without forcing padding/NaNs.
%     E1{s} : first estimates
%     E2{s} : second estimates
%     BR{s} : base rates
%     V{s}  : stimulus valence (±1)
%
% MODEL INVERSION
%   For each subject s, build U = [BR; E1; V] (3×T_s) and y = E2 (1×T_s)
%   and call upd_invert(y,U) to obtain subject posterior 'post{s}'.
%
% GROUP POOLING (FFX BPA)
%   BPA forms the precision-weighted mean in log space:
%     μ_BPA = (∑_s Σ_s^{-1})^{-1} (∑_s Σ_s^{-1} μ_s)
%   where μ_s, Σ_s are the subject post mean/covariance (log-scale).
%   We then report both log-scale and natural-scale (exp-transformed) 
%   moments.
%
% INFERENCE (one-tailed, Wald-type)
%   Parameters have lognormal priors and are exponentiated in g(·),
%   implying support on (0,∞). Testing “> 0” is the coherent hypothesis:
%     pp = 1 - normcdf(0, μ_BPA(i), sqrt(Σ_BPA(i,i)))
%   We flag “credible > 0” if pp ≥ 0.95 (heuristic threshold).
%
% DEPENDENCIES
%   upd_invert.m    — model inversion wrapper
%   upd_ffx_bpa.m   — FFX BPA helper (returns μ/Σ)
%
% AUTHOR / DATE
%   Matthew D. Greaves, University of Melbourne. 
%   Last updated: 25/10/2025.
%==========================================================================

% --------------------------- USER INPUTS ---------------------------------
% Provide these cell arrays in your workspace (e.g., loaded from .mat):
%   E1, E2, BR, V  (each S×1 cell; contents are column vectors)

assert(iscell(E1) && iscell(E2) && iscell(BR) && iscell(V), ...
    'E1, E2, BR, V must be S×1 cell arrays.');

S = numel(E1);

% ------------------------- INVERSION PER SUBJECT -------------------------
posts = cell(S,1);
for s = 1:S
    % Basic shape checks / coercions
    y = E2{s}(:)';                % 1×T
    U = [BR{s}(:)'; E1{s}(:)'; V{s}(:)'];   % 3×T  (rows: BR; E1; V)
    % Invert
    [post, ~, ~] = upd_invert(y, U);        % uses default lognormal priors
    posts{s} = post;
end

% ------------------ FIXED-EFFECTS BPA (log space) ------------------------
p = upd_ffx_bpa(posts, 'Phi');

% ------------- Use group moments on the natural scale --------------------
mu_nat = p.muNatPhi;          % [alpha; beta] means in natural space
sd_nat = p.sdNatPhi;          % [alpha; beta] SDs in natural space

% Posterior probability param > 0 (Gaussian approx on natural space)
pp_alpha_gt0 = 1 - normcdf(0, mu_nat(1), sd_nat(1));
pp_beta_gt0  = 1 - normcdf(0,  mu_nat(2), sd_nat(2));

% (Optional) also report 95% Gaussian CIs on natural scale
z = 1.96;
ci_alpha_nat = mu_nat(1) + z*[-1 1]*sd_nat(1);
ci_beta_nat  = mu_nat(2) + z*[-1 1]*sd_nat(2);

% ------------------------------ OUTPUT -----------------------------------
% Print
fprintf(['Natural scale (group): mean = [%.3f, %.3f], ',...
    'SD = [%.3f, %.3f]\n'], ...
        mu_nat(1), mu_nat(2), sd_nat(1), sd_nat(2));
fprintf('P(alpha > 0) = %.3f,  P(bias > 0) = %.3f\n', ...
        pp_alpha_gt0, pp_beta_gt0);
fprintf('95%% CI alpha = [%.3f, %.3f], beta = [%.3f, %.3f]\n', ...
        ci_alpha_nat(1), ci_alpha_nat(2), ci_beta_nat(1), ci_beta_nat(2));




