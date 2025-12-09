%==========================================================================
% EMPIRICAL ANALYSIS: STATIC LEARNING MODEL + MODEL COMPARISON
%  - Fixed-effects Bayesian Parameter Averaging (BPA) on full model
%  - One-tailed posterior inference for positivity-constrained parameters
%  - Random-effects Bayesian Model Comparison (RFX BMC) over 3 models
%
% PURPOSE
%   Given empirical belief-updating data, this script:
%     (i)   inverts three parameterisations of the static learning model
%           for each participant:
%              M1: full model          (α, β free)
%              M2: learning-only       (α free, β ≈ 0)
%              M3: valence-only        (β free, α ≈ 0)
%     (ii)  pools the posteriors from the full model (M1) using FFX BPA
%           and reports group-level learning parameters; and
%     (iii) compares the three models at the group level using
%           random-effects Bayesian model comparison (VBA_groupBMC),
%           printing expected model frequencies (Ef) and exceedance
%           probabilities (ep).
%
% DATA LAYOUT
%   We assume S subjects with subject-varying trial counts, missing trials,
%   or pre-cleaned exclusions. Cell arrays let each subject keep their own
%   [T_s × 1] vectors without forcing padding/NaNs.
%     E1{s} : first estimates
%     E2{s} : second estimates
%     BR{s} : base rates
%     V{s}  : stimulus valence (±1)
%
% MODEL INVERSION
%   For each subject s, build U = [BR; E1; V] (3×T_s) and y = E2 (1×T_s),
%   then call upd_invert(y, U, ...) under 3 different priors:
%     M1: default lognormal priors on [log α; log β]
%     M2: same, but β tightly shrunk to ~0
%     M3: same, but α tightly shrunk to ~0
%
% GROUP POOLING (FFX BPA on full model)
%   BPA forms the precision-weighted mean in log space for the full model:
%     μ_BPA = (∑_s Σ_s^{-1})^{-1} (∑_s Σ_s^{-1} μ_s)
%   where μ_s, Σ_s are the subject post mean/covariance (log-scale).
%   We then report both log-scale and natural-scale (exp-transformed) 
%   moments, plus one-tailed probabilities that α, β > 0.
%
% RFX BMC (VBA_groupBMC)
%   We pass the subject-wise free energies F(s,m) into VBA_groupBMC(F'),
%   which returns:
%     Ef : expected model frequencies
%     ep : model exceedance probabilities
%   summarising which parameterisation best explains the data at the
%   population level.
%
% DEPENDENCIES
%   upd_invert.m      — model inversion wrapper
%   upd_ffx_bpa.m     — FFX BPA helper (returns μ/Σ in log/natural space)
%   VBA_groupBMC.m    — random-effects BMC (from the VBA toolbox)
%
% AUTHOR / DATE
%   Matthew D. Greaves, University of Melbourne. 
%   Last updated: 26/11/2025.
%==========================================================================

% --------------------------- USER INPUTS ---------------------------------
% Provide these cell arrays in your workspace (e.g., loaded from .mat):
%   E1, E2, BR, V  (each S×1 cell; contents are column vectors)

assert(iscell(E1) && iscell(E2) && iscell(BR) && iscell(V), ...
    'E1, E2, BR, V must be S×1 cell arrays.');

S = numel(E1);

% ---------------------- PRIORS FOR THREE MODELS --------------------------
% Parameters are Phi = [log α; log β]
% Full-model priors (log-normal on α, β)
pE_full = log([1/128; 1/128]);   % prior mean on log(α), log(β)
pC_full = diag([1/4, 1/4]);      % prior covariance

M = 3;                           % number of models

muPhi    = cell(M,1);
SigmaPhi = cell(M,1);

% M1: full (α, β free)
muPhi{1}    = pE_full;
SigmaPhi{1} = pC_full;

% M2: α-only (β ≈ 0, tightly shrunk)
muPhi{2}    = pE_full;
SigmaPhi{2} = pC_full;
muPhi{2}(2)      = log(1e-6);    % β prior mean near 0
SigmaPhi{2}(2,2) = 1e-6;         % very small variance on β

% M3: β-only (α ≈ 0, tightly shrunk)
muPhi{3}    = pE_full;
SigmaPhi{3} = pC_full;
muPhi{3}(1)      = log(1e-6);    % α prior mean near 0
SigmaPhi{3}(1,1) = 1e-6;         % very small variance on α

% ------------------ INVERSION PER SUBJECT × MODEL ------------------------
posts = cell(S, M);          % posterior structs per subject × model
F     = nan(S, M);           % free energy per subject × model

for s = 1:S
    % Basic shape checks / coercions
    y = E2{s}(:)';                            % 1×T
    U = [BR{s}(:)'; E1{s}(:)'; V{s}(:)'];     % 3×T (rows: BR; E1; V)

    for m = 1:M
        [post, out, ~] = upd_invert(y, U, ...
            'muPhi',    muPhi{m}, ...
            'SigmaPhi', SigmaPhi{m});

        posts{s,m} = post;
        F(s,m)     = out.F;
    end
end

% ------------------ FIXED-EFFECTS BPA (full model only) ------------------
% Use the subject posteriors from the full model (M1) for BPA
posts_full = posts(:,1);              % S×1 cell, each entry = post for M1
p = upd_ffx_bpa(posts_full, 'Phi');   % BPA in log space

% ------------- Use group moments on the natural scale --------------------
mu_nat = p.muNatPhi;          % [alpha; beta] means in natural space
sd_nat = p.sdNatPhi;          % [alpha; beta] SDs in natural space

% Posterior probability param > 0 (Gaussian approx on natural scale)
pp_alpha_gt0 = 1 - normcdf(0, mu_nat(1), sd_nat(1));
pp_beta_gt0  = 1 - normcdf(0, mu_nat(2), sd_nat(2));

% (Optional) also report 95% Gaussian CIs on natural scale
z = 1.96;
ci_alpha_nat = mu_nat(1) + z*[-1 1]*sd_nat(1);
ci_beta_nat  = mu_nat(2) + z*[-1 1]*sd_nat(2);

% ---------------------- FFX SUMMARY: PRINT TO TERMINAL -------------------
fprintf('=============================================================\n');
fprintf('Static learning model: group-level parameters (full model)\n');
fprintf('-------------------------------------------------------------\n');
fprintf(['Natural scale (group): mean = [α=%.3f, β=%.3f], ',...
         'SD = [%.3f, %.3f]\n'], ...
        mu_nat(1), mu_nat(2), sd_nat(1), sd_nat(2));
fprintf('P(α > 0) = %.3f,  P(β > 0) = %.3f\n', ...
        pp_alpha_gt0, pp_beta_gt0);
fprintf('95%% CI α = [%.3f, %.3f], β = [%.3f, %.3f]\n', ...
        ci_alpha_nat(1), ci_alpha_nat(2), ci_beta_nat(1), ci_beta_nat(2));
fprintf('=============================================================\n');

% -------------------- FIXED-EFFECTS MODEL COMPARISON ---------------------
[~, winModel] = max(F, [], 2);          % S×1: winning model per subject
win_counts    = accumarray(winModel, 1, [M 1]);
F_group       = sum(F, 1);              % 1×M: sum log-evidence per model

F_shift       = F_group - max(F_group); % numerical stability
model_post_FFX = exp(F_shift) ./ sum(exp(F_shift));

fprintf('Fixed-effects model comparison (by summed free energy)\n');
fprintf('-------------------------------------------------------------\n');
fprintf('Number of subjects: %d\n', S);
fprintf('Models:\n');
fprintf('  1: Full (α, β free)\n');
fprintf('  2: α-only (β shrunk to ~0)\n');
fprintf('  3: β-only (α shrunk to ~0)\n\n');

for m = 1:M
    fprintf(['Model %d (FFX): n_wins = %d, sum F = %.2f,',...
        ' p(model|Y)_FFX ≈ %.3f\n'], ...
        m, win_counts(m), F_group(m), model_post_FFX(m));
end
fprintf('-------------------------------------------------------------\n');

% ---------------- RFX BAYESIAN MODEL COMPARISON (VBA_groupBMC) ----------
optionsBMC = struct;
optionsBMC.DisplayWin = 0;       % no GUI by default
optionsBMC.verbose    = 1;
optionsBMC.families   = [];      % no model families here
optionsBMC.figName    = 'Static learning: group BMC';
optionsBMC.modelNames = { ...
    '(\alpha,\beta)', ...   % Model 1: full
    '(\alpha)',       ...   % Model 2: α-only
    '(\beta)'         ...   % Model 3: β-only
    };

[postRFX, outRFX] = VBA_groupBMC(F', optionsBMC);   % F': models × subjects

Ef = outRFX.Ef;        % expected model frequencies (K×1)
ep = outRFX.ep;        % exceedance probabilities (K×1)

fprintf('Random-effects model comparison (RFX, VBA_groupBMC)\n');
fprintf('-------------------------------------------------------------\n');
for m = 1:M
    fprintf(['Model %d (RFX): Ef = %.3f, ep = %.3f ', ...
             '(n_wins = %d, sum F = %.2f)\n'], ...
        m, Ef(m), ep(m), win_counts(m), F_group(m));
end
fprintf('-------------------------------------------------------------\n');
fprintf(['Note: Ef = expected model frequency;',...
    ' ep = exceedance probability.\n']);
fprintf('=============================================================\n');
