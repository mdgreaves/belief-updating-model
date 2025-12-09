%==========================================================================
% SIMULATION–RECOVERY & MODEL CONFUSION ANALYSIS (Static Updating Model)
%
% This script (i) assesses parameter recoverability for a two-parameter
% static belief-updating model via Monte Carlo simulations across a grid of
% true parameters, and (ii) evaluates model discriminability using a model
% confusion (simulation → inversion → BMS) procedure. It also generates the
% supplementary figure with panels a–f (example trajectories, RMSE heatmap,
% and confusion matrix).
%
% WHAT IT DOES
%   1) Simulation–recovery:
%      - Simulate T trials per grid point (alpha, beta ∈ [0,1])
%      - Invert each synthetic dataset with variational Laplace
%      - Pool subject-level posteriors using FFX Bayesian parameter
%        averaging (BPA) in log-space (upd_ffx_bpa), transform to natural 
%        space, and compute RMSE over both parameters
%   2) Model confusion:
%      - Simulate data under three models: full (α,β), α-only, β-only
%      - Invert each dataset under all candidates
%      - Select by free energy (BMS by argmax F) and tabulate frequencies
%   3) Visualization:
%      - Example trial trajectories (a–d)
%      - RMSE heatmap over the parameter grid (e)
%      - Confusion matrix (f)
%
% INPUTS (script parameters below)
%   N, T             : subjects per grid point; trials per subject
%   a, b             : grids for true α and β (0…1)
%   esd              : observation noise SD used in recovery sims
%
% OUTPUTS (files / variables)
%   idnt.mat         : cell array `dat` holding posteriors & outputs
%   rmse             : RMSE(α,β) matrix after FFX BPA (natural scale)
%   confmat          : M×M confusion matrix (proportion selected)
%   on-screen figure : panels a–f (publication-style)
%
% DEPENDENCIES
%   - upd_sim.m          : generates synthetic data and inputs U
%   - upd_invert.m       : wraps VBA_NLStateSpaceModel for this model
%   - upd_ffx_bpa.m      : fixed-effects BPA in log-space (returns muNat)
%   - upd_fig.m          : helper for trial-level example panels
%   - VBA toolbox        : Variational Bayes for state-space models
%
% SOFTWARE
%   Parallelism: uses parpool + parfor (set pool profile as needed)
%
% REPRODUCIBILITY
%   Random seeds:
%     - Grid sims: rng(S(i)) within parfor for subject-specific seeds
%     - Confusion sims: rng(n) per Monte Carlo repetition
%
% USAGE
%   - Ensure all dependencies are on the MATLAB path
%   - (Optional) adjust grids N, T, a, b, esd
%   - Run the script; it will save `idnt.mat` and render the figure
%
% NOTES
%   - Priors for α,β are log-normal (see Supplementary Results)
%   - FFX BPA is used for recovery summary (precision-weighted pooling)
%   - In confusion analysis, “single-parameter” candidates are enforced by
%     tight priors on the omitted parameter (near-zero mean, tiny variance)
%
% AUTHOR / DATE
%   Matthew D. Greaves, University of Melbourne. 
%   Last updated: 25/11/2025.
%==========================================================================

% -------------------------------------------------------------------------
% Simulation-recovery analysis
% -------------------------------------------------------------------------

% Simulation parameters
N           = 20;                       % Simulations (subjects)
T           = 40;                       % Trials
[a, b]      = deal(linspace(0, 1, 10)); % True "learning" parameters
[A, B, S]   = ndgrid(a, b, 1:N);        % Parameter grid
esd         = 5;                        % Observation standard deviation
dat         = cell(size(A));            % Simulation data

if exist('idnt.mat', 'file')

    % Load previously saved output
    load('idnt.mat', 'dat');
else

% Start parallel pool
delete(gcp('nocreate'));
parpool('local');

parfor i    = 1:numel(A)
    % Set subject-specific seed for reproducibility
    rng(S(i));

    % Simulate data
    [y, u]  = upd_sim(T, log([A(i); B(i)]), 1/esd^2);

    % Invert model with simulated data
    [post, out, ~] = upd_invert(y, u);
    dat{i} = struct('post', post, 'out', out);
end

% Shut down the parallel pool
delete(gcp('nocreate'));

% Save output
save('idnt', 'dat');
end

% Use fixed-effects (FFX) Bayesian parameter-average (BPA) result used to 
% calculate RMSE across paramater grid
rmse = nan(size(A,1), size(A,2));
for i = 1:size(A,1)
    for j = 1:size(A,2)
        posts = arrayfun(@(k) dat{i,j,k}.post, 1:size(A,3), 'uni', 0)';
        p     = upd_ffx_bpa(posts, 'Phi');  % FFX pooling in log-space
        muNat = p.muNatPhi;                 % Transformed mean
        rmse(i,j) = sqrt(mean((muNat - [A(i,j,1); B(i,j,1)]).^2));
    end
end

% -------------------------------------------------------------------------
% Model confusion analysis
% -------------------------------------------------------------------------

% Identify the paramaters corresponding to the poorest paramater recovery
[~, indx] = max(rmse(:));
P         = [A(indx); B(indx)];     % True "learning" parameters
pE        = log([1/128; 1/128]);    % Prior expectation 
pC        = diag([1/4,1/4]);        % Prior covariance
M         = 3;                      % Number of models (M1, M2, M3)
px        = flip(1:M);              % Parameter selector
N         = 20;                     % Monte Carlo repetitions per model
F         = zeros(N, M, M);         % [rep, sim_model, cand_model]

for sm = 1:M
    for n = 1:N
        % Simulation model
        Psm = P;
        if px(sm) <= numel(P)
            Psm(px(sm)) = 0;
        end
        
        % Set simulation-specific seed
        rng(n);

        % Simulate with more realistic observation noise
        [y, u]  = upd_sim(T, log(Psm), 1);
        for cm = 1:M
            % Candidate model
            pEcm = pE;
            pCcm = pC;
            if px(cm) <= numel(pE)
                pEcm(px(cm)) = log(1e-6);
                pCcm(px(cm),px(cm)) = 1e-6;
            end
            
            % Invert model with simulated data
            [~, out, ~] = upd_invert(y, u, 'muPhi', pEcm,...
                'SigmaPhi', pCcm);

            % Store free energy
            F(n, sm, cm) = out.F;
        end
    end
end

% Create confusion matrix
confmat = zeros(M,M);
for sm = 1:M
    for n = 1:N
        [~, best] = max(squeeze(F(n,sm,:)));
        confmat(sm, best) = confmat(sm, best) + 1;
    end
end
confmat = confmat ./ N;   % Normalize to get frequencies

%--------------------------------------------------------------------------
% Visualisation
%--------------------------------------------------------------------------

% Parameters for panels a–d
Tv      = 10;                               % Number of example trials 
Pv      = [0, 0; 1, 0; 1/2, 0; 1/2, 1/4];   % Learning parameters
E2      = nan(size(Pv,1), Tv);              % Varible for second estimates
% Generate examples
for i = 1:size(Pv,1)
    rng(2);
    [E2(i,:), U] = upd_sim(Tv, log(Pv(i,:)), inf);
end
BR  = U(1,:)';                  % Base rate
E1  = U(2,:)';                  % Initial estimate
V   = U(3,:)';                  % Stimulus valenceaddpath
N   = sign(BR - E1) .* V;       % Feedback valence

% Font parameters
axis_font_size  = 18;
annot_font_size = 21;
note_text       = axis_font_size-3;

% Figure
fig = figure;
set(fig, 'Color', 'white', 'Units', 'normalized',...
    'Position', [0, 0, 1, 1]);

% Subplot a
ax_a = subplot(2, 3, 1);
upd_fig(BR, E1, E2(1,:)', V, N, true);
axesObjects = findall(gcf, 'Type', 'axes');
set(axesObjects, 'FontSize', axis_font_size);

% Subplot b
ax_b = subplot(2, 3, 2);
upd_fig(BR, E1, E2(2,:)', V, N, false);
axesObjects = findall(gcf, 'Type', 'axes');
set(axesObjects, 'FontSize', axis_font_size);

% Subplot c
ax_c = subplot(2, 3, 4);
upd_fig(BR, E1, E2(3,:)', V, N, false);
axesObjects = findall(gcf, 'Type', 'axes');
set(axesObjects, 'FontSize', axis_font_size);

% Subplot d
ax_d = subplot(2, 3, 5);
upd_fig(BR, E1, E2(4,:)', V, N, false);
axesObjects = findall(gcf, 'Type', 'axes');
set(axesObjects, 'FontSize', axis_font_size);

% Subplot e
ax_e = subplot(2, 3, 3);
c = [min(rmse(:)) max(rmse(:))];
imagesc(a, b, rmse'); 
axis image; 
set(gca,'YDir','normal');
colormap(gray); 
clim(c);
cb = colorbar; 
cb.Ticks = c; 
cb.TickLabels = compose('%.1g', c);
cb.FontSize = axis_font_size;
ycb = ylabel(cb, 'RMSE', 'FontSize', axis_font_size);
ycb.Position(1) = ycb.Position(1) - 2;  
xlabel('True \alpha', 'Interpreter', 'tex');
ylabel('True \beta', 'Interpreter', 'tex');
% title('RMSE over both parameters (BPA in log-space)');
xticks([0,1]); yticks([0,1]);
xticklabels([0,1]); yticklabels([0,1]);
axesObjects = findall(gcf, 'Type', 'axes');
set(axesObjects, 'FontSize', axis_font_size);

% Subplot f
ax_f = subplot(2, 3, 6);
clims = [0 1];
imagesc(confmat);
colormap(gray);
axis square;
clim(clims);
cb = colorbar;
cb.Ticks = [0,1]; 
cb.FontSize = axis_font_size;
ylabel(cb, 'Proportion selected', 'FontSize', axis_font_size);

xlabel('Selected model');
ylabel('Simulated model');
% title('Model confusion matrix');
xticks(1:M); yticks(1:M);
xticklabels({'(\alpha,\beta)','(\alpha)','(\beta)'});
yticklabels({'(\alpha,\beta)','(\alpha)','(\beta)'});
axesObjects = findall(gcf, 'Type', 'axes');
set(axesObjects, 'FontSize', axis_font_size);

% Tidy up the figure
bw = 0.05;
set(ax_a, 'Position', [bw*1.25+bw/6, 1/2+bw*1.5, 2.62/10, 3/10])
set(ax_a.Legend, 'Position', [0.0761 0.785 1/10 0.05])
set(ax_a.Legend, 'FontSize', note_text)
set(ax_b, 'Position', [1/3+bw*1.25+bw/6, 1/2+bw*1.5, 2.62/10, 3/10])
set(ax_c, 'Position', [bw*1.25+bw/6, 0+bw*3.25, 2.62/10, 3/10])
set(ax_d, 'Position', [1/3+bw*1.25+bw/6, 0+bw*3.25, 2.62/10, 3/10])
set(ax_e, 'Position', [2/3+bw+bw/6, 1/2+bw*1.5, 2/10, 3/10])
set(ax_f, 'Position', [2/3+bw+bw/6, 0+bw*3.25, 2/10, 3/10])

% Draw a dashed rectangle around the selected subplots
% Get positions of the two axes
posE = get(ax_e, 'Position');
posF = get(ax_f, 'Position');

% Compute bounding box
x0 = min(posE(1), posF(1)) - bw*0.75;
y0 = min(posE(2), posF(2)) - bw*1.5;
x1 = max(posE(1)+posE(3), posF(1)+posF(3)) + bw*0.825;
y1 = max(posE(2)+posE(4), posF(2)+posF(4)) + bw*1.2;

% Rectangle position [x y width height]
rectPos = [x0, y0, x1 - x0, y1 - y0];

% Draw rectangle
annotation('rectangle', rectPos, ...
    'LineStyle', '--', ...
    'LineWidth', 1.5, ...
    'Color', [0.5, 0.5, 0.5]);

% Small annotations
text(ax_a, 0.7, 1.08, ['\it{\alpha}\rm = ', num2str(0), ', ',...
    '\it{\beta}\rm = ', num2str(0)],...
    'Units', 'normalized', 'FontSize',...
    axis_font_size, 'FontWeight', 'normal');
text(ax_b, 0.7, 1.08, ['\it{\alpha}\rm = ', num2str(1), ', ',...
    '\it{\beta}\rm = ', num2str(0)],...
    'Units', 'normalized', 'FontSize',...
    axis_font_size, 'FontWeight', 'normal');
text(ax_c, 0.661, 1.08, ['\it{\alpha}\rm = ', num2str(1/2), ', ',...
    '\it{\beta}\rm = ', num2str(0)],...
    'Units', 'normalized', 'FontSize',...
    axis_font_size, 'FontWeight', 'normal');
text(ax_d, 0.6, 1.08, ['\it{\alpha}\rm = ', num2str(1/2), ', ',...
    '\it{\beta}\rm = ', num2str(1/4)],...
    'Units', 'normalized', 'FontSize',...
    axis_font_size, 'FontWeight', 'normal');

% Annotation
text(ax_a, -0.1, 1.125, 'a', 'Units', 'normalized', 'FontSize',...
    annot_font_size, 'FontWeight', 'bold');
text(ax_b, -0.1, 1.125, 'b', 'Units', 'normalized', 'FontSize',...
    annot_font_size, 'FontWeight', 'bold');
text(ax_c, -0.1, 1.125, 'c', 'Units', 'normalized', 'FontSize',...
    annot_font_size, 'FontWeight', 'bold');
text(ax_d, -0.1, 1.125, 'd', 'Units', 'normalized', 'FontSize',...
    annot_font_size, 'FontWeight', 'bold');
text(ax_e, -0.1, 1.125, 'e', 'Units', 'normalized', 'FontSize',...
    annot_font_size, 'FontWeight', 'bold');
text(ax_f, -0.1, 1.125, 'f', 'Units', 'normalized', 'FontSize',...
    annot_font_size, 'FontWeight', 'bold');
