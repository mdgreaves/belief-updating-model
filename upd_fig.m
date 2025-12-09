function upd_fig(BR, E1, E2, V, N, lgd_flag)
%==========================================================================
% upd_fig  —  Trial-wise visualization for the static updating model
%
% PURPOSE
%   Plot first estimates (E1), base rates (BR), and second estimates (E2)
%   across trials, with overlaid markers indicating stimulus valence (V)
%   and feedback valence (N). Designed for panels a–d in the supplement.
%
% WHAT IT DRAWS
%   • Lines (with circle markers) for:
%       - E1  : medium gray [0.5 0.5 0.5]
%       - BR  : dark   gray [0.3 0.3 0.3]
%       - E2  : dashed, magenta-leaning [0.90 0.10 0.60]
%   • Filled circles at y=100 for stimulus valence V∈{-1,+1}:
%       - V = +1 (positive) → red   [1.0 0.0 0.0]
%       - V = −1 (negative) → blue  [0.1 0.2 0.6]
%   • Filled triangles at y=95 for feedback valence N:
%       - N = +1 (favourable) → red   [1.0 0.0 0.0]
%       - N = −1 (unfavourable) → blue [0.1 0.2 0.6]
%       - N =  NaN (undefined/neutral) → orange [1.0 0.5 0.0]
%
% SIGN CONVENTIONS (for reference)
%   Feedback valence N is typically computed as:
%       N = sign(BR - E1) .* V
%   (i.e., whether the base rate implies “better/worse than expected”
%    given the stimulus valence). This function accepts N precomputed.
%
% INPUTS
%   BR      : T×1 (or 1×T) vector of base rates (0–100)
%   E1      : T×1 (or 1×T) vector of first estimates (0–100)
%   E2      : T×1 (or 1×T) vector of second estimates (0–100)
%   V       : T×1 (or 1×T) vector of stimulus valence in {-1,+1}
%   N       : T×1 (or 1×T) vector of feedback valence in {-1,+1} (NaN ok)
%   lgd_flag: logical, if true shows legend for E1/BR/E2
%
% OUTPUTS
%   None. Plots into the current axes and sets basic axis/label formatting.
%
% AXES / STYLE
%   • y-limits fixed to [0,100], y-ticks at [0,100]
%   • x-limits span the trial range with half-padding
%   • y-label is nudged rightward for readability
%   • If lgd_flag is true, legend is placed at 'northwest'
%
% USAGE EXAMPLE
%   T  = 10;
%   BR = 100*rand(T,1);
%   E1 = 100*rand(T,1);
%   E2 = 100*rand(T,1);
%   V  = randsrc(T,1,[-1 1; 0.5 0.5]);                 % ±1
%   N  = sign(BR - E1) .* V;                           % feedback valence
%   figure; upd_fig(BR, E1, E2, V, N, true);
%
% NOTES
%   • Assumes percentage scale [0,100]. Adjust ylim/markers if not.
%   • If you prefer grayscale-only figures, replace the RGBs accordingly.
%   • This function does not call hold off; it leaves the axes in its
%     current hold state after plotting.
%
% AUTHOR / DATE
%   Matthew D. Greaves, University of Melbourne. 
%   Last updated: 25/11/2025.
%==========================================================================

% Define colours for stimulus valence
Cs               = zeros(length(V), 3);
Cs(V == -1, :)   = repmat([0.1, 0.2, 0.6], sum(V == -1), 1); % Blue for -1
Cs(V == 1, :)    = repmat([1, 0, 0], sum(V == 1), 1);        % Red for 1

% Define colours for feedback valence
Cf               = zeros(length(N), 3);
Cf(N == -1, :)   = repmat([0.1, 0.2, 0.6], sum(N == -1), 1); % Blue for -1
Cf(N == 1, :)    = repmat([1, 0, 0], sum(N == 1), 1);        % Red for 1
Cf(isnan(N), :)  = repmat([1, 0.5, 0], sum(isnan(N)), 1);    % Orange for 0

% Create figure
h1 = plot(1:length(N), [E1, BR, E2], 'LineWidth', 2.75, 'Marker', 'o',...
    'MarkerSize', 7.5, 'MarkerFaceColor','auto');
set(h1(1), 'Color', [0.5, 0.5, 0.5]);
set(h1(2), 'Color', [0.3, 0.3, 0.3]);
set(h1(end), 'LineStyle', '--', 'Color', [0.90, 0.10, 0.60]); 
hold on;
scatter(1:length(V), 100 * ones(length(V), 1), 120, Cs, 'filled', 'o');
hold on;
scatter(1:length(N), 95 * ones(length(N), 1), 120, Cf, 'filled', 'v');
ylim([0, 100]);
yticks([0, 100]);
xlim([1/2, length(N)+1/2]);
xticks(1:length(N));
hY = ylabel('Belief (%)', 'FontSize', 12);
hY.Position(1) = hY.Position(1) + 1/2;   
xlabel('Trials', 'FontSize', 12);
set(gca, 'FontSize', 12);
if lgd_flag
legend(h1, ...
    {'First estimate',...
    'Base rate',...
    'Second estimate'},...
    'Location', 'northwest', 'Interpreter', 'tex');
end


end
