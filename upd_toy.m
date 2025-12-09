%==========================================================================
% TOY SYNTHETIC DATA GENERATION FOR EMPIRICAL MODEL INVERSION
%
% PURPOSE
%   This helper script generates toy belief-updating data for S subjects
%   over T trials, using the static learning model with user-defined
%   “true” parameters. The outputs are stored in cell arrays formatted
%   exactly as required by the empirical analysis script (upd_data_bmc.m).
%
% DESCRIPTION
%   - Simulates trial-by-trial base rates (BR), first estimates (E1),
%     second estimates (E2), and stimulus valence (V) using upd_sim.
%   - Uses log-space parameters to match the model inversion routine.
%   - Facilitates quick testing of BPA and posterior inference pipelines
%     without needing real experimental data.
%
% OUTPUT VARIABLES
%   E1{s} : [T×1] first estimates for subject s
%   E2{s} : [T×1] second estimates for subject s
%   BR{s} : [T×1] base rates for subject s
%   V{s}  : [T×1] stimulus valence (±1) for subject s
%
% EXAMPLE USAGE
%   >> run this script
%   >> upd_data   % run the inversion and group inference
%
% AUTHOR / DATE
%   Matthew D. Greaves, University of Melbourne. 
%   Last updated: 25/11/2025.
%==========================================================================

% Toy synthetic data for S=5 subjects, T=40 trials
S = 5; T = 40; esd = 1;                 
true_a = log(1/2); 
true_b = log(2/10);

E1 = cell(S,1);
E2 = cell(S,1);
BR = cell(S,1);
V  = cell(S,1);

for s=1:S
    rng(100+s);
    % build inputs via upd_sim using true params in log space
    [y, U] = upd_sim(T, [true_a; true_b], 1/esd^2);
    
    % Unpack to cell arrays expected by the script
    BR{s} = U(1,:)';        % base rate
    E1{s} = U(2,:)';        % first estimate
    V{s}  = U(3,:)';        % valence (±1)
    E2{s} = y(:);           % second estimate (column)
end

% ------------------------------------------------------------------------- 
% Now run upd_empirical.m
S = 5; T = 40; esd = 1;                 
true_a = log(1/2); true_b = log(2/10);

E1 = cell(S,1);
E2 = cell(S,1);
BR = cell(S,1);
V  = cell(S,1);

for s=1:S
    rng(100+s);
    % build inputs via upd_sim using true params in log space
    [y, U] = upd_sim(T, [true_a; true_b], 1/esd^2);

    % Unpack to cell arrays expected by the script
    BR{s} = U(1,:)';        % base rate
    E1{s} = U(2,:)';        % first estimate
    V{s}  = U(3,:)';        % valence (±1)
    E2{s} = y(:);           % second estimate (column)
end

% ------------------------------------------------------------------------- 
% Now run upd_data.m

