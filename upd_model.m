function y = upd_model(~, P, u, ~)
% =========================================================================
% Observation function for a simple belief-updating model, designed for use
% with VBA_simulate and VBA_NLStateSpaceModel.
%
% The model computes an updated estimate (y) as a function of:
%   - Base rate (BR)
%   - First estimate (E1)
%   - Stimulus valence (V)
%   - Parameters: learning rate intercept (a) and bias (b)
%
% Model:
%   N  = sign(BR - E1) * V
%   LR = exp(a) + N * exp(b)
%   y  = E1 + LR * (BR - E1)
%
% Called internally by:
%   y_t = g(x, phi, u_t, inG)
%
% Inputs:
%   ~   : (unused) state vector (no hidden states in this model)
%   P   : [2×1] parameter vector, [a; b]
%   u   : [3×1] input vector at time t:
%           u(1) = BR  (base rate)
%           u(2) = E1  (initial estimate)
%           u(3) = V   (stimulus valence, ±1)
%   ~   : (unused) inG structure
%
% Output:
%   y   : [1×1] scalar, updated estimate
%
% AUTHOR / DATE
%   Matthew D. Greaves, University of Melbourne. 
%   Last updated: 25/11/2025.
% =========================================================================

% --- Unpack inputs -------------------------------------------------------
BR  = u(1);      % Base rate
E1  = u(2);      % Initial estimate
V   = u(3);      % Stimulus valence

% --- Unpack parameters ---------------------------------------------------
a   = P(1);       % Learning rate intercept
b   = P(2);       % Learning rate bias

% --- Compute valence-weighted learning rate ------------------------------
N   = sign(BR - E1) * V;    % Feedback valence
LR  = exp(a) + N * exp(b);  % Effective learning rate

% --- Update estimate -----------------------------------------------------
y   = E1 + LR * (BR - E1);  % Updated estimate (scalar)

end
