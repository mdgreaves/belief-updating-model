function [y, U] = upd_sim(T, P, sigma, U)
% upd_sim — simulate y for the static updating model via VBA_simulate.
% Inputs:
%   T      : number of time points
%   P      : [a; b] (log-scale params)
%   sigma  : observation precision
%   U      : (optional) inputs, 3×T or T×3 with rows [BR; E1; V]
% Outputs:
%   y      : 1×T simulated observations
%   U      : 3×T inputs used
%
% AUTHOR / DATE
%   Matthew D. Greaves, University of Melbourne. 
%   Last updated: 25/11/2025.

if nargin < 4 || isempty(U)
    BR = 100*rand(1,T);
    E1 = 100*rand(1,T);
    V  = sign(randn(1,T));      % ±1
    U  = [BR; E1; V];           % 3×T
else
    if size(U,1) ~= 3 && size(U,2) == 3, U = U.'; end
    if size(U,1) ~= 3 || size(U,2) ~= T
        error('U must be 3×T (rows: BR; E1; V).');
    end
end

y = VBA_simulate(T, [], @upd_model, [], P, U, Inf, sigma, struct,...
    []);

if iscolumn(y), y = y.'; end

end
