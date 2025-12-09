function [posterior, out, yhat] = upd_invert(y, U, varargin)
% upd_invert — invert the static updating model with VBA_NLStateSpaceModel.
% Inputs:
%   y         : 1×T (or T×1) observations
%   U         : 3×T (or T×3) inputs [BR; E1; V]
% Name–Value (all optional):
%   'muPhi'    : 2×1 prior mean for [a; b]         (default: [0; 0])
%   'SigmaPhi' : 2×2 prior covariance              (default: diag([4,4]))
%   'a_sigma'  : prior shape for obs precision     (default: 1)
%   'b_sigma'  : prior rate  for obs precision     (default: 1)
%   'g'        : observation fn handle             (default: @upd_model)
%
% Outputs:
%   posterior  : VBA posterior struct
%   out        : VBA output struct (use out.F for free energy, out.fit.R2)
%   yhat       : 1×T model predictions
%
% AUTHOR / DATE
%   Matthew D. Greaves, University of Melbourne. 
%   Last updated: 25/11/2025.

% ---- Parse options ------------------------------------------------------
p = inputParser;
addParameter(p,'muPhi', log([1/128; 1/128]));
addParameter(p,'SigmaPhi',diag([1/4,1/4]));
addParameter(p,'a_sigma',1);
addParameter(p,'b_sigma',1);
addParameter(p,'g',@upd_model);
parse(p,varargin{:});
opt = p.Results;

% ---- Shape checks -------------------------------------------------------
if iscolumn(y), y = y.'; end
if size(U,1) ~= 3 && size(U,2) == 3, U = U.'; end
assert(size(U,1)==3, 'U must be 3×T (rows: BR; E1; V).');
assert(size(U,2)==numel(y), 'U and y must have the same T.');

% ---- Dimensions & priors ------------------------------------------------
dim = struct('n',0,'n_theta',0,'n_phi',2,'p',1,'u',3);
options = struct();
options.priors.muPhi    = opt.muPhi(:);
options.priors.SigmaPhi = opt.SigmaPhi;
options.priors.a_sigma  = opt.a_sigma;
options.priors.b_sigma  = opt.b_sigma;

% ---- Turn off plotting & printing ---------------------------------------
options.DisplayWin = 0;   % turn off the GUI
options.verbose    = 1;   % silence command-line output
options.GnFigs     = 0;   % suppress figure generation

% ---- Inversion ----------------------------------------------------------
[posterior, out] = VBA_NLStateSpaceModel(y, U, [], opt.g, dim, options);

% ---- Predictions --------------------------------------------------------
if isfield(out,'suffStat') && isfield(out.suffStat,'gx')
    yhat = out.suffStat.gx;
else
    yhat = nan(size(y));
    for t = 1:numel(y)
        yhat(t) = opt.g([], posterior.muPhi, U(:,t), []);
    end
end
end
