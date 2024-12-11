function [uxW, uxT, uyW, uyT] = RSS(alphaW, alphaT)
% RSS - Computes confidence intervals on wing and tail Cl by propagating 
% uncertainties on measures of alpha and iT angles (see paper page 3).
% Uniform distributions are assumed for measures errors.
% 
% INPUTS: 
% - alphaW, float: uncertainty on wing alpha measure, default +-0.2° (deg)
% - alphaT, float: uncertainty on tail alpha measure, default +-0.4° (deg)
%
% OUTPUTS:
% - uyW, float: propagated uncertainty on wing cl
% - uyT, float: propagated uncertainty on tail cl
%
% CALLED FUNCTIONS: -
%
% REVISIONS:
% - #v0 11/12/24, Boscariol Jacopo
%               Changes: release.
    
    % using default paper values if not specified
    if ~nargin
        alphaW = 0.2;
        alphaT = 0.4;
    end

    % page 3 of the paper, uniform dist assumed
    uxW = deg2rad(alphaW)/sqrt(3);

    % page 3 considering also iT uncertainty, uniform dist assumed
    uxT = deg2rad(alphaT)/sqrt(3);
    
    % RSS derivatives (cl = 2*pi*alpha + ...)
    dydx = 2*pi;
    
    uyW = dydx*uxW;
    uyT = dydx*uxT;

end
