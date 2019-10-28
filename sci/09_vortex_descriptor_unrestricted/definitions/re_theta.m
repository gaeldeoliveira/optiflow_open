function [rt, rt_ue_over_nue ] = re_theta( t, ue_over_nue )
%RE_THETA returns the momentum thickness reynolds
%   rt                  reynolds theta
%   rt_ue_over_nue      derivative to ue_over_nue

% Reynolds Theta
rt = ue_over_nue .* t;
% Derivatives
% rt_ue_over_nue = t;                     % Original Derivative
rt_ue_over_nue = rt ./ ue_over_nue;       % Vectorization compatible format

end

