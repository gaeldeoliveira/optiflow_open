function [ d , d_h, d_hh, d_t ] = delta_hh( h, hh, t)
%delta_hh Returns BL thickness from h, hh and t
%   h   Shape Factor
%   hh  Heads H1 Shape Factor
%   t   Boundary Layer Thickness (theta)

% Boundary Layer Thickness
d = (h + hh) .* t;

% Derivatives
d_h      = t;           % This expression can be revised for vectorization
d_hh     = t;
d_t      = (h + hh);

end

