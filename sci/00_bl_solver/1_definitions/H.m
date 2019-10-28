function [ h, h_dstr, h_t] = H( dstr , t)
%H      Shape Factor
%   Derivatives outputed in following order
% h         Shape Factor
% h_dstr    Shape Factor derivative to displacement thickness (mass)
% h_t       Shape Factor derivative to momentum thickness
h       =   dstr ./  t;
h_dstr  =   1    ./  t;
h_t     = - dstr ./ (t.^2);

end

