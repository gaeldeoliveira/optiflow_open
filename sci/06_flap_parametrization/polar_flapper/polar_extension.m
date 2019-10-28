% Airfoil Coefficients at large angle of attack
%
%   Based on: ftp://ftp.ecn.nl/pub/www/library/report/2001/rx01004.pdf

AReff         = 10;             % Effective Aspect Ratio
alfa          = 5*pi/180;       % Angle of attack
r_nose_over_c = .2;             % Nose radius over chord
Cd90          =        

% Normal force coefficient
Cn = Cd90 * ( 1/(.56 + .44*sin(alfa)) - .41*(1-exp(-17/AReff)) ) * sin(alfa);
% Virtual Angle
gamma = 0.28 * sqrt(r_nose_over_c);
% Tangential force coefficient
Ct = 1/2 * 0.0075 * cos(alfa) + Cn * sin(gamma);

% Lift Coefficient
Cl = Cn * cos(alfa) - Ct * sin(alfa);
% Drag Coefficient
Cd = Cn * sin(alfa) + Ct * cos(alfa);