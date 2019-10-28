% Declare Variables
syms U_inf nu           % External Flow Speed, Kinematic Viscosity
syms x                  % Longitudinal Coordinate
syms y                  % Normal Coordinate (orthonormal to x)
syms f f1 f2 f3         % Define Blasius function and its eta derivatives


% Define Normal Scale

% Define scaled normal coordinate
eta         = sqrt(U_inf / nu) * y / sqrt(x);

% Compute First Derivatives
deta_dx     = diff(eta, x);
deta_dy     = diff(eta, y);

% Define Speed in terms of Blasius function
u = U_inf * f1;

% Define 






