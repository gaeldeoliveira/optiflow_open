function [u_x, u_y, u_mag_xy, u_r, u_z, u_mag] = vortex_ring_induction_alfa3_fun(Gamma, R, delta, ellipke_tol, t_x, t_r)

% Restate induction Point (P) in reference paper notation
% (embedded in XZ plane, coordinates given as P = [r, 0, z])
r     = abs(t_y);        % [m       ] - (vectorizable, tier 1) Corresponds to x coordinate of induction point (vectorizable)
z     = t_x     ;        % [m       ] - (vectorizable, tier 1) Corresponds to z coordinate of induction point (vectorizable)

% Compute Intermediate variables
A     =      (r-R).^2 + z.^2 + delta.^2 ;       % Big A from paper   (eq. 6)
a     = sqrt((r+R).^2 + z.^2 + delta.^2);       % Small a from paper (eq. 6)
m     = 4.*    r.*R ./ a.^2;                    % Small m from paper (eq. 6)
% Modify m definition to allow for negative x/r doesn't work on its own, 

% Compute Elliptic Integrals of the first (K of m) and second (E of m) kind
% (do this once and reuse result throughout expression, costly operation)
[K,E] = ellipke( m , ellipke_tol);

% Proceed to vortex ring induction, as given in reference paper notation
% Radial induction, equation 4
u_r   = Gamma ./ (2 * pi() .* a) .* (z./r) .* (   (r.^2 + R.^2 + z.^2 + delta.^2) .* (E./A) - K);
% Axial induction, equation 5
u_z   = Gamma ./ (2 * pi() .* a) .*           ( - (r.^2 - R.^2 + z.^2 + delta.^2) .* (E./A) + K);
% Velocity magnitude
u_mag = sqrt(u_r.^2 + u_z.^2);

% Restate velocity into (x_hat,y_hat) plane coordinates of external flow analysis
u_x   = u_z;
u_y   = u_r .* sign(t_y);
u_mag_xy = sqrt(u_x.^2 + u_y.^2);


end