function [u_z, u_r, u_mag] = vortex_ring_induction_alfa6_fun(Gamma, Zr, Rr, delta, ellipke_tol, Zt, Rt)

% Restate induction Point (P) in reference paper notation
% (embedded in XZ plane, coordinates given as P = [r, 0, z])
Rt0    = Rt;             % [m       ] - (vectorizable, tier 1) Corresponds to x coordinate of induction point (vectorizable)
Rt     = abs(Rt);        % [m       ] - (vectorizable, tier 1) Corresponds to x coordinate of induction point (vectorizable)
Zt     = Zt - Zr;        % [m       ] - (vectorizable, tier 1) Corresponds to z coordinate of induction point (vectorizable)

if Rr > 0
    % Compute Intermediate variables
    A     =      (Rt-Rr).^2 + Zt.^2 + delta.^2 ;       % Big A from paper   (eq. 6)
    a     = sqrt((Rt+Rr).^2 + Zt.^2 + delta.^2);       % Small a from paper (eq. 6)
    %a     = sqrt((Rt+Rr).^2 + Zt.^2);       % Small a from paper (eq. 6)
    m     = 4.*    Rt.*Rr ./ a.^2;                    % Small m from paper (eq. 6)
    % Modify m definition to allow for negative x/r doesn't work on its own,
    
    % Compute Elliptic Integrals of the first (K of m) and second (E of m) kind
    % (do this once and reuse result throughout expression, costly operation)
    try
        [K,E] = ellipke( m , ellipke_tol);
    catch
        disp('Error: m values innapropriate');
    end
    
    % Proceed to vortex ring induction, as given in reference paper notation
    % Radial induction, equation 4
    u_r   = Gamma ./ (2 * pi() .* a) .* (Zt./Rt) .* (   (Rt.^2 + Rr.^2 + Zt.^2 + delta.^2) .* (E./A) - K);
    % Axial induction, equation 5
    u_z   = Gamma ./ (2 * pi() .* a) .*           ( - (Rt.^2 - Rr.^2 + Zt.^2 + delta.^2) .* (E./A) + K);
    % Velocity magnitude
    u_mag = sqrt(u_r.^2 + u_z.^2);
    
    % Restate velocity into (x_hat,y_hat) plane coordinates of external flow analysis
    u_r   = u_r .* sign(Rt0);
else
    % Ignore 
    u_z     = zeros(size(Zt + Rt));
    u_r     = zeros(size(Zt + Rt));
    u_mag   = zeros(size(Zt + Rt));
end