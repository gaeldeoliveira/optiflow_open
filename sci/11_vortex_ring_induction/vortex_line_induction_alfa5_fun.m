function [u_z, u_r, u_mag] = vortex_line_induction_alfa5_fun(Gamma, Zr1, Rr1, Zr2, Rr2, delta, ellipke_tol, K, Zt, Rt)
% Make it over 4 points
% K =4;
% Move on to it over 16 points
% K =16;

% Initialize (for vectorization)
u_z   = zeros(size(Zt));
u_r   = zeros(size(Zt));
u_mag = zeros(size(Zt));


for k = 1:K
    % Make intermediate points
    Zr_k = Zr1 * (1-(k-1)/(K-1)) + Zr2 * (k-1)/(K-1);
    Rr_k = Rr1 * (1-(k-1)/(K-1)) + Rr2 * (k-1)/(K-1);
    % Compute induction from each intermediate point
    [u_z_k, u_r_k, u_mag_k] = vortex_ring_induction_alfa5_fun(Gamma/K, Zr_k, Rr_k, delta, ellipke_tol, Zt, Rt);
    % Add contribution of intermediate point to accumulator
    u_z = u_z + u_z_k;
    u_r = u_r + u_r_k;
end

% And compute velocity magnitude (vectorized)
u_mag = sqrt(u_z.^2 + u_r.^2);