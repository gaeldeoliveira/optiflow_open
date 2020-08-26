% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                         %
% Implement Vortex Ring Induction function described in Torque 2012:      %
%                                                                         %
%   de Vaal JB, Hansen MOL, Moan T. Validation of a vortex ring model     %
%   suited for aeroelastic simulations of floating wind turbines. Journal %
%   of Physics: Conference Series 555 (2014) 012025                       %
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                         %
% We are focusing on the case in which:                                   %
%   The vorticity vector is tangent to the vortex ring                    %
% It is assumed that:                                                     %
%   An orthonormal cartesian frame is used (XYZ)                          %
%   The vortex ring is centered on the origin and embedded in the XY plane%
%   Induced velocities are calculated at a point belonging to the XZ plane%
% Without loss of generality because:                                     %
%   These vortex rings are axisymmetric, so must be their induction.      %
% To compute induction with arbitrary orientations and targets, notice:   %
%   This kind of vortex ring induces no swirl, because the vorticity      %
%   vector is fully contained in the xy plane. In other words there are   %
%   no flow orbits around the z axis.                                     %
% Induction must threrefore be contained in a plane that contains both:   %
%   The vector between the ring center and the induction point            %
%   The vector normal to the vortex ring center                           %
%                                                                         %
% Purely swirling vortex rings are axisymmetric too but their induction is%
% not contained is the same plane. But keep in mind that non-axisymmetric %
% vortex rings also exist, its all about the orientation of the vorticity %
% vector.                                                                 %
%                                                                         %
% The code uses case-sensitive variable names. I understand this is       %
% ambiguous but it is compact and matches the notation of the ref paper.  %
%                                                                         %
%   Gael de Oliveira, August 2017                                         %
%                                                                         %
% % % % % % % % % % % % % % % % % % % %% % % % % % % % % % % % % % % % % % %

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                         %
% Stability challenges :                                                  %
%   Elliptic integrals are defined the range 0 <= M <= 1. Effectively this%
%   restricts applicability to x,r > 0. But this limitation can be        %
%   circumvented with a proper choice of vector orientations when making  %
%   the coordinate transformations for the general case!                  %
%                                                                         %
%   For M=1 we have K=Inf which is "sort" of an issue. I still don't grasp%
%   exactly where and when this happends, but it should be taken care of! %
%   A typical case is R=r with z=delta=0. Simple to show from definition  %
%   of m by reducing the equation m=1 to a polynomial. Other roots are    %
%   easy to assess and map for (z^2+delta^2) values.                      %
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Clean workspace
clear all; close all; clc; %#ok<CLALL>

% Set Numerical Accuracy
ellipke_tol = eps();    % Numerical tolerance of elliptic integral evaluation (default is eps, higher tolerance/lower precision is faster!)

% Define Vortex Ring
% (centered on origin, embedded in XY plane)
Gamma = 1;             % [s^-1 m^3] - (vectorizable, tier 2) Circulation of Vortex Ring (not completely sure about units!
Zr    = 1;             % [m       ] - (vectorizable, tier 2) Radius of Vortex Ring 
Rr    = 1;             % [m       ] - (vectorizable, tier 2) Radius of Vortex Ring 
delta = 0;          % [m       ] - (vectorizable, tier 2) Regularization parameter related to smoothing core/kernel size!
n_points = 100;        % [int.    ] - (sequential          ) Number of points in plotting mesh

% Define target point (T = [t_x, t_y]) in (x_hat,y_hat) plane coordinates of external flow analysis
Zt   = Zr;
Rt   = Rr;

% Compute Induction
[u_z, u_r, u_mag] = vortex_ring_induction_alfa5_fun(Gamma, Zr, Rr, delta, ellipke_tol, Zt, Rt);

% Display Results for diagnosis!
display(['z/R       : ' num2str(Zt/Rr)]);
display(['r/R       : ' num2str(Rt/Rr)]);
display(['u_z       : ' num2str(u_z)]);
display(['u_r       : ' num2str(u_r)]);
display(['u_mag     : ' num2str(u_mag)]);

% Now make a mesh
[t_Zt_mesh , t_Rt_mesh] = meshgrid(linspace(-abs(Zt),abs(Zt), n_points), linspace(-abs(Rt),abs(Rt), n_points));
% And compute induction over it
[u_z_mesh, u_r_mesh, u_mag_mesh] = vortex_ring_induction_alfa5_fun(Gamma, Zr, Rr, delta, ellipke_tol, t_Zt_mesh, t_Rt_mesh);
% Now plot
figure(1)
surf(t_Zt_mesh , t_Rt_mesh, u_mag_mesh)
view(2); shading flat;


% Or plot in (enlarged) quivers 
figure(2)
quiver(t_Zt_mesh, t_Rt_mesh, u_z_mesh, u_r_mesh, 10);

% Make plotting line accross ring (Zt sweep)
t = linspace(0, 1, 10*n_points+1);
t_Zt_range = Zr + 1*(t-0.05);
t_Rt_range = Rr + 0*(t);
[u_z_range, u_r_range, u_mag_range] = vortex_ring_induction_alfa5_fun(Gamma, Zr, Rr, delta, ellipke_tol, t_Zt_range, t_Rt_range);
figure(3)
plot(t, u_z_range, 'x-');  hold on;
plot(t, u_r_range, 'x-');

grid on

% Make plotting line accross ring (Rt sweep)
t = linspace(0, 1, 10*n_points+1);
t_Zt_range = Zr + 0*(t);
t_Rt_range = Rr + 1*(t-0.05);
[u_z_range, u_r_range, u_mag_range] = vortex_ring_induction_alfa5_fun(Gamma, Zr, Rr, delta, ellipke_tol, t_Zt_range, t_Rt_range);
figure(4)
plot(t, u_z_range, 'x-'); hold on;
plot(t, u_r_range, 'x-');
grid on

% Let us now look at the vortex segment
Rr1 =  1; Rr2 = 1;
Zr1 = -1; Zr2 = 1;
K = 128;

for k = 1:K
    % Make intermediate points
    Zr_k = Zr1 * (1-(k-1)/(K-1)) + Zr2 * (k-1)/(K-1);
    Rr_k = Rr1 * (1-(k-1)/(K-1)) + Rr2 * (k-1)/(K-1);
end

Zc = 0.5 * (Zr1 + Zr2);
Rc = 0.5 * (Rr1 + Rr2);

[u_z, u_r, u_mag] = vortex_line_induction_alfa5_fun(Gamma, Zr1, Rr1, Zr2, Rr2, delta, ellipke_tol, K, Zc, Rc);

% And now for some plotting
[t_Zt_mesh , t_Rt_mesh] = meshgrid(linspace(-1,5, n_points), linspace(-3,3, n_points));
% And compute induction over it
[u_z_mesh, u_r_mesh, u_mag_mesh] = vortex_line_induction_alfa5_fun(Gamma, Zr1, Rr1, Zr2, Rr2, delta, ellipke_tol, K, t_Zt_mesh, t_Rt_mesh);
% Now plot
figure(21)
surf(t_Zt_mesh , t_Rt_mesh, u_mag_mesh)
view(2); shading flat;
caxis([0 1]);


% Or plot in (enlarged) quivers 
figure(22)
[t_Zt_mesh , t_Rt_mesh] = meshgrid(linspace(-2,2, 50), linspace(-2,2, 50));
% And compute induction over it
K = 128;
[u_z_mesh, u_r_mesh, u_mag_mesh] = vortex_line_induction_alfa5_fun(Gamma, Zr1, Rr1, Zr2, Rr2, delta, ellipke_tol, K, t_Zt_mesh, t_Rt_mesh);
quiver(t_Zt_mesh, t_Rt_mesh, u_z_mesh, u_r_mesh, 1);


[t_Zt_mesh , t_Rt_mesh] = meshgrid(linspace(-2,2, 100), linspace(-2,2, 100));
K = 128; delta = 0;
[u_z_mesh_delta0, u_r_mesh_delta0, u_mag_mesh_delta0] = vortex_line_induction_alfa5_fun(Gamma, Zr1, Rr1, Zr2, Rr2, delta, ellipke_tol, K, t_Zt_mesh, t_Rt_mesh);
K = 128; delta = 1e-2;
[u_z_mesh_delta1, u_r_mesh_delta1, u_mag_mesh_delta1] = vortex_line_induction_alfa5_fun(Gamma, Zr1, Rr1, Zr2, Rr2, delta, ellipke_tol, K, t_Zt_mesh, t_Rt_mesh);
figure(31)
surf(t_Zt_mesh, t_Rt_mesh, u_r_mesh_delta1 - u_r_mesh_delta0); view(2);
figure(32)
surf(t_Zt_mesh, t_Rt_mesh, u_z_mesh_delta1 - u_z_mesh_delta0); view(2);



