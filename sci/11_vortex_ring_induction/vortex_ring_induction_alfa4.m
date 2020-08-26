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

% Set Numerical Accuracy
ellipke_tol = eps();    % Numerical tolerance of elliptic integral evaluation (default is eps, higher tolerance/lower precision is faster!)

% Define Vortex Ring
% (centered on origin, embedded in XY plane)
Gamma = 1;              % [s^-1 m^3] - (vectorizable, tier 2) Circulation of Vortex Ring (not completely sure about units!
R     = 1;              % [m       ] - (vectorizable, tier 2) Radius of Vortex Ring 
delta = 0.1;           % [m       ] - (vectorizable, tier 2) Regularization parameter related to smoothing core/kernel size!

% Define target point (T = [t_x, t_y]) in (x_hat,y_hat) plane coordinates of external flow analysis
t_x   = 2;
t_y   = 2;

% Compute Induction
[u_x, u_y, u_r, u_z, u_mag] = vortex_ring_induction_alfa4_fun(Gamma, R, delta, ellipke_tol, t_x, t_y);

% Display Results for diagnosis!
display(['z/R       : ' num2str(t_x/R)]);
display(['r/R       : ' num2str(t_y/R)]);
display(['u_z       : ' num2str(u_z)]);
display(['u_r       : ' num2str(u_r)]);
display(['u_mag     : ' num2str(u_mag)]);
display([' ']); %#ok<NBRAK,DISPLAYPROG>
display(['x         : ' num2str(t_x)]);
display(['y         : ' num2str(t_y)]);
display(['u_x       : ' num2str(u_x)]);
display(['u_y       : ' num2str(u_y)]);

% Now make a mesh
n_points = 1000;
[t_x_mesh , t_y_mesh] = meshgrid(linspace(-abs(t_x),abs(t_x), n_points), linspace(-abs(t_y),abs(t_y), n_points));
% And compute induction over it
[u_x_mesh, u_y_mesh, u_r_mesh, u_z_mesh, u_mag_mesh] = vortex_ring_induction_alfa4_fun(Gamma, R, delta, ellipke_tol, t_x_mesh, t_y_mesh);
% Now plot
figure(1)
surf(t_x_mesh , t_y_mesh, u_mag_mesh)
view(2); shading flat;


% Or plot in (enlarged) quivers 
figure(2)
quiver(t_x_mesh, t_y_mesh, u_x_mesh, u_y_mesh, 10);





