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
delta = 0.01;           % [m       ] - (vectorizable, tier 2) Regularization parameter related to smoothing core/kernel size!

% Define target point (T = [t_x, t_y]) in (x_hat,y_hat) plane coordinates of external flow analysis
t_x   = 1;
t_y   = -1;

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
[K,E] = ellipke( m , eps );

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

% Display Results for diagnosis!
display(['z/R       : ' num2str(z/R)]);
display(['r/R       : ' num2str(r/R)]);
display(['u_z       : ' num2str(u_z)]);
display(['u_r       : ' num2str(u_r)]);
display(['u_mag     : ' num2str(u_mag)]);
display([' ']); %#ok<NBRAK,DISPLAYPROG>
display(['x         : ' num2str(t_x)]);
display(['y         : ' num2str(t_y)]);
display(['u_x       : ' num2str(u_x)]);
display(['u_y       : ' num2str(u_y)]);
display(['u_mag_xy  : ' num2str(u_mag_xy)]);


% Plot Elliptic Integrals
tic; [K,E] = ellipke(linspace(0,1),eps) ; toc;
figure(1); plot(linspace(0,1), K); grid on; ylabel('K'); figure(2); plot(linspace(0,1), E); grid on; ylabel('E');










