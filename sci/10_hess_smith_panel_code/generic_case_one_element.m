% A simple modular Panel Code based on a free interpretation of the
% Hess-Smith method in velocity components
% Gael de Oliveira
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%       ActiHS  :   A simple modular Panel Code based on a free           %
%                   interpretation of the Hess-Smith method in velocity   %
%                   components                                            %
%                                                                         %
%       Usage   :   Standalone with script, for sail optimization, within %
%                   the kirikou-dogoro actuator codes or other codes and  %
%                   derivatives from the author                           %
%                                                                         %
%       Date    :   April 2011 to March 2017                              %
%       Author  :   Gael de Oliveira                                      %
%                                                                         %
%       License :   MIT, as the rest of this repository 		    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Start with a clean environment
clear all; close all; clc;

%% Inputs
airfoil_file = 'airfoils/naca0012_160.air';  % Airfoil shape
px_c         = 1/4;                 % Airfoil Rotation Point (x coordinate, in airfoil units)
py_c         = 0;                   % Airfoil Rotation Point (y coordinate, in airfoil units)
alpha_geo    = 10*pi()/180;         % Airfoil Rotation Angle (in radians, corresponds to angle of attack if free-stream set to 0)
% Inflow: Straight Free-Stream Definition
u_inf        = 1;                   % Free-Stream Magnitude
alpha        = 0*pi/180;            % Set inflow angle of attack in radians
% Inflow: Still air Rotation Definition (only relevant for VAWT airfoils, set rotation to false otherwise)
c_over_R     = 0.1;                 % Now specify rotational effects
x_cp         = 1/4;                 % Chordwise coordinate of point on longest chord line that lies on radius of rotation around mean point
rotation     = false;               % Is effect of rotation accounted for ? ()

%% Preprocessing
% Load Airfoil Geometry
coord        = load(airfoil_file);

% Extract Airfoil Coordinates
px           = coord(:,1);
py           = coord(:,2);
% Rotate Airfoil around (px_c,py_c) point with alpha_rot angle
px_rotated   =   cos(alpha_geo) * (px-px_c) + sin(alpha_geo) * (py - py_c) + px_c;
py_rotated   = - sin(alpha_geo) * (px-px_c) + cos(alpha_geo) * (py - py_c) + py_c;

%% Create Case and Solve it!
% Create Case
ipc = inviscid_panel_case(px_rotated, py_rotated);      % Create a new inviscid_panel_case object instance, holding all the data and controlling execution flow
% Describe Straight Inflow
ipc.u_inf    = u_inf;
ipc.alpha    = alpha;
% Describe Still-Air Rotation
ipc.c_over_R = c_over_R;
ipc.x_cp     = x_cp;
ipc.rotation = rotation;
% Force, Solve and Postprocess
ipc.generate_solution;

%% Plot Results

% Plot Geometry
figure(1)
plot(px        , py         , 'o-'); hold on;
plot(px_rotated, py_rotated , 'o-');
xlabel('x/c'); ylabel('y/c')
grid on; axis equal;

% Plot results!
figure(1)
plot(ipc.px_middle , ipc.cp_plot)
grid on;
xlabel('x/c'); ylabel('Cp')
title('Inviscid Cp Plot')

% Forces are computed and available in ipc fields!
disp(['AoA (Geometric + Inflow  ) : ' num2str((alpha + alpha_geo)*180/pi()) ' deg'])
disp(['Cl  (Pressure Integration) : ' num2str(ipc.L_Cp)])
disp(['Cl  (Total Circulation   ) : ' num2str(ipc.L_gamma)])
disp(['Cl  (Theoretical         ) : ' num2str(2*pi()*sin(alpha + alpha_geo)+0.09*(max(py)-min(py)))])

%% Now plot velocity field
tic
% Define Velocity Field Generation Window Bounds
% Bounds
x_min = -0.25; x_max = 1.25;
y_min = -0.50; y_max = 0.50;
% Resolution
x_res = 120;
y_res =  80;
% Ranges
x_range = linspace(x_min, x_max, x_res);
y_range = linspace(y_min, y_max, y_res);
             
% % Make Mesh For Point Computation
[x_mesh,y_mesh] = meshgrid(x_range, y_range);
             
% % Compute Induced Speeds on Mesh Points 
[u_induced_mesh , v_induced_mesh] = ipc.induced_speed_on_many_points(x_mesh,y_mesh);           
% Add Free-Stream to Induced Speeds
u_mesh = u_induced_mesh + u_inf;
v_mesh = v_induced_mesh;
% And compute velocity norm
u_norm_mesh = sqrt(u_mesh.^2 + v_mesh.^2);
toc
% 
figure(3)
surf(x_mesh, y_mesh, u_norm_mesh); view(2); shading flat;
xlabel('x/c'); ylabel('y/c'); colorbar;
title('Velocity Magnitude')
figure(4)
surf(x_mesh, y_mesh, u_induced_mesh); view(2); shading flat;
xlabel('x/c'); ylabel('y/c'); colorbar;
title('Velocity Induced in x-Direction')

figure(5)
surf(x_mesh, y_mesh, 0.5*u_norm_mesh.^2); view(2); shading flat;
xlabel('x/c'); ylabel('y/c'); colorbar;
title('Pressure Coefficient')

