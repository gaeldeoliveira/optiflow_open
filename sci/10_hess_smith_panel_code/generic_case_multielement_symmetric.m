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
% Prediction: Single Element Code on Naca 0012
%AoA (Geometric + Inflow  ) : 10 deg
%Cl  (Pressure Integration) : 1.1063
%Cl  (Total Circulation   ) : -1.1052
%Cl  (Theoretical         ) : 1.1019
%Elapsed time is 43.915679 seconds.
%
%                      L_Cp: [1.1162 -1.1162]
%                   L_gamma: [-1.1158 1.1158]
%                      D_Cp: [-0.0026 -0.0026]
%
%                      L_Cp: [-1.1172 1.1172]
%                   L_gamma: [1.1163 -1.1163]
%                      D_Cp: [-0.0026 -0.0026]
%

%% Start with a clean environment
clear all; close all; clc; %#ok<CLALL>

%% Inputs
% Geometry
airfoil_file = 'airfoils/naca0012.air';  % Filename: Xfoil Plain Format (counterclockwise ordering, TE->LE->LE, adimensionalized by chord)
px_c         = 1/4;             % Airfoil Rotation Point (x coordinate, in airfoil units)
py_c         = 0;               % Airfoil Rotation Point (y coordinate, in airfoil units)
x_duct       = 0;               % Translation of airfoil rotation point from origin (x direction)
y_duct       = 2.0;           % Translation of airfoil rotation point from origin (y direction)
alpha_duct   = 10*pi()/180;     % Duct Airfoils Geometric Angle of Attack (in radians)
% Inflow: Straight Free-Stream Definition
u_inf        = 1;               % Free-Stream Magnitude
alpha        = 0*pi/180;        % Set angle of attack in radians
% Inflow: Still air Rotation Definition (only relevant for VAWT airfoils, set rotation to false otherwise)
c_over_R     = 0.1;             % Now specify rotational effects
x_cp         = 1/4;             % Chordwise coordinate of point on longest chord line that lies on radius of rotation around mean point
rotation     = false;           % Is effect of rotation accounted for ? ()
% Postprocessing Options
graph_viz    = false;           % Make a beautiful picture of the velocity field? (setting to true takes and awful lot of time!)


%%  Preprocessing
% Load Airfoil Geometry
coord        = load(airfoil_file);
% Extract single airfoil coordinates
px           = coord(:,1);
py           = coord(:,2);

% Make Duct Top Airfoil
%   Rotate    (+alpha around px_c,py_c)
%   Translate (to x_duct,+y_duct)
px_top       =   cos(alpha_duct) * (px-px_c) + sin(alpha_duct) * (py - py_c) + px_c + x_duct;
py_top       = - sin(alpha_duct) * (px-px_c) + cos(alpha_duct) * (py - py_c) + py_c + y_duct;
% Make Duct Bottom Airfoil
%   Reflect   (-py over chordline)
%   Rotate    (-alpha around px_c,py_c)
%   Translate (to x_duct,-y_duct)
px_bot       =   cos(-alpha_duct) * (px-px_c) + sin(-alpha_duct) * (-(py - py_c)) +   px_c  + x_duct;
py_bot       = - sin(-alpha_duct) * (px-px_c) + cos(-alpha_duct) * (-(py - py_c)) + (-py_c) - y_duct;
%   Reorder   (to maintain counterclockwise ordering after horizontal reflection)
px_bot       = px_bot(end:-1:1);
py_bot       = py_bot(end:-1:1);

% Make Coordinates for Duct Geometry from airfoil geometry
px_duct      = [px_top ; px_bot];
py_duct      = [py_top ; py_bot];

% Make Index of Trailing Edge Start and Ending points
i_start      = [1 ;length(px_top)+1];
i_end        = [length(px_top) ;length(px_top)+length(px_bot)];

%% Case Creation, Definition and Solution
% Create Case Object
ipc          = inviscid_panel_case_multielement(px_duct , py_duct, i_start , i_end);
% Describe Straight Inflow
ipc.u_inf    = u_inf;
ipc.alpha    = alpha;
% Describe Still-Air Rotation
ipc.c_over_R = c_over_R;
ipc.x_cp     = x_cp;
ipc.rotation = rotation;
% Force, Solve and Postprocess
ipc.generate_solution;

%% Extract Results By Element
% First Element (Top)
n_element = 1;
element(n_element).indices   = ipc.te_first_indices(n_element):ipc.te_last_indices(n_element);
element(n_element).px_middle = ipc.px_middle(element(n_element).indices);
element(n_element).py_middle = ipc.py_middle(element(n_element).indices);
element(n_element).cp_plot   = ipc.cp_plot(  element(n_element).indices);
% Second Element (Bottom)
n_element = 2;
element(n_element).indices   = ipc.te_first_indices(n_element):ipc.te_last_indices(n_element);
element(n_element).px_middle = ipc.px_middle(element(n_element).indices);
element(n_element).py_middle = ipc.py_middle(element(n_element).indices);
element(n_element).cp_plot   = ipc.cp_plot(  element(n_element).indices);

%% Plot
% Plot Geometry!
figure(1)
plot(element(1).px_middle , element(1).py_middle, 'o-'); hold on;
plot(element(2).px_middle , element(2).py_middle, '^-'); hold on;
plot(px_duct, py_duct, 'xk')
grid on; axis equal;
xlabel('x/c'); ylabel('y/c');
legend('Element 1 (panel centers)', 'Element 2 (panel centers)', 'Panel Edges');

% Plot Cp Distribution!
figure(2)
plot(element(1).px_middle , element(1).cp_plot, 'o-'); hold on;
plot(element(2).px_middle , element(2).cp_plot, '^-'); 
grid on;
xlabel('x/c'); ylabel('Cp');
title('Inviscid Cp Plot')
legend('Element 1', 'Element 2');


%% Plot streamwise along centerline
% Define Centerline
x_centerline = linspace(-2,2);
y_centerline = zeros(size(x_centerline));
% Compute Induced Speeds on Centerline
[u_induced_centerline , v_induced_centerline] = ipc.induced_speed_on_many_points(x_centerline,y_centerline);
% Add Free-Stream to Induced Speeds
u_centerline = u_induced_centerline + u_inf;
v_centerline = v_induced_centerline;
figure(5)
subplot(211)
plot(x_centerline, u_centerline)
xlabel('x/c'); ylabel('u/u_inf'); grid on;
title('Streamwise Velocity along Centerline')
subplot(212)
plot(x_centerline, v_centerline)
xlabel('x/c'); ylabel('v/u_inf'); grid on;
title('Transverse Velocity along Centerline')


%% Now plot velocity field
if graph_viz == true
    tic
    % Define Velocity Field Generation Window Bounds
    % Bounds
    x_min = -0.25; x_max = 1.25;
    y_min = -1; y_max = 1;
    % Resolution
    x_res = 480;
    y_res = 320;
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
end

