% Script to study evolution of vortex pair (with a heuristically
% parabolized system!)
% Sames as alfa_1, but coordinate system is changed!

%addpath induction_functions/
%addpath data_functions/


close all; clear all; clc

% Define Vortex System Invariant Characteristics
n_cells = 5;
S       = 1.0;


y_v0    = 0.4;
eta0    = 0.2;
gamma_v = 0.2;
sigma_v = 0.1;

w_conv  = 1.0;  % Free stream (z,w)
dx      = 0.1;

% Initial Conditions
x_v = 0   ;     % x is aligned with free-stream and wall
y_v = y_v0;     % y is normal to wall and free-stream
z_v = eta0;     % z is aligned with wall and normal to free-stream

figure(1)
plot3( 0,  y_v0 ,  eta0, 'xr'); hold on; grid on; axis equal;
plot3( 0,  y_v0 , -eta0, 'xr');
figure(2)
plot(  eta0, y_v0 , 'xr'); hold on; grid on; axis equal;
plot( -eta0, y_v0 , 'xr'); 

% Construct Vortex Descriptor Object
VD = vortex_descriptor(S, n_cells);

% Start of Hypothetical Cycle (Marching Solution)
for n = 1:100
[w_v, v_v] = VD.induction_replicated_duplicated_vortex_system(z_v , y_v, z_v, y_v, gamma_v, sigma_v);
dt = dx / w_conv;

x_v = x_v + dx;
y_v = y_v + dt * v_v;
z_v = z_v + dt * w_v;

figure(1)
plot3( x_v, y_v, z_v, 'x'); hold on; grid on; axis equal;
plot3( x_v, y_v,-z_v, 'x');
xlabel('x'); ylabel('y'); zlabel('z');
axis([0 x_v 0 2 -S S])
figure(2)
plot(  z_v, y_v, 'x');
plot( -z_v, y_v, 'x');
ylabel('y'); xlabel('z');
axis([-S S 0 2])

end


figure(1)
% Make Plotting Ranges
n_range = 500;
% x_range = ((n_cells+1)*2*S) * linspace(-1, 1, n_range );
% y_range = linspace(-1, 1, n_range );

z_range = S * linspace(-1, 1, n_range );
y_range = linspace(0, 2, n_range );

% Make Plotting Mesh
[z_mesh, y_mesh] = meshgrid(z_range, y_range);

x_mesh = x_v*ones(size(z_mesh));

[w_mesh, v_mesh] = VD.induction_replicated_duplicated_vortex_system(z_mesh, y_mesh, z_v, y_v, gamma_v, sigma_v);

mag_mesh = sqrt(w_mesh.^2 + v_mesh.^2);

surf(x_mesh, y_mesh, z_mesh, mag_mesh)
shading interp

% Start of Hypothetical Cycle (Marching Solution)
for n = 1:100
[w_v, v_v] = VD.induction_replicated_duplicated_vortex_system(z_v , y_v, z_v, y_v, gamma_v, sigma_v);
dt = dx / w_conv;

x_v = x_v + dx;
y_v = y_v + dt * v_v;
z_v = z_v + dt * w_v;

figure(1)
plot3( x_v, y_v, z_v, 'x'); hold on; grid on; axis equal;
plot3( x_v, y_v,-z_v, 'x');
xlabel('x'); ylabel('y'); zlabel('z');
axis([0 x_v 0 2 -S S])

figure(2)
plot(  z_v, y_v, 'x');
plot( -z_v, y_v, 'x');
ylabel('y'); xlabel('z');
axis([-S S 0 2])

end


figure(1)
% Make Plotting Ranges
n_range = 500;
% x_range = ((n_cells+1)*2*S) * linspace(-1, 1, n_range );
% y_range = linspace(-1, 1, n_range );

z_range = S * linspace(-1, 1, n_range );
y_range = linspace(0, 2, n_range );

% Make Plotting Mesh
[z_mesh, y_mesh] = meshgrid(z_range, y_range);

x_mesh = x_v*ones(size(z_mesh));

[w_mesh, v_mesh] = VD.induction_replicated_duplicated_vortex_system(z_mesh, y_mesh, z_v, y_v, gamma_v, sigma_v);

mag_mesh = sqrt(w_mesh.^2 + v_mesh.^2);

surf(x_mesh, y_mesh, z_mesh, mag_mesh)
shading interp

% view([6.3963   -0.0000    0.7687   -3.2267; ...
%       2.7013    0.9362   -0.2248   -1.6117; ...
%       7.1966   -0.3514   -0.5988   46.1364; ...
%           0         0         0    1.0000])