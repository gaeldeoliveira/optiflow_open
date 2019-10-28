% Clean WorkSpace
clear all; close all; clc;

% Define Filename
nd          = 10;
data_folder = 'data/';
fname_root  = 'ZPG_yaw0_VGrect_d';
fname_ext   = '.dat';
if nd > 0
    filename    = [data_folder  fname_root, num2str(nd), fname_ext ];
else
     filename    = [data_folder  fname_root(1:end-1), 'o', fname_ext ];
end

% Experimental Conditions
exp_cond.h     =  5.0e-3;                 % [m   ] height of vane
exp_cond.d     = 12.5e-3;                 % [m   ] distance between pair trailing edges
exp_cond.l     = 12.5e-3;                 % [m   ] vane chord lenght
exp_cond.D     = 30.0e-3;                 % [m   ] pair separation
exp_cond.AoA   = 18;                      % [m   ] pair separation
exp_cond.u_inf = 15.0   ;                 % [m/s ] edge velocity

case_cond.S = exp_cond.D / 2;     % [m   ] half-width of vortex system cell
case_cond.h = exp_cond.h;         % [m   ] height of vane
case_cond.x = exp_cond.h*nd;      % [m   ] distance behing vane trailing edge

% Read data
data_struct = read_dat_file_into_data_struct(filename);

% Now, try plotting something out!
figure(1)
surf(data_struct.z , data_struct.y, data_struct.u); 
view(2); shading flat;
xlabel('x'); ylabel('y'); colorbar();
title('Streamwise Velocity')

% Perturbation Component
figure(2)
u_pert = zeros(size(data_struct.u));
for i_line = 1:data_struct.I
    u_pert(i_line,:) = data_struct.u(i_line,:) - mean(data_struct.u(i_line,:));
end
surf(data_struct.z , data_struct.y, -u_pert); 
view(2); shading flat;
xlabel('x'); ylabel('y'); colorbar();
title('Streamwise Velocity Perturbation');

figure(3)
[curl_yz cav]= curl(data_struct.z , data_struct.y, data_struct.w , data_struct.v);
pcolor(data_struct.z , data_struct.y,curl_yz); 
shading interp; hold on; 
quiver(data_struct.z , data_struct.y, data_struct.w , data_struct.v)
title('Vorticity')

figure(4)
surf(data_struct.z , data_struct.y, data_struct.v); 
view(2); shading flat;
xlabel('x'); ylabel('y'); colorbar();
title('Y Velocity')

% Which means we could reshape all arrays from the start, in principle! The
% data is already pregridded! Youppi! vz corresponds to the streamwise
% component of the velocity


% Now reinterpolate to get a velocity profile!
target_y = linspace(min(min(data_struct.y)), max(max(data_struct.y)));
target_z = 0; %exp_cond.D/2;% * ones(size(target_y));

% This gives an error
%target_vz = interp2(data_struct.x, data_struct.y, data_struct.vz, target_x , target_y );
% But this works (so we need to transpose the data, not important for surf (as of today), but important for interp2 (ah, meshgrid subtleties!) )
target_u = interp2(data_struct.z, data_struct.y, data_struct.u, target_z , target_y );

% Now plot
figure(5)
subplot(121)
plot(target_u , target_y)
ylabel('y'); zlabel('vz'); grid on;
% And check law of the wall
subplot(122)
plot(target_u , log10(target_y))
ylabel('log_{10}(y)'); zlabel('vz'); grid on;


save([filename(1:end-3),'mat'], 'data_struct')

% 
% d5
z_min = -0.007549;
z_max =  0.00944;
y_min =  0.0134; 
y_max =  0.02756;
% d6
z_min = -0.007454;
z_max =  0.009157;
y_min =  0.01708; 
y_max =  0.03152;
% d7
z_min = -0.007171;
z_max =  0.009723;
y_min =  0.02152; 
y_max =  0.03388;
% d8
z_min = -0.007454;
z_max =  0.009629;
y_min =  0.0252; 
y_max =  0.03388;
% d9
z_min = -0.007738;
z_max =  0.009912;
y_min =  0.02973; 
y_max =  0.03388;
% d10
z_min = -0.007738;
z_max =  0.009912;
y_min =  0.02973; 
y_max =  0.03388;

%
z_min_list = [-0.007549, -0.007454, -0.007171, -0.007454, -0.007738, -0.007738];
z_max_list = [ 0.00944 ,  0.009157,  0.009723,  0.009629,  0.009912,  0.009912];
y_min_list = [ 0.0134  ,  0.01708 ,  0.02152 ,  0.0252  ,  0.02973 ,  0.02973 ];
y_max_list = [ 0.02756 ,  0.03152 ,  0.03388 ,  0.03388 ,  0.03388 ,  0.03388 ];

fieldname = 'u';
data_struct_f = filter_data_struct_field( data_struct, fieldname, z_min, z_max, y_min, y_max);


figure(6)
surf(data_struct.z , data_struct.y, data_struct_f.u); 
view(2); shading flat;
xlabel('x'); ylabel('y'); colorbar();
title('Streamwise Velocity (Filtered)')
