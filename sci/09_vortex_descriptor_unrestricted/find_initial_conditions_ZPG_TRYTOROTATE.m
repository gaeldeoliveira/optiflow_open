% Start from clean sheet
clear all; close all;
clc

% Add necessary stuff
addpath closure_relations
addpath definitions/

% Base Flow
fname_root_bf     = 'ZPG_yaw0_base_d';    % Filename Root
d_list_bf         =  [0 50];              % Stances (in x_over_h)

% Experimental Conditions
exp_cond.h     =  5.0e-3;                 % [m   ] height of vane
exp_cond.d     = 12.5e-3;                 % [m   ] distance between pair trailing edges
exp_cond.l     = 12.5e-3;                 % [m   ] vane chord lenght
exp_cond.D     = 30.0e-3;                 % [m   ] pair separation
exp_cond.AoA   = 18;                      % [m   ] pair separation
exp_cond.u_inf = 15.16   ;                % [m/s ] edge velocity
exp_cond.nu_inf= 1.461e-5;                % [m2 s-1]  Kinematic Viscosity (nu = mu / rho)
% Define Crossflow scaling variables
S              = exp_cond.D / 2;          % [m   ] half-width of vortex system cell
h              = exp_cond.h;              % [m   ] height of vane
x              = exp_cond.h*d_list_bf;    % [m   ] distance behing vane trailing edge
u_inf          = exp_cond.u_inf;          % [m/s ] edge velocity
nu_inf         = exp_cond.nu_inf;         % [m2 s-1]  Kinematic Viscosity (nu = mu / rho)
delta_ref      =  ... 0.0150 ;            % [m   ] ref boundary layer height (for gauge, need not be accurate!)
    delta_99(1.35,2467,u_inf/nu_inf);

% Now Read and Scale Data
% Read (no filtering needed here, because no vane reflections)
ED_bf             = BLT_experiment_data_reader(fname_root_bf, d_list_bf, h);
% Set Crossflow scaling variables
ED_bf.set_crossflow_scaling_variables(S, delta_ref);


% Now show results on 

target_y = linspace(0,2.1)';
z_ref    = -1;
u_bf_0h  = ED_bf.interp_velocities(0 , target_y,z_ref) / u_inf;
u_bf_50h = ED_bf.interp_velocities(50, target_y,z_ref) / u_inf; 
plot(u_bf_0h , target_y, u_bf_50h, target_y); grid on; hold on;
xlabel('Streamwise Speed'); ylabel('y/\delta - normal coordinate');

[dstr_m10, theta_m10, Hk_m10, target_y_m10, u_target_y_m10] = ED_bf.integral_parameters(0, -1.0, 2.1, u_inf);
[dstr_m05, theta_m05, Hk_m05, target_y_m05, u_target_y_m05] = ED_bf.integral_parameters(0, -0.5, 2.1, u_inf);
[dstr_c00, theta_c00, Hk_c00, target_y_c00, u_target_y_c00] = ED_bf.integral_parameters(0,  0.0, 2.1, u_inf);
[dstr_M05, theta_M05, Hk_M05, target_y_M05, u_target_y_M05] = ED_bf.integral_parameters(0,  0.5, 2.1, u_inf);
[dstr_M10, theta_M10, Hk_M10, target_y_M10, u_target_y_M10] = ED_bf.integral_parameters(0,  1.0, 2.1, u_inf);
plot([-1, -0.5, 0, 0.5, 1], [Hk_m10, Hk_m05, Hk_c00, Hk_M05, Hk_M10]); grid on

target_x = 0; 
target_z = linspace(-1,1);          % List of z/S stances
Hk_z     = zeros(size(target_z ));  % Shape factor at z/S stance
theta_z  = zeros(size(target_z ));  % Momentum thicknes at z/S stance (over delta_ref)
for n_z = 1:length(target_z)
    [~, theta_z(n_z), Hk_z(n_z)] = ED_bf.integral_parameters(0,  target_z(n_z), 2.0, u_inf);
end
rt_z     = theta_z*u_inf/nu_inf;    % Reynolds theta at z/S stance

% Take reference values as averages
rt_0    = mean(rt_z);
Hk_0    = mean(Hk_z);

figure(1)
subplot(321)
plot(target_z , rt_z); grid on; hold on;
ph1 = plot([-1 1], rt_0*[1 1], 'k:'); 
axis([-1.0000    1.0000  1800 2800])
xlabel('z/S - Spanwise position');
ylabel('Re_{\theta} - Mom. Reynolds')
legend(ph1, ['Reference Re_{\theta}^0=' , num2str(rt_0,4)], 'Location', 'South')
subplot(323)
plot(target_z , Hk_z); grid on; hold on;
ph2 = plot([-1 1], Hk_0*[1 1], 'k:'); 
xlabel('z/S - Spanwise position');
ylabel('H_k - Shape Factor')
legend(ph2 , ['Reference  H_k^0=' , num2str(Hk_0,3)], 'Location', 'South')
axis([-1.0000    1.0000    1.15 1.45])

% set(gcf , 'PaperType', 'A5');
% orient landscape;
% print ('-dpdf', ['initial_conditions_ZPG_subplot.pdf']) %#ok<NBRAK>


%% Now look at velocity profiles
[dstr_l, theta_l, Hk_l, target_y_l, u_target_y_l] = ED_bf.integral_parameters(0, -1.0, 2.0, u_inf);
[dstr_c, theta_c, Hk_c, target_y_c, u_target_y_c] = ED_bf.integral_parameters(0,  0.0, 2.0, u_inf);
[dstr_r, theta_r, Hk_r, target_y_r, u_target_y_r] = ED_bf.integral_parameters(0,  1.0, 2.0, u_inf);

% Make a Swafford profile
SP = swafford_profile();
SP.update_hk_rt_pair(Hk_0, rt_0);
u_target_y_c_swafford = SP.evaluate_profile(target_y_c * delta_ref/ theta_c);

% Plot Comparison
figure(2)
subplot(221)
exp_l = plot(u_target_y_l,          target_y_l, 'Color', [0 0.4470 0.7410]); hold on; grid on;
exp_c = plot(u_target_y_c,          target_y_c, 'Color', [0 0.4470 0.7410]); 
exp_r = plot(u_target_y_r,          target_y_r, 'Color', [0 0.4470 0.7410]); 
swaf  = plot(u_target_y_c_swafford, target_y_c, 'k'); 
axis([0 1.05 min(target_y_l) max(target_y_l)]);
xlabel('u/u_inf'); ylabel('y/\delta_{ref}');

subplot(222)
exp_l = semilogy(u_target_y_l,          target_y_l, 'Color', [0 0.4470 0.7410]); hold on; grid on;
exp_c = semilogy(u_target_y_c,          target_y_c, 'Color', [0 0.4470 0.7410]); 
exp_r = semilogy(u_target_y_r,          target_y_r, 'Color', [0 0.4470 0.7410]); 
swaf  = semilogy(u_target_y_c_swafford, target_y_c, 'k'); 
axis([0 1.05 min(target_y_l) max(target_y_l)]);
xlabel('u/u_inf'); ylabel('log(y/\delta_{ref})');
legend([exp_l , swaf],  'Experimental Baseflow', ['Swafford (H_k =' num2str(Hk_0,3) '  Re_{\theta} =' num2str(rt_0,4) ])

% set(gcf , 'PaperType', 'A5');
% orient landscape;
% print ('-dpdf', ['swafford_velocity_profile_ZPG.pdf']) %#ok<NBRAK>


%% Now try to make an offset for each velocity

% Make refined Swafford Profile
SP                      = swafford_profile();
SP.update_hk_rt_pair(Hk_0, rt_0);
target_y_hf             =  linspace(0,2,4000);
u_target_y_hf_swafford  = SP.evaluate_profile(target_y_hf * delta_ref/ theta_c);

% Define number of points on added profile
N_added_stances = 5;


% For left
% Find starting velocity of profile
u_target_y_l_start = u_target_y_l(1);
% Equivalent corresponding height on Swafford profile
y_l_start_offset = interp1(u_target_y_hf_swafford, target_y_hf, u_target_y_l_start);
% Make vector of added stances in y
y_l_added_stances = y_l_start_offset * linspace(0,1,N_added_stances)';
% Make vector of added stances in u
u_l_added_stances = interp1(target_y_hf, u_target_y_hf_swafford, y_l_added_stances);
% Make extended vector of stances in y 
y_l_extended      = [y_l_added_stances(1:end-1); target_y_l(:) + y_l_start_offset];
% Make extended vector of stances in u
u_l_extended      = [u_l_added_stances(1:end-1); u_target_y_l(:)];
% Now compare results
plot(u_l_extended, y_l_extended); hold on
plot(u_target_y_l, target_y_l  ); grid on


% For center
% Find starting velocity of profile
u_target_y_c_start = u_target_y_c(1);
% Equivalent corresponding height on Swafford profile
y_c_start_offset = interp1(u_target_y_hf_swafford, target_y_hf, u_target_y_c_start);
% Make vector of added stances in y
y_c_added_stances = y_c_start_offset * linspace(0,1,N_added_stances)';
% Make vector of added stances in u
u_c_added_stances = interp1(target_y_hf, u_target_y_hf_swafford, y_c_added_stances);
% Make extended vector of stances in y 
y_c_extended      = [y_c_added_stances(1:end-1); target_y_c(:) + y_c_start_offset];
% Make extended vector of stances in u
u_c_extended      = [u_c_added_stances(1:end-1); u_target_y_c(:)];
% Now compare results
plot(u_c_extended, y_c_extended); hold on
plot(u_target_y_c, target_y_c  ); grid on



% For right
% Find starting velocity of profile
u_target_y_r_start = u_target_y_r(1);
% Equivalent corresponding height on Swafford profile
y_r_start_offset = interp1(u_target_y_hf_swafford, target_y_hf, u_target_y_r_start);
% Make vector of added stances in y
y_r_added_stances = y_r_start_offset * linspace(0,1,N_added_stances)';
% Make vector of added stances in u
u_r_added_stances = interp1(target_y_hf, u_target_y_hf_swafford, y_r_added_stances);
% Make extended vector of stances in y 
y_r_extended      = [y_c_added_stances(1:end-1); target_y_r(:) + y_r_start_offset];
% Make extended vector of stances in u
u_r_extended      = [u_r_added_stances(1:end-1); u_target_y_r(:)];
% Now compare results
plot(u_r_extended, y_r_extended); hold on
plot(u_target_y_r, target_y_r  ); grid on



% Now plot comparison
figure(102)
subplot(221)
exp_l = plot(u_l_extended,        y_l_extended, 'Color', [0 0.4470 0.7410]); hold on; grid on;
exp_c = plot(u_c_extended,        y_c_extended, 'Color', [0 0.4470 0.7410]); 
exp_r = plot(u_r_extended,        y_r_extended, 'Color', [0 0.4470 0.7410]); 
swaf  = plot(u_target_y_c_swafford, target_y_c, 'k'); 
axis([0 1.05 min(target_y_l) max(target_y_l)]);
xlabel('u/u_inf'); ylabel('y/\delta_{ref}');

subplot(222)
exp_l = semilogy(u_l_extended,        y_l_extended, 'Color', [0 0.4470 0.7410]); hold on; grid on;
exp_c = semilogy(u_c_extended,        y_c_extended, 'Color', [0 0.4470 0.7410]); 
exp_r = semilogy(u_r_extended,        y_r_extended, 'Color', [0 0.4470 0.7410]); 
swaf  = semilogy(u_target_y_c_swafford, target_y_c, 'k'); 
axis([0 1.05 min(target_y_l) max(target_y_l)]);
xlabel('u/u_inf'); ylabel('log(y/\delta_{ref})');
legend([exp_l , swaf],  'Experimental Baseflow', ['Swafford (H_k =' num2str(Hk_0,3) '  Re_{\theta} =' num2str(rt_0,4) ])




figure(17)
plot(u_target_y_l,          target_y_l); hold on; grid on;
plot(u_target_y_c,          target_y_c); 
plot(u_target_y_r,          target_y_r); 

figure(85)
plot(u_l_extended, y_l_extended); hold on
plot(u_c_extended, y_c_extended); hold on
plot(u_r_extended, y_r_extended); hold on


dst_l_extended = trapz(y_l_extended, 1 - u_l_extended)
dst_c_extended = trapz(y_c_extended, 1 - u_c_extended)
dst_r_extended = trapz(y_r_extended, 1 - u_r_extended)

dst_l_original = trapz(target_y_l, 1 -u_target_y_l)
dst_c_original = trapz(target_y_c, 1 -u_target_y_c)
dst_r_original = trapz(target_y_r, 1 -u_target_y_r)


tet_l_extended = trapz(y_l_extended, u_l_extended .* (1 - u_l_extended))
tet_c_extended = trapz(y_c_extended, u_c_extended .* (1 - u_c_extended))
tet_r_extended = trapz(y_r_extended, u_r_extended .* (1 - u_r_extended))

tet_l_original = trapz(target_y_l, u_target_y_l .*(1 -u_target_y_l))
tet_c_original = trapz(target_y_c, u_target_y_c .*(1 -u_target_y_c))
tet_r_original = trapz(target_y_r, u_target_y_r .*(1 -u_target_y_r))


hk_l_extended = dst_l_extended / tet_l_extended
hk_c_extended = dst_c_extended / tet_c_extended
hk_r_extended = dst_r_extended / tet_r_extended

hk_l_original = dst_l_original / tet_l_original
hk_c_original = dst_c_original / tet_c_original
hk_r_original = dst_r_original / tet_r_original


% Now try to apply rotation trick to an entire "frame of the boundary layer"
% Input arguments 
u_inf
nu_inf

% Set baseflow ED (ED_bf) as current ED for code prototyping
ED = ED_bf;

% Select current frame (this will become a loop counter)
nd_list = 1; 

% Get boundary layer parameters for this frame
target_x = 0; 
target_z = linspace(-1,1);          % List of z/S stances
Hk_z     = zeros(size(target_z ));  % Shape factor at z/S stance
theta_z  = zeros(size(target_z ));  % Momentum thicknes at z/S stance (over delta_ref)
for n_z = 1:length(target_z)
    [~, theta_z(n_z), Hk_z(n_z)] = ED_bf.integral_parameters(ED.d_list(nd_list),  target_z(n_z), 2.0, u_inf);
end
rt_z     = theta_z*u_inf/nu_inf;    % Reynolds theta at z/S stance

% Take reference values as averages
rt_mean    = mean(rt_z);
Hk_mean    = mean(Hk_z);
theta_mean = mean(theta_z);

% Make high fidelity swafford profile
SP                      = swafford_profile();
SP.update_hk_rt_pair(Hk_mean, rt_mean);
target_y_hf             =  linspace(0,2,4000);
u_target_y_hf_swafford  = SP.evaluate_profile(target_y_hf * delta_ref/ theta_mean);

% Now get original u and y _meshes
y_mesh_original_scaled = ED.data_struct_cell{nd_list}.y / ED.delta_ref;
u_mesh_original_scaled = ED.data_struct_cell{nd_list}.u / u_inf;

% And allocate space for extended u mesh
u_mesh_extrapolated    = zeros(size(u_mesh_original_scaled));


% Now make an extrapolation for each row (this will become a nested loop)
for n_z_index = 1:size(y_mesh_original, 2)

% Extract current row
y_line_original_scaled = y_mesh_original(:,n_z_index);
u_line_original_scaled = u_mesh_original(:,n_z_index);

% Reorder for extrapolation
% [target_y_l , sort_index]  = sort(y_line_original_scaled);
% u_target_y_l               = u_line_original_scaled(sort_index);

% Reorder for extrapolation (including zero filtering)
% Filter
y_line_original_scaled_filtered = y_line_original_scaled(u_line_original_scaled>eps());
u_line_original_scaled_filtered = u_line_original_scaled(u_line_original_scaled>eps());
% Then reorder
[target_y_l , sort_index]  = sort(y_line_original_scaled_filtered);
u_target_y_l               = u_line_original_scaled_filtered(sort_index);

% Add struct cell
% Find starting velocity of profile
u_target_y_l_start = u_target_y_l(1);
% Equivalent corresponding height on Swafford profile
y_l_start_offset = interp1(u_target_y_hf_swafford, target_y_hf, u_target_y_l_start);
% Make vector of added stances in y
y_l_added_stances = y_l_start_offset * linspace(0,1,N_added_stances)';
% Make vector of added stances in u
u_l_added_stances = interp1(target_y_hf, u_target_y_hf_swafford, y_l_added_stances);
% Make extended vector of stances in y 
y_l_extended      = [y_l_added_stances(1:end-1); target_y_l(:) + y_l_start_offset];
% Make extended vector of stances in u
u_l_extended      = [u_l_added_stances(1:end-1); u_target_y_l(:)];
% Now compare results
plot(u_l_extended, y_l_extended); hold on
plot(u_target_y_l, target_y_l  ); grid on
% Now re-interpolate onto extrapolated u_mesh (unscaled)
u_mesh_extrapolated(:,n_z_index) = interp1(y_l_extended, u_l_extended, y_line_original_scaled) * u_inf;


end

% Feed result back into extrapolated mesh!
ED.data_struct_cell{nd_list}.u_mesh_extrapolated = u_mesh_extrapolated;

figure(20)
surf(ED.data_struct_cell{1}.z, ED.data_struct_cell{1}.y, u_mesh_extrapolated); view(2); shading flat
figure(30)
surf(ED.data_struct_cell{1}.z, ED.data_struct_cell{1}.y, ED.data_struct_cell{1}.u); view(2); shading flat
