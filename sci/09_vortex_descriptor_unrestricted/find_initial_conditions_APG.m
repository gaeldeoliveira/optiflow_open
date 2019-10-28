% Start from clean sheet
clear all; close all; clc; %#ok<CLALL>

% Add necessary stuff
addpath closure_relations
addpath definitions/

% Base Flow
fname_root_bf     = 'APG_yaw0_base_d';    % Filename Root
d_list_bf         =  [0 50];              % Stances (in x_over_h)

% Experimental Conditions
exp_cond.h     =  5.0e-3;                 % [m   ] height of vane
exp_cond.d     = 12.5e-3;                 % [m   ] distance between pair trailing edges
exp_cond.l     = 12.5e-3;                 % [m   ] vane chord lenght
exp_cond.D     = 30.0e-3;                 % [m   ] pair separation
exp_cond.AoA   = 18;                      % [m   ] pair separation
exp_cond.u_inf = 14.73   ;                % [m/s ] edge velocity
exp_cond.nu_inf= 1.461e-5;                % [m2 s-1]  Kinematic Viscosity (nu = mu / rho)
% Define Crossflow scaling variables
S              = exp_cond.D / 2;          % [m   ] half-width of vortex system cell
h              = exp_cond.h;              % [m   ] height of vane
x              = exp_cond.h*d_list_bf;    % [m   ] distance behing vane trailing edge
u_inf          = exp_cond.u_inf;          % [m/s ] edge velocity
nu_inf         = exp_cond.nu_inf;         % [m2 s-1]  Kinematic Viscosity (nu = mu / rho)
delta_ref      =  ... 0.0150 ;            % [m   ] ref boundary layer height (for gauge, need not be accurate! determined by manual iteration: 1st guess->feed new rt-> feed new rt-> until it is consistent enough!, then do the same with (hk,rt) pairs)
    delta_99(1.3567, 2851.6,u_inf/nu_inf); % delta_99(1.3567, 2851.9,u_inf/nu_inf);  %delta_99(1.3564,2854.7,u_inf/nu_inf); % delta_99(1.35,2854.9,u_inf/nu_inf); % delta_99(1.35,2857.6,u_inf/nu_inf); %delta_99(1.35,3582.7,u_inf/nu_inf);   % For ZPG we had: delta_99(1.35,2467,u_inf/nu_inf);

% Now Read and Scale Data
% Read (no filtering needed here, because no vane reflections)
ED_bf             = BLT_experiment_data_reader(fname_root_bf, d_list_bf, h);
% Set Crossflow scaling variables
ED_bf.set_crossflow_scaling_variables(S, delta_ref);


% Now show results on a first diagnostic plot
y_over_delta_first_range = 2; % This goas until above the data! (which is fine for diagnostic purposes!)
target_y = linspace(0,y_over_delta_first_range)';
z_ref    = -1;
u_bf_0h  = ED_bf.interp_velocities(0 , target_y,z_ref) / u_inf;
u_bf_50h = ED_bf.interp_velocities(50, target_y,z_ref) / u_inf; 
figure(3)
plot(u_bf_0h , target_y, u_bf_50h, target_y); grid on; hold on;
xlabel('Streamwise Speed'); ylabel('y/\delta - normal coordinate');

% Now make second computation of boundary layer parameters over selected
% range (based on observation of previous plot!
y_over_delta_ref_range = 1.6;

% Make a first computation at select places!
[dstr_m10, theta_m10, Hk_m10, target_y_m10, u_target_y_m10] = ED_bf.integral_parameters(0, -1.0, y_over_delta_ref_range , u_inf);
[dstr_m05, theta_m05, Hk_m05, target_y_m05, u_target_y_m05] = ED_bf.integral_parameters(0, -0.5, y_over_delta_ref_range , u_inf);
[dstr_c00, theta_c00, Hk_c00, target_y_c00, u_target_y_c00] = ED_bf.integral_parameters(0,  0.0, y_over_delta_ref_range , u_inf);
[dstr_M05, theta_M05, Hk_M05, target_y_M05, u_target_y_M05] = ED_bf.integral_parameters(0,  0.5, y_over_delta_ref_range , u_inf);
[dstr_M10, theta_M10, Hk_M10, target_y_M10, u_target_y_M10] = ED_bf.integral_parameters(0,  1.0, y_over_delta_ref_range , u_inf);
figure(4)
plot([-1, -0.5, 0, 0.5, 1], [Hk_m10, Hk_m05, Hk_c00, Hk_M05, Hk_M10]); grid on
xlabel('z/S'); ylabel('H_k');
% And a second computation everywhere!
target_x = 0; 
target_z = linspace(-1,1);          % List of z/S stances
Hk_z     = zeros(size(target_z ));  % Shape factor at z/S stance
theta_z  = zeros(size(target_z ));  % Momentum thicknes at z/S stance (over delta_ref)
for n_z = 1:length(target_z)
    % Reduce integration range a bit, to handle a bit of bleed (but not too much, trim 0.1*delta99 reference at top!)
    [~, theta_z(n_z), Hk_z(n_z)] = ED_bf.integral_parameters(0,  target_z(n_z), y_over_delta_ref_range , u_inf);
end
rt_z     = theta_z*u_inf/nu_inf;    % Reynolds theta at z/S stance

% Take reference values as averages
rt_0    = mean(rt_z);
Hk_0    = mean(Hk_z);

figure(1)
subplot(321)
plot(target_z , rt_z); grid on; hold on;
ph1 = plot([-1 1], rt_0*[1 1], 'k:'); 
axis([-1.0000    1.0000  2200 3200])
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

set(gcf , 'PaperType', 'A5');
orient landscape;
print ('-dpdf', ['initial_conditions_APG_subplot.pdf']) %#ok<NBRAK>


%% Now look at velocity profiles
[dstr_l, theta_l, Hk_l, target_y_l, u_target_y_l] = ED_bf.integral_parameters(0, -1.0, y_over_delta_ref_range, u_inf);
[dstr_c, theta_c, Hk_c, target_y_c, u_target_y_c] = ED_bf.integral_parameters(0,  0.0, y_over_delta_ref_range, u_inf);
[dstr_r, theta_r, Hk_r, target_y_r, u_target_y_r] = ED_bf.integral_parameters(0,  1.0, y_over_delta_ref_range, u_inf);

% Make a Swafford profile
SP = swafford_profile();
SP.update_hk_rt_pair(Hk_0, rt_0);
u_target_y_c_swafford = SP.evaluate_profile(target_y_c * delta_ref/ theta_c);

% Plot Comparison
figure(2)
subplot(241)
exp_l = plot(u_target_y_l,          target_y_l, 'Color', [0 0.4470 0.7410]); hold on; grid on;
exp_c = plot(u_target_y_c,          target_y_c, 'Color', [0 0.4470 0.7410]); 
exp_r = plot(u_target_y_r,          target_y_r, 'Color', [0 0.4470 0.7410]); 
swaf  = plot(u_target_y_c_swafford, target_y_c, 'k'); 
axis([0 1.05 min(target_y_l) max(target_y_l)]);
xlabel('u/u_inf'); ylabel('y/\delta_{ref}');

subplot(242)
exp_l = semilogy(u_target_y_l,          target_y_l, 'Color', [0 0.4470 0.7410]); hold on; grid on;
exp_c = semilogy(u_target_y_c,          target_y_c, 'Color', [0 0.4470 0.7410]); 
exp_r = semilogy(u_target_y_r,          target_y_r, 'Color', [0 0.4470 0.7410]); 
swaf  = semilogy(u_target_y_c_swafford, target_y_c, 'k'); 
axis([0 1.05 min(target_y_l) max(target_y_l)]);
xlabel('u/u_inf'); ylabel('log(y/\delta_{ref})');
legend([exp_l , swaf],  'Experimental Baseflow', ['Swafford (H_k =' num2str(Hk_0,3) '  Re_{\theta} =' num2str(rt_0,4) ])

set(gcf , 'PaperType', 'A5');
orient landscape;
print ('-dpdf', ['swafford_velocity_profile_APG.pdf']) %#ok<NBRAK>






