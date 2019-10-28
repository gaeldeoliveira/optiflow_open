% Script to study the accumulation of a passive scalar from the convection
% of a forcing field
%
% This is the equation for the transport of a passive scalar:
%   dcdt + (u*ddx + v*ddy)c = D*lap(c)
% Decompose c into c = u_bar + u_tilde, so that:
%   ddt(u_bar+u_tilde) + (u*ddx + v*ddy)(u_bar+u_tilde) = D*lap(u_bar+u_tilde)
% Where u_tilde = 0 for all x,y at t=0 and u_bar(x,y,t) is prescribed in
% space and time. Reworking the equation, leads to:
%   ddt(u_tilde)= - ddt(u_bar) - (u*ddx + v*ddy)(u_bar+u_tilde)
% Which is a field (x,y )of ODEs, and can therefore be integrated in time.
% This is what we will do  with a simple scheme in this  script.
% The equation can be expanded with the distributive property
%   ddt(u_tilde)= - ddt(u_bar) - (u*ddx + v*ddy)(u_bar+u_tilde) + D*lap(u_bar+u_tilde)
%   ddt(u_tilde)= - ddt(u_bar) - (u*ddx + v*ddy)(u_bar) - (u*ddx + v*ddy)(u_tilde) + D*lap(u_bar) + D*lap(u_tilde)
%
% Here are our variables:
%   The convecting field is provided by the induction of a vortex system:
%   u       is the x component of the convecting field
%   v       is the y component of the convecting field
%   The forcing field u_bar is constant in time, and has a simple form for
%   development purposes: 
%   u_bar   = @(x, y, t) alpha * y
%   Where alpha is a prescribed constant (slope). In the future u_bar will
%   be provided by the Swafford profile, and its parameters will be allowed
%   to vary in time.
%   Furthermore, in later versions, dt will be related to u_bar and allowed
%   to vary in space, but for now, let us simply fix it.
%
%   Same as alfa3, except that now the time step is not constant anyu
%
%   set(findall(gcf,'-property','FontSize'),'FontSize',18)
%
%TEMPORARY
%Regular      : state.y_v     = h_VG;             % [m /s ] distance of filament core center to wall
%Current trial: state.y_v     = 0.6 * h_VG;             % [m /s ] distance of filament core center to wall


%% Preparation
clear all; close all; clc; 

% Add Paths
addpath closure_relations/
addpath definitions/
addpath intersections/
addpath datagen/


%% Set inputs for BLT ZPG straight case (Baldacchinno)
vg_case.rho_inf         = 1.225;                  % [kg/m3]   Free Stream Density
vg_case.nu_inf          = 1.461e-5;               % [m2 s-1]  Kinematic Viscosity (nu = mu / rho)
vg_case.u_inf           = 15.16;                  % [m/s]     Unperturbed Free Stream Speed (U_inf!=Ue)

vg_case.hk0             = 1.41;                   % [adim.] Shape Factor
vg_case.rt0             = 2498;                   % [adim.] Reynolds Theta of Starting Flow

vg_case.D_VG            =  30.0e-3;                 % [m    ] distance between center of vortex pairs
vg_case.d_VG            =  12.5e-3;                 % [m    ] distance between VGs in a single pair
vg_case.c_VG            =  12.5e-3;                 % [m    ] chord of VG  (5cm)
vg_case.h_VG            =   5.0e-3;                 % [m    ] height of VG
vg_case.AoA_VG          = -18     ;                 % [deg  ] VG angle of attack (geometric, to free-stream)

vg_case.x0              =    0     * vg_case.h_VG;  % Starting Time (0, or -h_vg/4 , or /3/4*h_vg to enable predevelopment, vortex does not start at TE of Vg) % Exploration = - (3/4) * h_VG;
vg_case.x_end           =    50    * vg_case.h_VG;  % End Time (0.25 = 50 * h_vg)
vg_case.dx              =    0.0000125/2;           % TimeStep Weird stuff happens for dt = 0.1 (time step constraint?)
vg_case.dx_info         =    0.25 * vg_case.h_VG;   % Time between information events
vg_case.friendly_name   =    'BLT ZPG Straight Inflow, Rectangular Vane, Baldacchino 15';

% % Run case
snapshots = run_VG_case_alfa41(vg_case);
% % Plot results
n_snapshot = 1 * floor(length(snapshots.x_data) / 4); subplot(221); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 2 * floor(length(snapshots.x_data) / 4); subplot(222); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 3 * floor(length(snapshots.x_data) / 4); subplot(223); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 4 * floor(length(snapshots.x_data) / 4); subplot(224); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
% % Save data
save('datagen/snapshotmats/B001.mat', 'snapshots')

%% Set inputs for Wendt 2001 (Table 5 1st case)
vg_case.rho_inf         = 1.2256;                   % [kg/m3]   Free Stream Density
vg_case.nu_inf          = 1.48e-5;                  % [m2 s-1]  Kinematic Viscosity (nu = mu / rho)
vg_case.u_inf           = 85.00;                    % [m/s]     Unperturbed Free Stream Speed (U_inf!=Ue)

vg_case.hk0             = 1.41;                     % [adim.] Shape Factor (guesstimate based on fact that IFMF has ZPG section in region of VG application!)
vg_case.rt0             = 2500;                     % [adim.] Reynolds Theta of Starting Flow (temporary, initial guess for solver that will determine Rt from delta and hk guess)

vg_case.D_VG            =  0.053407075;             % [m    ] distance between center of vortex pairs
vg_case.d_VG            =  0.005363332;             % [m    ] distance between VGs in a single pair
vg_case.c_VG            =  0.0136;                 % [m    ] chord of VG  (5cm)
vg_case.h_VG            =  0.0051;                 % [m    ] height of VG
vg_case.AoA_VG          =  16    ;                 % [deg  ] VG angle of attack (geometric, to free-stream)

vg_case.x0              =    0     * vg_case.h_VG;  % Starting Time (0, or -h_vg/4 , or /3/4*h_vg to enable predevelopment, vortex does not start at TE of Vg) % Exploration = - (3/4) * h_VG;
vg_case.x_end           =    50    * vg_case.h_VG;  % End Time (0.25 = 50 * h_vg)
vg_case.dx              =    0.0000125/2;           % TimeStep Weird stuff happens for dt = 0.1 (time step constraint?)
vg_case.dx_info         =    0.25 * vg_case.h_VG;   % Time between information events
vg_case.friendly_name   =    'Wendt 2001, Table 5, case 1 per excel file';

% Now determine rt0 from delta99 (as provided in Wendt)
% BL is flat (ZPG) for IFMF cases so assume Hk=1.4
delta_expdata    = 0.01764314;
delta_99_of_rt = @(rt) delta_99(vg_case.hk0,rt,vg_case.u_inf/vg_case.nu_inf);
vg_case.rt0    = fzero(@(rt) delta_99_of_rt(rt)-delta_expdata, vg_case.rt0);
% % Run case
snapshots = run_VG_case_alfa41(vg_case);
% % Plot results
n_snapshot = 1 * floor(length(snapshots.x_data) / 4); subplot(221); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 2 * floor(length(snapshots.x_data) / 4); subplot(222); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 3 * floor(length(snapshots.x_data) / 4); subplot(223); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 4 * floor(length(snapshots.x_data) / 4); subplot(224); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
% % Save data
save('datagen/snapshotmats/W001.mat', 'snapshots')


%% Set inputs for Wendt 2001 (Table 5 2nd case)
vg_case.rho_inf         = 1.2256;                   % [kg/m3]   Free Stream Density
vg_case.nu_inf          = 1.48e-5;                  % [m2 s-1]  Kinematic Viscosity (nu = mu / rho)
vg_case.u_inf           = 85.00;                    % [m/s]     Unperturbed Free Stream Speed (U_inf!=Ue)

vg_case.hk0             = 1.41;                     % [adim.] Shape Factor (guesstimate based on fact that IFMF has ZPG section in region of VG application!)
vg_case.rt0             = 2500;                     % [adim.] Reynolds Theta of Starting Flow (temporary, initial guess for solver that will determine Rt from delta and hk guess)

vg_case.D_VG            =  0.053407075;             % [m    ] distance between center of vortex pairs
vg_case.d_VG            =  0.005363332;             % [m    ] distance between VGs in a single pair
vg_case.c_VG            =  0.0136;                 % [m    ] chord of VG  (5cm)
vg_case.h_VG            =  0.0102;                 % [m    ] height of VG
vg_case.AoA_VG          =  16    ;                 % [deg  ] VG angle of attack (geometric, to free-stream)

vg_case.x0              =    0     * vg_case.h_VG;  % Starting Time (0, or -h_vg/4 , or /3/4*h_vg to enable predevelopment, vortex does not start at TE of Vg) % Exploration = - (3/4) * h_VG;
vg_case.x_end           =    5    * vg_case.h_VG;  % End Time (0.25 = 50 * h_vg)
vg_case.dx              =    0.0000125/2;           % TimeStep Weird stuff happens for dt = 0.1 (time step constraint?)
vg_case.dx_info         =    0.25 * vg_case.h_VG;   % Time between information events
vg_case.friendly_name   =    'Wendt 2001, Table 5, case 2 per excel file';

% Now determine rt0 from delta99 (as provided in Wendt)
% BL is flat (ZPG) for IFMF cases so assume Hk=1.4
delta_expdata    = 0.01764314;
delta_99_of_rt = @(rt) delta_99(vg_case.hk0,rt,vg_case.u_inf/vg_case.nu_inf);
vg_case.rt0    = fzero(@(rt) delta_99_of_rt(rt)-delta_expdata, vg_case.rt0);
% % Run case
snapshots = run_VG_case_alfa41(vg_case);
% Filter out unconverged results (when filaments get out ot domain)
snapshots.N_snapshots = size(snapshots.u_tilde_data, 3);
u_invalidity_threshold = 3 * vg_case.u_inf;
n_snapshot = 1;
n_snapshot_max = snapshots.N_snapshots;
while n_snapshot < snapshots.N_snapshots
    if max(max(abs(snapshots.u_tilde_data(:,:,n_snapshot)))) > u_invalidity_threshold
        n_snapshot_max = n_snapshot-1;
    end
    n_snapshot = n_snapshot + 1;
end
snapshots_filtered = snapshots;
snapshots_filtered.u_tilde_data = snapshots_filtered.u_tilde_data(:,:,1:n_snapshot_max);
snapshots_filtered.x_data = snapshots_filtered.x_data(1:n_snapshot_max);
snapshots_filtered.hk_data = snapshots_filtered.hk_data(1:n_snapshot_max);
snapshots_filtered.theta_data = snapshots_filtered.theta_data(1:n_snapshot_max);
snapshots_filtered.ctau_data = snapshots_filtered.ctau_data(1:n_snapshot_max);
snapshots_filtered.ue_data = snapshots_filtered.ue_data(1:n_snapshot_max);
snapshots_filtered.rt_data = snapshots_filtered.rt_data(1:n_snapshot_max);
snapshots_filtered.y_v_over_d0_data = snapshots_filtered.y_v_over_d0_data(1:n_snapshot_max);
snapshots_filtered.z_v_over_S_data = snapshots_filtered.z_v_over_S_data(1:n_snapshot_max);
snapshots_filtered.gamma_v_data = snapshots_filtered.gamma_v_data(1:n_snapshot_max);
snapshots_filtered.sigma_v_data = snapshots_filtered.sigma_v_data(1:n_snapshot_max);
snapshots_filtered.t_v_data = snapshots_filtered.t_v_data(1:n_snapshot_max);
snapshots_filtered.N_snapshots = size(snapshots_filtered.u_tilde_data, 3);
% Now replace original structure
snapshots = snapshots_filtered;
% % Plot results
n_snapshot = 1 * floor(length(snapshots.x_data) / 4); subplot(221); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 2 * floor(length(snapshots.x_data) / 4); subplot(222); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 3 * floor(length(snapshots.x_data) / 4); subplot(223); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 4 * floor(length(snapshots.x_data) / 4); subplot(224); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
% % Save data
save('datagen/snapshotmats/W002.mat', 'snapshots')

%% Set inputs for Wendt 2001 (Table 5 3rd case)
vg_case.rho_inf         = 1.2256;                   % [kg/m3]   Free Stream Density
vg_case.nu_inf          = 1.48e-5;                  % [m2 s-1]  Kinematic Viscosity (nu = mu / rho)
vg_case.u_inf           = 85.00;                    % [m/s]     Unperturbed Free Stream Speed (U_inf!=Ue)

vg_case.hk0             = 1.41;                     % [adim.] Shape Factor (guesstimate based on fact that IFMF has ZPG section in region of VG application!)
vg_case.rt0             = 2500;                     % [adim.] Reynolds Theta of Starting Flow (temporary, initial guess for solver that will determine Rt from delta and hk guess)

vg_case.D_VG            =  0.10681415;              % [m    ] distance between center of vortex pairs
vg_case.d_VG            =  0.016011123;             % [m    ] distance between VGs in a single pair
vg_case.c_VG            =  0.0406;                  % [m    ] chord of VG  (5cm)
vg_case.h_VG            =  0.0051;                  % [m    ] height of VG
vg_case.AoA_VG          =  16    ;                  % [deg  ] VG angle of attack (geometric, to free-stream)

vg_case.x0              =    0     * vg_case.h_VG;  % Starting Time (0, or -h_vg/4 , or /3/4*h_vg to enable predevelopment, vortex does not start at TE of Vg) % Exploration = - (3/4) * h_VG;
vg_case.x_end           =    50    * vg_case.h_VG;  % End Time (0.25 = 50 * h_vg)
vg_case.dx              =    0.0000125/2;           % TimeStep Weird stuff happens for dt = 0.1 (time step constraint?)
vg_case.dx_info         =    0.25 * vg_case.h_VG;   % Time between information events
vg_case.friendly_name   =    'Wendt 2001, Table 5, case 3 per excel file';

% Now determine rt0 from delta99 (as provided in Wendt)
% BL is flat (ZPG) for IFMF cases so assume Hk=1.4
delta_expdata    = 0.01764314;
delta_99_of_rt = @(rt) delta_99(vg_case.hk0,rt,vg_case.u_inf/vg_case.nu_inf);
vg_case.rt0    = fzero(@(rt) delta_99_of_rt(rt)-delta_expdata, vg_case.rt0);
% % Run case
snapshots = run_VG_case_alfa41(vg_case);
% Filter out unconverged results (when filaments get out ot domain)
snapshots.N_snapshots = size(snapshots.u_tilde_data, 3);
u_invalidity_threshold = 3 * vg_case.u_inf;
n_snapshot = 1;
n_snapshot_max = snapshots.N_snapshots;
while n_snapshot < snapshots.N_snapshots
    if max(max(abs(snapshots.u_tilde_data(:,:,n_snapshot)))) > u_invalidity_threshold
        n_snapshot_max = n_snapshot-1;
    end
    n_snapshot = n_snapshot + 1;
end
snapshots_filtered = snapshots;
snapshots_filtered.u_tilde_data = snapshots_filtered.u_tilde_data(:,:,1:n_snapshot_max);
snapshots_filtered.x_data = snapshots_filtered.x_data(1:n_snapshot_max);
snapshots_filtered.hk_data = snapshots_filtered.hk_data(1:n_snapshot_max);
snapshots_filtered.theta_data = snapshots_filtered.theta_data(1:n_snapshot_max);
snapshots_filtered.ctau_data = snapshots_filtered.ctau_data(1:n_snapshot_max);
snapshots_filtered.ue_data = snapshots_filtered.ue_data(1:n_snapshot_max);
snapshots_filtered.rt_data = snapshots_filtered.rt_data(1:n_snapshot_max);
snapshots_filtered.y_v_over_d0_data = snapshots_filtered.y_v_over_d0_data(1:n_snapshot_max);
snapshots_filtered.z_v_over_S_data = snapshots_filtered.z_v_over_S_data(1:n_snapshot_max);
snapshots_filtered.gamma_v_data = snapshots_filtered.gamma_v_data(1:n_snapshot_max);
snapshots_filtered.sigma_v_data = snapshots_filtered.sigma_v_data(1:n_snapshot_max);
snapshots_filtered.t_v_data = snapshots_filtered.t_v_data(1:n_snapshot_max);
snapshots_filtered.N_snapshots = size(snapshots_filtered.u_tilde_data, 3);
% Now replace original structure
snapshots = snapshots_filtered;
% % Plot results
n_snapshot = 1 * floor(length(snapshots.x_data) / 4); subplot(221); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 2 * floor(length(snapshots.x_data) / 4); subplot(222); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 3 * floor(length(snapshots.x_data) / 4); subplot(223); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 4 * floor(length(snapshots.x_data) / 4); subplot(224); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
% % Save data
save('datagen/snapshotmats/W003.mat', 'snapshots')


%% Set inputs for Wendt 2001 (Table 5 4th case)
vg_case.rho_inf         = 1.2256;                   % [kg/m3]   Free Stream Density
vg_case.nu_inf          = 1.48e-5;                  % [m2 s-1]  Kinematic Viscosity (nu = mu / rho)
vg_case.u_inf           = 85.00;                    % [m/s]     Unperturbed Free Stream Speed (U_inf!=Ue)

vg_case.hk0             = 1.41;                     % [adim.] Shape Factor (guesstimate based on fact that IFMF has ZPG section in region of VG application!)
vg_case.rt0             = 2500;                     % [adim.] Reynolds Theta of Starting Flow (temporary, initial guess for solver that will determine Rt from delta and hk guess)

vg_case.D_VG            =  0.10681415;              % [m    ] distance between center of vortex pairs
vg_case.d_VG            =  0.043213123;             % [m    ] distance between VGs in a single pair
vg_case.c_VG            =  0.0406;                  % [m    ] chord of VG  (5cm)
vg_case.h_VG            =  0.0051;                  % [m    ] height of VG
vg_case.AoA_VG          =  16    ;                  % [deg  ] VG angle of attack (geometric, to free-stream)

vg_case.x0              =    0     * vg_case.h_VG;  % Starting Time (0, or -h_vg/4 , or /3/4*h_vg to enable predevelopment, vortex does not start at TE of Vg) % Exploration = - (3/4) * h_VG;
vg_case.x_end           =    50    * vg_case.h_VG;  % End Time (0.25 = 50 * h_vg)
vg_case.dx              =    0.0000125/2;           % TimeStep Weird stuff happens for dt = 0.1 (time step constraint?)
vg_case.dx_info         =    0.25  * vg_case.h_VG;   % Time between information events
vg_case.friendly_name   =    'Wendt 2001, Table 5, case 4 per excel file';

% Now determine rt0 from delta99 (as provided in Wendt)
% BL is flat (ZPG) for IFMF cases so assume Hk=1.4
delta_expdata    = 0.01764314;
delta_99_of_rt = @(rt) delta_99(vg_case.hk0,rt,vg_case.u_inf/vg_case.nu_inf);
vg_case.rt0    = fzero(@(rt) delta_99_of_rt(rt)-delta_expdata, vg_case.rt0);
% % Run case
snapshots = run_VG_case_alfa41(vg_case);
% Filter out unconverged results (when filaments get out ot domain)
snapshots.N_snapshots = size(snapshots.u_tilde_data, 3);
u_invalidity_threshold = 3 * vg_case.u_inf;
n_snapshot = 1;
n_snapshot_max = snapshots.N_snapshots;
while n_snapshot < snapshots.N_snapshots
    if max(max(abs(snapshots.u_tilde_data(:,:,n_snapshot)))) > u_invalidity_threshold
        n_snapshot_max = n_snapshot-1;
    end
    n_snapshot = n_snapshot + 1;
end
snapshots_filtered = snapshots;
snapshots_filtered.u_tilde_data = snapshots_filtered.u_tilde_data(:,:,1:n_snapshot_max);
snapshots_filtered.x_data = snapshots_filtered.x_data(1:n_snapshot_max);
snapshots_filtered.hk_data = snapshots_filtered.hk_data(1:n_snapshot_max);
snapshots_filtered.theta_data = snapshots_filtered.theta_data(1:n_snapshot_max);
snapshots_filtered.ctau_data = snapshots_filtered.ctau_data(1:n_snapshot_max);
snapshots_filtered.ue_data = snapshots_filtered.ue_data(1:n_snapshot_max);
snapshots_filtered.rt_data = snapshots_filtered.rt_data(1:n_snapshot_max);
snapshots_filtered.y_v_over_d0_data = snapshots_filtered.y_v_over_d0_data(1:n_snapshot_max);
snapshots_filtered.z_v_over_S_data = snapshots_filtered.z_v_over_S_data(1:n_snapshot_max);
snapshots_filtered.gamma_v_data = snapshots_filtered.gamma_v_data(1:n_snapshot_max);
snapshots_filtered.sigma_v_data = snapshots_filtered.sigma_v_data(1:n_snapshot_max);
snapshots_filtered.t_v_data = snapshots_filtered.t_v_data(1:n_snapshot_max);
snapshots_filtered.N_snapshots = size(snapshots_filtered.u_tilde_data, 3);
% Now replace original structure
snapshots = snapshots_filtered;
% % Plot results
n_snapshot = 1 * floor(length(snapshots.x_data) / 4); subplot(221); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 2 * floor(length(snapshots.x_data) / 4); subplot(222); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 3 * floor(length(snapshots.x_data) / 4); subplot(223); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 4 * floor(length(snapshots.x_data) / 4); subplot(224); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
% % Save data
save('datagen/snapshotmats/W004.mat', 'snapshots')


%% Set inputs for Wendt 2001 (Table 6 5th case)
vg_case.rho_inf         = 1.2256;                   % [kg/m3]   Free Stream Density
vg_case.nu_inf          = 1.48e-5;                  % [m2 s-1]  Kinematic Viscosity (nu = mu / rho)
vg_case.u_inf           = 85.00;                    % [m/s]     Unperturbed Free Stream Speed (U_inf!=Ue)

vg_case.hk0             = 1.41;                     % [adim.] Shape Factor (guesstimate based on fact that IFMF has ZPG section in region of VG application!)
vg_case.rt0             = 2500;                     % [adim.] Reynolds Theta of Starting Flow (temporary, initial guess for solver that will determine Rt from delta and hk guess)

vg_case.D_VG            =  0.053407075;              % [m    ] distance between center of vortex pairs
vg_case.d_VG            =  0.012860668;             % [m    ] distance between VGs in a single pair
vg_case.c_VG            =  0.0136;                  % [m    ] chord of VG  (5cm)
vg_case.h_VG            =  0.0051;                  % [m    ] height of VG
vg_case.AoA_VG          =  -16;                  % [deg  ] VG angle of attack (geometric, to free-stream)

vg_case.x0              =    0     * vg_case.h_VG;  % Starting Time (0, or -h_vg/4 , or /3/4*h_vg to enable predevelopment, vortex does not start at TE of Vg) % Exploration = - (3/4) * h_VG;
vg_case.x_end           =    50    * vg_case.h_VG;  % End Time (0.25 = 50 * h_vg)
vg_case.dx              =    0.0000125/2;           % TimeStep Weird stuff happens for dt = 0.1 (time step constraint?)
vg_case.dx_info         =    0.25  * vg_case.h_VG;   % Time between information events
vg_case.friendly_name   =    'Wendt 2001, Table 6, case 5 per excel file';

% Now determine rt0 from delta99 (as provided in Wendt)
% BL is flat (ZPG) for IFMF cases so assume Hk=1.4
delta_expdata    = 0.01764314;
delta_99_of_rt = @(rt) delta_99(vg_case.hk0,rt,vg_case.u_inf/vg_case.nu_inf);
vg_case.rt0    = fzero(@(rt) delta_99_of_rt(rt)-delta_expdata, vg_case.rt0);
% % Run case
snapshots = run_VG_case_alfa41(vg_case);
% Filter out unconverged results (when filaments get out ot domain)
snapshots.N_snapshots = size(snapshots.u_tilde_data, 3);
u_invalidity_threshold = 3 * vg_case.u_inf;
n_snapshot = 1;
n_snapshot_max = snapshots.N_snapshots;
while n_snapshot < snapshots.N_snapshots
    if max(max(abs(snapshots.u_tilde_data(:,:,n_snapshot)))) > u_invalidity_threshold
        n_snapshot_max = n_snapshot-1;
    end
    n_snapshot = n_snapshot + 1;
end
snapshots_filtered = snapshots;
snapshots_filtered.u_tilde_data = snapshots_filtered.u_tilde_data(:,:,1:n_snapshot_max);
snapshots_filtered.x_data = snapshots_filtered.x_data(1:n_snapshot_max);
snapshots_filtered.hk_data = snapshots_filtered.hk_data(1:n_snapshot_max);
snapshots_filtered.theta_data = snapshots_filtered.theta_data(1:n_snapshot_max);
snapshots_filtered.ctau_data = snapshots_filtered.ctau_data(1:n_snapshot_max);
snapshots_filtered.ue_data = snapshots_filtered.ue_data(1:n_snapshot_max);
snapshots_filtered.rt_data = snapshots_filtered.rt_data(1:n_snapshot_max);
snapshots_filtered.y_v_over_d0_data = snapshots_filtered.y_v_over_d0_data(1:n_snapshot_max);
snapshots_filtered.z_v_over_S_data = snapshots_filtered.z_v_over_S_data(1:n_snapshot_max);
snapshots_filtered.gamma_v_data = snapshots_filtered.gamma_v_data(1:n_snapshot_max);
snapshots_filtered.sigma_v_data = snapshots_filtered.sigma_v_data(1:n_snapshot_max);
snapshots_filtered.t_v_data = snapshots_filtered.t_v_data(1:n_snapshot_max);
snapshots_filtered.N_snapshots = size(snapshots_filtered.u_tilde_data, 3);
% Now replace original structure
snapshots = snapshots_filtered;
% % Plot results
n_snapshot = 1 * floor(length(snapshots.x_data) / 4); subplot(221); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 2 * floor(length(snapshots.x_data) / 4); subplot(222); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 3 * floor(length(snapshots.x_data) / 4); subplot(223); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 4 * floor(length(snapshots.x_data) / 4); subplot(224); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
% % Save data
save('datagen/snapshotmats/W005.mat', 'snapshots')


%% Set inputs for Wendt 2001 (Table 6 6th case)
vg_case.rho_inf         = 1.2256;                   % [kg/m3]   Free Stream Density
vg_case.nu_inf          = 1.48e-5;                  % [m2 s-1]  Kinematic Viscosity (nu = mu / rho)
vg_case.u_inf           = 85.00;                    % [m/s]     Unperturbed Free Stream Speed (U_inf!=Ue)

vg_case.hk0             = 1.41;                     % [adim.] Shape Factor (guesstimate based on fact that IFMF has ZPG section in region of VG application!)
vg_case.rt0             = 2500;                     % [adim.] Reynolds Theta of Starting Flow (temporary, initial guess for solver that will determine Rt from delta and hk guess)

vg_case.D_VG            =  0.053407075;              % [m    ] distance between center of vortex pairs
vg_case.d_VG            =  0.012860668;             % [m    ] distance between VGs in a single pair
vg_case.c_VG            =  0.0136;                  % [m    ] chord of VG  (5cm)
vg_case.h_VG            =  0.0102;                  % [m    ] height of VG
vg_case.AoA_VG          =  -16;                  % [deg  ] VG angle of attack (geometric, to free-stream)

vg_case.x0              =    0     * vg_case.h_VG;  % Starting Time (0, or -h_vg/4 , or /3/4*h_vg to enable predevelopment, vortex does not start at TE of Vg) % Exploration = - (3/4) * h_VG;
vg_case.x_end           =    50    * vg_case.h_VG;  % End Time (0.25 = 50 * h_vg)
vg_case.dx              =    0.0000125/2;           % TimeStep Weird stuff happens for dt = 0.1 (time step constraint?)
vg_case.dx_info         =    0.25  * vg_case.h_VG;   % Time between information events
vg_case.friendly_name   =    'Wendt 2001, Table 6, case 6 per excel file';

% Now determine rt0 from delta99 (as provided in Wendt)
% BL is flat (ZPG) for IFMF cases so assume Hk=1.4
delta_expdata    = 0.01764314;
delta_99_of_rt = @(rt) delta_99(vg_case.hk0,rt,vg_case.u_inf/vg_case.nu_inf);
vg_case.rt0    = fzero(@(rt) delta_99_of_rt(rt)-delta_expdata, vg_case.rt0);
% % Run case
snapshots = run_VG_case_alfa41(vg_case);
% Filter out unconverged results (when filaments get out ot domain)
snapshots.N_snapshots = size(snapshots.u_tilde_data, 3);
u_invalidity_threshold = 3 * vg_case.u_inf;
n_snapshot = 1;
n_snapshot_max = snapshots.N_snapshots;
while n_snapshot < snapshots.N_snapshots
    if max(max(abs(snapshots.u_tilde_data(:,:,n_snapshot)))) > u_invalidity_threshold
        n_snapshot_max = n_snapshot-1;
    end
    n_snapshot = n_snapshot + 1;
end
snapshots_filtered = snapshots;
snapshots_filtered.u_tilde_data = snapshots_filtered.u_tilde_data(:,:,1:n_snapshot_max);
snapshots_filtered.x_data = snapshots_filtered.x_data(1:n_snapshot_max);
snapshots_filtered.hk_data = snapshots_filtered.hk_data(1:n_snapshot_max);
snapshots_filtered.theta_data = snapshots_filtered.theta_data(1:n_snapshot_max);
snapshots_filtered.ctau_data = snapshots_filtered.ctau_data(1:n_snapshot_max);
snapshots_filtered.ue_data = snapshots_filtered.ue_data(1:n_snapshot_max);
snapshots_filtered.rt_data = snapshots_filtered.rt_data(1:n_snapshot_max);
snapshots_filtered.y_v_over_d0_data = snapshots_filtered.y_v_over_d0_data(1:n_snapshot_max);
snapshots_filtered.z_v_over_S_data = snapshots_filtered.z_v_over_S_data(1:n_snapshot_max);
snapshots_filtered.gamma_v_data = snapshots_filtered.gamma_v_data(1:n_snapshot_max);
snapshots_filtered.sigma_v_data = snapshots_filtered.sigma_v_data(1:n_snapshot_max);
snapshots_filtered.t_v_data = snapshots_filtered.t_v_data(1:n_snapshot_max);
snapshots_filtered.N_snapshots = size(snapshots_filtered.u_tilde_data, 3);
% Now replace original structure
snapshots = snapshots_filtered;
% % Plot results
n_snapshot = 1 * floor(length(snapshots.x_data) / 4); subplot(221); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 2 * floor(length(snapshots.x_data) / 4); subplot(222); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 3 * floor(length(snapshots.x_data) / 4); subplot(223); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 4 * floor(length(snapshots.x_data) / 4); subplot(224); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
% % Save data
save('datagen/snapshotmats/W006.mat', 'snapshots')

%% Set inputs for Wendt 2001 (Table 6 7th case)
vg_case.rho_inf         = 1.2256;                   % [kg/m3]   Free Stream Density
vg_case.nu_inf          = 1.48e-5;                  % [m2 s-1]  Kinematic Viscosity (nu = mu / rho)
vg_case.u_inf           = 85.00;                    % [m/s]     Unperturbed Free Stream Speed (U_inf!=Ue)

vg_case.hk0             = 1.41;                     % [adim.] Shape Factor (guesstimate based on fact that IFMF has ZPG section in region of VG application!)
vg_case.rt0             = 2500;                     % [adim.] Reynolds Theta of Starting Flow (temporary, initial guess for solver that will determine Rt from delta and hk guess)

vg_case.D_VG            =  0.10681415;              % [m    ] distance between center of vortex pairs
vg_case.d_VG            =  0.038392877;             % [m    ] distance between VGs in a single pair
vg_case.c_VG            =  0.0406;                  % [m    ] chord of VG  (5cm)
vg_case.h_VG            =  0.0051;                  % [m    ] height of VG
vg_case.AoA_VG          =  -16;                  % [deg  ] VG angle of attack (geometric, to free-stream)

vg_case.x0              =    0     * vg_case.h_VG;  % Starting Time (0, or -h_vg/4 , or /3/4*h_vg to enable predevelopment, vortex does not start at TE of Vg) % Exploration = - (3/4) * h_VG;
vg_case.x_end           =    50    * vg_case.h_VG;  % End Time (0.25 = 50 * h_vg)
vg_case.dx              =    0.0000125/2;           % TimeStep Weird stuff happens for dt = 0.1 (time step constraint?)
vg_case.dx_info         =    0.25  * vg_case.h_VG;   % Time between information events
vg_case.friendly_name   =    'Wendt 2001, Table 6, case 7 per excel file';

% Now determine rt0 from delta99 (as provided in Wendt)
% BL is flat (ZPG) for IFMF cases so assume Hk=1.4
delta_expdata    = 0.01764314;
delta_99_of_rt = @(rt) delta_99(vg_case.hk0,rt,vg_case.u_inf/vg_case.nu_inf);
vg_case.rt0    = fzero(@(rt) delta_99_of_rt(rt)-delta_expdata, vg_case.rt0);
% % Run case
snapshots = run_VG_case_alfa41(vg_case);
% Filter out unconverged results (when filaments get out ot domain)
snapshots.N_snapshots = size(snapshots.u_tilde_data, 3);
u_invalidity_threshold = 3 * vg_case.u_inf;
n_snapshot = 1;
n_snapshot_max = snapshots.N_snapshots;
while n_snapshot < snapshots.N_snapshots
    if max(max(abs(snapshots.u_tilde_data(:,:,n_snapshot)))) > u_invalidity_threshold
        n_snapshot_max = n_snapshot-1;
    end
    n_snapshot = n_snapshot + 1;
end
snapshots_filtered = snapshots;
snapshots_filtered.u_tilde_data = snapshots_filtered.u_tilde_data(:,:,1:n_snapshot_max);
snapshots_filtered.x_data = snapshots_filtered.x_data(1:n_snapshot_max);
snapshots_filtered.hk_data = snapshots_filtered.hk_data(1:n_snapshot_max);
snapshots_filtered.theta_data = snapshots_filtered.theta_data(1:n_snapshot_max);
snapshots_filtered.ctau_data = snapshots_filtered.ctau_data(1:n_snapshot_max);
snapshots_filtered.ue_data = snapshots_filtered.ue_data(1:n_snapshot_max);
snapshots_filtered.rt_data = snapshots_filtered.rt_data(1:n_snapshot_max);
snapshots_filtered.y_v_over_d0_data = snapshots_filtered.y_v_over_d0_data(1:n_snapshot_max);
snapshots_filtered.z_v_over_S_data = snapshots_filtered.z_v_over_S_data(1:n_snapshot_max);
snapshots_filtered.gamma_v_data = snapshots_filtered.gamma_v_data(1:n_snapshot_max);
snapshots_filtered.sigma_v_data = snapshots_filtered.sigma_v_data(1:n_snapshot_max);
snapshots_filtered.t_v_data = snapshots_filtered.t_v_data(1:n_snapshot_max);
snapshots_filtered.N_snapshots = size(snapshots_filtered.u_tilde_data, 3);
% Now replace original structure
snapshots = snapshots_filtered;
% % Plot results
n_snapshot = 1 * floor(length(snapshots.x_data) / 4); subplot(221); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 2 * floor(length(snapshots.x_data) / 4); subplot(222); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 3 * floor(length(snapshots.x_data) / 4); subplot(223); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 4 * floor(length(snapshots.x_data) / 4); subplot(224); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
% % Save data
save('datagen/snapshotmats/W007.mat', 'snapshots')


%% Set inputs for Wendt 2001 (Table 6 8th case)
vg_case.rho_inf         = 1.2256;                   % [kg/m3]   Free Stream Density
vg_case.nu_inf          = 1.48e-5;                  % [m2 s-1]  Kinematic Viscosity (nu = mu / rho)
vg_case.u_inf           = 85.00;                    % [m/s]     Unperturbed Free Stream Speed (U_inf!=Ue)

vg_case.hk0             = 1.41;                     % [adim.] Shape Factor (guesstimate based on fact that IFMF has ZPG section in region of VG application!)
vg_case.rt0             = 2500;                     % [adim.] Reynolds Theta of Starting Flow (temporary, initial guess for solver that will determine Rt from delta and hk guess)

vg_case.D_VG            =  0.10681415;              % [m    ] distance between center of vortex pairs
vg_case.d_VG            =  0.065594877;             % [m    ] distance between VGs in a single pair
vg_case.c_VG            =  0.0406;                  % [m    ] chord of VG  (5cm)
vg_case.h_VG            =  0.0051;                  % [m    ] height of VG
vg_case.AoA_VG          =  -16;                     % [deg  ] VG angle of attack (geometric, to free-stream)

vg_case.x0              =    0     * vg_case.h_VG;  % Starting Time (0, or -h_vg/4 , or /3/4*h_vg to enable predevelopment, vortex does not start at TE of Vg) % Exploration = - (3/4) * h_VG;
vg_case.x_end           =    50    * vg_case.h_VG;  % End Time (0.25 = 50 * h_vg)
vg_case.dx              =    0.0000125/2;           % TimeStep Weird stuff happens for dt = 0.1 (time step constraint?)
vg_case.dx_info         =    0.25  * vg_case.h_VG;   % Time between information events
vg_case.friendly_name   =    'Wendt 2001, Table 6, case 9 per excel file';

% Now determine rt0 from delta99 (as provided in Wendt)
% BL is flat (ZPG) for IFMF cases so assume Hk=1.4
delta_expdata    = 0.01764314;
delta_99_of_rt = @(rt) delta_99(vg_case.hk0,rt,vg_case.u_inf/vg_case.nu_inf);
vg_case.rt0    = fzero(@(rt) delta_99_of_rt(rt)-delta_expdata, vg_case.rt0);
% % Run case
snapshots = run_VG_case_alfa41(vg_case);
% Filter out unconverged results (when filaments get out ot domain)
snapshots.N_snapshots = size(snapshots.u_tilde_data, 3);
u_invalidity_threshold = 3 * vg_case.u_inf;
n_snapshot = 1;
n_snapshot_max = snapshots.N_snapshots;
while n_snapshot < snapshots.N_snapshots
    if max(max(abs(snapshots.u_tilde_data(:,:,n_snapshot)))) > u_invalidity_threshold
        n_snapshot_max = n_snapshot-1;
    end
    n_snapshot = n_snapshot + 1;
end
snapshots_filtered = snapshots;
snapshots_filtered.u_tilde_data = snapshots_filtered.u_tilde_data(:,:,1:n_snapshot_max);
snapshots_filtered.x_data = snapshots_filtered.x_data(1:n_snapshot_max);
snapshots_filtered.hk_data = snapshots_filtered.hk_data(1:n_snapshot_max);
snapshots_filtered.theta_data = snapshots_filtered.theta_data(1:n_snapshot_max);
snapshots_filtered.ctau_data = snapshots_filtered.ctau_data(1:n_snapshot_max);
snapshots_filtered.ue_data = snapshots_filtered.ue_data(1:n_snapshot_max);
snapshots_filtered.rt_data = snapshots_filtered.rt_data(1:n_snapshot_max);
snapshots_filtered.y_v_over_d0_data = snapshots_filtered.y_v_over_d0_data(1:n_snapshot_max);
snapshots_filtered.z_v_over_S_data = snapshots_filtered.z_v_over_S_data(1:n_snapshot_max);
snapshots_filtered.gamma_v_data = snapshots_filtered.gamma_v_data(1:n_snapshot_max);
snapshots_filtered.sigma_v_data = snapshots_filtered.sigma_v_data(1:n_snapshot_max);
snapshots_filtered.t_v_data = snapshots_filtered.t_v_data(1:n_snapshot_max);
snapshots_filtered.N_snapshots = size(snapshots_filtered.u_tilde_data, 3);
% Now replace original structure
snapshots = snapshots_filtered;
% % Plot results
n_snapshot = 1 * floor(length(snapshots.x_data) / 4); subplot(221); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 2 * floor(length(snapshots.x_data) / 4); subplot(222); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 3 * floor(length(snapshots.x_data) / 4); subplot(223); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 4 * floor(length(snapshots.x_data) / 4); subplot(224); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
% % Save data
save('datagen/snapshotmats/W008.mat', 'snapshots')


%% Set inputs for Logdberg 2009 (JFM paper 1st case)
vg_case.rho_inf         = 1.2256;                   % [kg/m3]   Free Stream Density
vg_case.nu_inf          = 1.48e-5;                  % [m2 s-1]  Kinematic Viscosity (nu = mu / rho)
vg_case.u_inf           = 26.50;                    % [m/s]     Unperturbed Free Stream Speed (U_inf!=Ue)

vg_case.hk0             = 1.41;                     % [adim.] Shape Factor (guesstimate based on fact that IFMF has ZPG section in region of VG application!)
vg_case.rt0             = 2500;                     % [adim.] Reynolds Theta of Starting Flow (temporary, initial guess for solver that will determine Rt from delta and hk guess)

vg_case.D_VG            =  0.05;                    % [m    ] distance between center of vortex pairs
vg_case.d_VG            =  0.017323085;             % [m    ] distance between VGs in a single pair
vg_case.c_VG            =  0.018634971;             % [m    ] chord of VG  (5cm)
vg_case.h_VG            =  0.006;                   % [m    ] height of VG
vg_case.AoA_VG          =  -15;                     % [deg  ] VG angle of attack (geometric, to free-stream)

vg_case.x0              =    0     * vg_case.h_VG;  % Starting Time (0, or -h_vg/4 , or /3/4*h_vg to enable predevelopment, vortex does not start at TE of Vg) % Exploration = - (3/4) * h_VG;
vg_case.x_end           =    50    * vg_case.h_VG;  % End Time (0.25 = 50 * h_vg)
vg_case.dx              =    0.0000125/8;           % TimeStep Weird stuff happens for dt = 0.1 (time step constraint?)
vg_case.dx_info         =    0.25  * vg_case.h_VG;   % Time between information events
vg_case.friendly_name   =    'Logdberg 2009, JFM paper, 1st case';

% Now determine rt0 from delta99 (as provided in Wendt)
% BL is flat (ZPG) for IFMF cases so assume Hk=1.4
delta_expdata           = 0.027580938;
delta_99_of_rt          = @(rt) delta_99(vg_case.hk0,rt,vg_case.u_inf/vg_case.nu_inf);
vg_case.rt0             = fzero(@(rt) delta_99_of_rt(rt)-delta_expdata, vg_case.rt0);
% % Run case
snapshots = run_VG_case_alfa41(vg_case);
% Filter out unconverged results (when filaments get out ot domain)
snapshots.N_snapshots = size(snapshots.u_tilde_data, 3);
u_invalidity_threshold = 3 * vg_case.u_inf;
n_snapshot = 1;
n_snapshot_max = snapshots.N_snapshots;
while n_snapshot < snapshots.N_snapshots
    if max(max(abs(snapshots.u_tilde_data(:,:,n_snapshot)))) > u_invalidity_threshold
        n_snapshot_max = n_snapshot-1;
    end
    n_snapshot = n_snapshot + 1;
end
snapshots_filtered = snapshots;
snapshots_filtered.u_tilde_data = snapshots_filtered.u_tilde_data(:,:,1:n_snapshot_max);
snapshots_filtered.x_data = snapshots_filtered.x_data(1:n_snapshot_max);
snapshots_filtered.hk_data = snapshots_filtered.hk_data(1:n_snapshot_max);
snapshots_filtered.theta_data = snapshots_filtered.theta_data(1:n_snapshot_max);
snapshots_filtered.ctau_data = snapshots_filtered.ctau_data(1:n_snapshot_max);
snapshots_filtered.ue_data = snapshots_filtered.ue_data(1:n_snapshot_max);
snapshots_filtered.rt_data = snapshots_filtered.rt_data(1:n_snapshot_max);
snapshots_filtered.y_v_over_d0_data = snapshots_filtered.y_v_over_d0_data(1:n_snapshot_max);
snapshots_filtered.z_v_over_S_data = snapshots_filtered.z_v_over_S_data(1:n_snapshot_max);
snapshots_filtered.gamma_v_data = snapshots_filtered.gamma_v_data(1:n_snapshot_max);
snapshots_filtered.sigma_v_data = snapshots_filtered.sigma_v_data(1:n_snapshot_max);
snapshots_filtered.t_v_data = snapshots_filtered.t_v_data(1:n_snapshot_max);
snapshots_filtered.N_snapshots = size(snapshots_filtered.u_tilde_data, 3);
% Now replace original structure
snapshots = snapshots_filtered;
% % Plot results
n_snapshot = 1 * floor(length(snapshots.x_data) / 4); subplot(221); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 2 * floor(length(snapshots.x_data) / 4); subplot(222); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 3 * floor(length(snapshots.x_data) / 4); subplot(223); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 4 * floor(length(snapshots.x_data) / 4); subplot(224); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
% % Save data
save('datagen/snapshotmats/L001.mat', 'snapshots')


%% Set inputs for Logdberg 2009 (JFM paper 2nd case)
vg_case.rho_inf         = 1.2256;                   % [kg/m3]   Free Stream Density
vg_case.nu_inf          = 1.48e-5;                  % [m2 s-1]  Kinematic Viscosity (nu = mu / rho)
vg_case.u_inf           = 26.50;                    % [m/s]     Unperturbed Free Stream Speed (U_inf!=Ue)

vg_case.hk0             = 1.41;                     % [adim.] Shape Factor (guesstimate based on fact that IFMF has ZPG section in region of VG application!)
vg_case.rt0             = 2500;                     % [adim.] Reynolds Theta of Starting Flow (temporary, initial guess for solver that will determine Rt from delta and hk guess)

vg_case.D_VG            =  0.083;                   % [m    ] distance between center of vortex pairs
vg_case.d_VG            =  0.029038476;             % [m    ] distance between VGs in a single pair
vg_case.c_VG            =  0.031058285;             % [m    ] chord of VG  (5cm)
vg_case.h_VG            =  0.010;                   % [m    ] height of VG
vg_case.AoA_VG          =  -15;                     % [deg  ] VG angle of attack (geometric, to free-stream)

vg_case.x0              =    0     * vg_case.h_VG;  % Starting Time (0, or -h_vg/4 , or /3/4*h_vg to enable predevelopment, vortex does not start at TE of Vg) % Exploration = - (3/4) * h_VG;
vg_case.x_end           =    50    * vg_case.h_VG;  % End Time (0.25 = 50 * h_vg)
vg_case.dx              =    0.0000125/8;           % TimeStep Weird stuff happens for dt = 0.1 (time step constraint?)
vg_case.dx_info         =    0.25  * vg_case.h_VG;   % Time between information events
vg_case.friendly_name   =    'Logdberg 2009, JFM paper, 2nd case';

% Now determine rt0 from delta99 (as provided in Wendt)
% BL is flat (ZPG) for IFMF cases so assume Hk=1.4
delta_expdata           = 0.027580938;
delta_99_of_rt          = @(rt) delta_99(vg_case.hk0,rt,vg_case.u_inf/vg_case.nu_inf);
vg_case.rt0             = fzero(@(rt) delta_99_of_rt(rt)-delta_expdata, vg_case.rt0);
% % Run case
snapshots = run_VG_case_alfa41(vg_case);
% Filter out unconverged results (when filaments get out ot domain)
snapshots.N_snapshots = size(snapshots.u_tilde_data, 3);
u_invalidity_threshold = 3 * vg_case.u_inf;
n_snapshot = 1;
n_snapshot_max = snapshots.N_snapshots;
while n_snapshot < snapshots.N_snapshots
    if max(max(abs(snapshots.u_tilde_data(:,:,n_snapshot)))) > u_invalidity_threshold
        n_snapshot_max = n_snapshot-1;
    end
    n_snapshot = n_snapshot + 1;
end
snapshots_filtered = snapshots;
snapshots_filtered.u_tilde_data = snapshots_filtered.u_tilde_data(:,:,1:n_snapshot_max);
snapshots_filtered.x_data = snapshots_filtered.x_data(1:n_snapshot_max);
snapshots_filtered.hk_data = snapshots_filtered.hk_data(1:n_snapshot_max);
snapshots_filtered.theta_data = snapshots_filtered.theta_data(1:n_snapshot_max);
snapshots_filtered.ctau_data = snapshots_filtered.ctau_data(1:n_snapshot_max);
snapshots_filtered.ue_data = snapshots_filtered.ue_data(1:n_snapshot_max);
snapshots_filtered.rt_data = snapshots_filtered.rt_data(1:n_snapshot_max);
snapshots_filtered.y_v_over_d0_data = snapshots_filtered.y_v_over_d0_data(1:n_snapshot_max);
snapshots_filtered.z_v_over_S_data = snapshots_filtered.z_v_over_S_data(1:n_snapshot_max);
snapshots_filtered.gamma_v_data = snapshots_filtered.gamma_v_data(1:n_snapshot_max);
snapshots_filtered.sigma_v_data = snapshots_filtered.sigma_v_data(1:n_snapshot_max);
snapshots_filtered.t_v_data = snapshots_filtered.t_v_data(1:n_snapshot_max);
snapshots_filtered.N_snapshots = size(snapshots_filtered.u_tilde_data, 3);
% Now replace original structure
snapshots = snapshots_filtered;
% % Plot results
n_snapshot = 1 * floor(length(snapshots.x_data) / 4); subplot(221); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 2 * floor(length(snapshots.x_data) / 4); subplot(222); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 3 * floor(length(snapshots.x_data) / 4); subplot(223); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 4 * floor(length(snapshots.x_data) / 4); subplot(224); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
% % Save data
save('datagen/snapshotmats/L002.mat', 'snapshots')


%% Set inputs for Logdberg 2009 (JFM paper 3rd case)
vg_case.rho_inf         = 1.2256;                   % [kg/m3]   Free Stream Density
vg_case.nu_inf          = 1.48e-5;                  % [m2 s-1]  Kinematic Viscosity (nu = mu / rho)
vg_case.u_inf           = 26.50;                    % [m/s]     Unperturbed Free Stream Speed (U_inf!=Ue)

vg_case.hk0             = 1.41;                     % [adim.] Shape Factor (guesstimate based on fact that IFMF has ZPG section in region of VG application!)
vg_case.rt0             = 2500;                     % [adim.] Reynolds Theta of Starting Flow (temporary, initial guess for solver that will determine Rt from delta and hk guess)

vg_case.D_VG            =  0.15;                   % [m    ] distance between center of vortex pairs
vg_case.d_VG            =  0.051969256;             % [m    ] distance between VGs in a single pair
vg_case.c_VG            =  0.055904914;             % [m    ] chord of VG  (5cm)
vg_case.h_VG            =  0.018;                   % [m    ] height of VG
vg_case.AoA_VG          =  -15;                     % [deg  ] VG angle of attack (geometric, to free-stream)

vg_case.x0              =    0     * vg_case.h_VG;  % Starting Time (0, or -h_vg/4 , or /3/4*h_vg to enable predevelopment, vortex does not start at TE of Vg) % Exploration = - (3/4) * h_VG;
vg_case.x_end           =    50    * vg_case.h_VG;  % End Time (0.25 = 50 * h_vg)
vg_case.dx              =    0.0000125/16;           % TimeStep Weird stuff happens for dt = 0.1 (time step constraint?)
vg_case.dx_info         =    0.25  * vg_case.h_VG;   % Time between information events
vg_case.friendly_name   =    'Logdberg 2009, JFM paper, 3rd case';

% Now determine rt0 from delta99 (as provided in Wendt)
% BL is flat (ZPG) for IFMF cases so assume Hk=1.4
delta_expdata           = 0.027580938;
delta_99_of_rt          = @(rt) delta_99(vg_case.hk0,rt,vg_case.u_inf/vg_case.nu_inf);
vg_case.rt0             = fzero(@(rt) delta_99_of_rt(rt)-delta_expdata, vg_case.rt0);
% % Run case
snapshots = run_VG_case_alfa41(vg_case);
% Filter out unconverged results (when filaments get out ot domain)
snapshots.N_snapshots = size(snapshots.u_tilde_data, 3);
u_invalidity_threshold = 3 * vg_case.u_inf;
n_snapshot = 1;
n_snapshot_max = snapshots.N_snapshots;
while n_snapshot < snapshots.N_snapshots
    if max(max(abs(snapshots.u_tilde_data(:,:,n_snapshot)))) > u_invalidity_threshold
        n_snapshot_max = n_snapshot-1;
    end
    n_snapshot = n_snapshot + 1;
end
snapshots_filtered = snapshots;
snapshots_filtered.u_tilde_data = snapshots_filtered.u_tilde_data(:,:,1:n_snapshot_max);
snapshots_filtered.x_data = snapshots_filtered.x_data(1:n_snapshot_max);
snapshots_filtered.hk_data = snapshots_filtered.hk_data(1:n_snapshot_max);
snapshots_filtered.theta_data = snapshots_filtered.theta_data(1:n_snapshot_max);
snapshots_filtered.ctau_data = snapshots_filtered.ctau_data(1:n_snapshot_max);
snapshots_filtered.ue_data = snapshots_filtered.ue_data(1:n_snapshot_max);
snapshots_filtered.rt_data = snapshots_filtered.rt_data(1:n_snapshot_max);
snapshots_filtered.y_v_over_d0_data = snapshots_filtered.y_v_over_d0_data(1:n_snapshot_max);
snapshots_filtered.z_v_over_S_data = snapshots_filtered.z_v_over_S_data(1:n_snapshot_max);
snapshots_filtered.gamma_v_data = snapshots_filtered.gamma_v_data(1:n_snapshot_max);
snapshots_filtered.sigma_v_data = snapshots_filtered.sigma_v_data(1:n_snapshot_max);
snapshots_filtered.t_v_data = snapshots_filtered.t_v_data(1:n_snapshot_max);
snapshots_filtered.N_snapshots = size(snapshots_filtered.u_tilde_data, 3);
% Now replace original structure
snapshots = snapshots_filtered;
% % Plot results
n_snapshot = 1 * floor(length(snapshots.x_data) / 4); subplot(221); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 2 * floor(length(snapshots.x_data) / 4); subplot(222); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 3 * floor(length(snapshots.x_data) / 4); subplot(223); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
n_snapshot = 4 * floor(length(snapshots.x_data) / 4); subplot(224); surf(snapshots.z_mesh/snapshots.S, snapshots.y_mesh/snapshots.delta0, snapshots.u_tilde_data(:,:,n_snapshot)); view(2); shading flat;
% % Save data
save('datagen/snapshotmats/L003.mat', 'snapshots')





