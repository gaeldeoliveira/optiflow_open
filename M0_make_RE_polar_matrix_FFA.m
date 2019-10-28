%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       OptiFLOW - Parallel MultiObjective Airfoil Optimization System
%           Gael de Oliveira, Ricardo Pereira, 2010-2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       File: 
%           Make a matrix of airfoil polars (alpha, Re)
%           
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Clean environment
close all; clear all; clc; %#ok<CLALL>

%% Polar conditions
case_overview.polar_range    = [-5 20 0.2 0];       % Make polar from -15 to 20 degrees in 0.5 degree increments, in two sweeps from 0 to 20 and -0.5 to -15
case_overview.polar_range2   = [-4 18 0.2 0];       % Make polar from -15 to 20 degrees in 0.5 degree increments, in two sweeps from 0 to 20 and -0.5 to -15
case_overview.Re_list        = [3e6, 9e6, 18e6];    % Reynolds Number (NaN means inviscid)
case_overview.N_crit         =   9   ;              % critical amplification factor for Tollmien-Schlichting waves
case_overview.thick_only     = true

%% Make 21% airfoil polar_mesh
if not(case_overview.thick_only == true)
    case_overview.filename       = 'rotor_integration/airfoil_families/FFA/FFAw3211.txt';   % Name of airfoil file (can also have folders, absolute or from cd (./))
    % Make rough polar
    case_overview.xtr_top        =   0.05;              % set top    transition position for free
    case_overview.xtr_bot        =   0.10;              % set bottom transition position for free
    polar_mesh_trip = make_polar_mesh(case_overview); polar_mesh = polar_mesh_trip; %#ok<NASGU>
    save([case_overview.filename(1:end-4), '_polar_mesh_trip.mat'], 'polar_mesh');
    % Make free transition polar
    case_overview.xtr_top        =   0.99;              % set top    transition position for free
    case_overview.xtr_bot        =   0.99;              % set bottom transition position for free
    polar_mesh_free = make_polar_mesh(case_overview); polar_mesh = polar_mesh_free; %#ok<NASGU>
    save([case_overview.filename(1:end-4), '_polar_mesh_free.mat'], 'polar_mesh');
end

%% Make 24% airfoil polar_mesh
if not(case_overview.thick_only == true)
    case_overview.filename       = 'rotor_integration/airfoil_families/FFA/FFAw3241.txt';   % Name of airfoil file (can also have folders, absolute or from cd (./))
    % Make rough polar
    case_overview.xtr_top        =   0.05;              % set top    transition position for free
    case_overview.xtr_bot        =   0.10;              % set bottom transition position for free
    polar_mesh_trip = make_polar_mesh(case_overview); polar_mesh = polar_mesh_trip; %#ok<NASGU>
    save([case_overview.filename(1:end-4), '_polar_mesh_trip.mat'], 'polar_mesh');
    % Make free transition polar
    case_overview.xtr_top        =   0.99;              % set top    transition position for free
    case_overview.xtr_bot        =   0.99;              % set bottom transition position for free
    polar_mesh_free = make_polar_mesh(case_overview); polar_mesh = polar_mesh_free; %#ok<NASGU>
    save([case_overview.filename(1:end-4), '_polar_mesh_free.mat'], 'polar_mesh');
end

%% Make 30% airfoil polar_mesh
if not(case_overview.thick_only == true)
    case_overview.filename       = 'rotor_integration/airfoil_families/FFA/FFAw3301.txt';   % Name of airfoil file (can also have folders, absolute or from cd (./))
    % Make rough polar
    case_overview.xtr_top        =   0.05;              % set top    transition position for free
    case_overview.xtr_bot        =   0.10;              % set bottom transition position for free
    polar_mesh_trip = make_polar_mesh(case_overview); polar_mesh = polar_mesh_trip; %#ok<NASGU>
    save([case_overview.filename(1:end-4), '_polar_mesh_trip.mat'], 'polar_mesh');
    % Make free transition polar
    case_overview.xtr_top        =   0.99;              % set top    transition position for free
    case_overview.xtr_bot        =   0.99;              % set bottom transition position for free
    polar_mesh_free = make_polar_mesh(case_overview); polar_mesh = polar_mesh_free; %#ok<NASGU>
    save([case_overview.filename(1:end-4), '_polar_mesh_free.mat'], 'polar_mesh');
end

%% Make 36% airfoil polar_mesh
if not(case_overview.thick_only == true)
    case_overview.filename       = 'rotor_integration/airfoil_families/FFA/FFAw3360.txt';   % Name of airfoil file (can also have folders, absolute or from cd (./))
    % Make rough polar
    case_overview.xtr_top        =   0.05;              % set top    transition position for free
    case_overview.xtr_bot        =   0.10;              % set bottom transition position for free
    polar_mesh_trip = make_polar_mesh(case_overview); polar_mesh = polar_mesh_trip; %#ok<NASGU>
    save([case_overview.filename(1:end-4), '_polar_mesh_trip.mat'], 'polar_mesh');
    % Make free transition polar
    case_overview.xtr_top        =   0.99;              % set top    transition position for free
    case_overview.xtr_bot        =   0.99;              % set bottom transition position for free
    polar_mesh_free = make_polar_mesh(case_overview); polar_mesh = polar_mesh_free; %#ok<NASGU>
    save([case_overview.filename(1:end-4), '_polar_mesh_free.mat'], 'polar_mesh');
end

%% Make 48% airfoil polar_mesh
if case_overview.thick_only == true
    case_overview.filename       = 'rotor_integration/airfoil_families/FFA/FFAw3480.txt';   % Name of airfoil file (can also have folders, absolute or from cd (./))
    % Make rough polar
    case_overview.xtr_top        =   0.05;              % set top    transition position for free
    case_overview.xtr_bot        =   0.10;              % set bottom transition position for free
    polar_mesh_trip = make_polar_mesh(case_overview); polar_mesh = polar_mesh_trip; %#ok<NASGU>
    save([case_overview.filename(1:end-4), '_polar_mesh_trip.mat'], 'polar_mesh');
    % Make free transition polar
    case_overview.xtr_top        =   0.99;              % set top    transition position for free
    case_overview.xtr_bot        =   0.99;              % set bottom transition position for free
    polar_mesh_free = make_polar_mesh(case_overview); polar_mesh = polar_mesh_free; %#ok<NASGU>
    save([case_overview.filename(1:end-4), '_polar_mesh_free.mat'], 'polar_mesh');
end

%% Make 60% airfoil polar_mesh
if case_overview.thick_only == true
    case_overview.filename       = 'rotor_integration/airfoil_families/FFA/FFAw3600.txt';   % Name of airfoil file (can also have folders, absolute or from cd (./))
    % Make rough polar
    case_overview.xtr_top        =   0.05;              % set top    transition position for free
    case_overview.xtr_bot        =   0.10;              % set bottom transition position for free
    polar_mesh_trip = make_polar_mesh(case_overview);   % fix convergence issue manually
    polar_mesh_trip.cl_mesh(:,2) = polar_mesh_trip.cl_mesh(:,1);
    polar_mesh_trip.cl_mesh(:,3) = polar_mesh_trip.cl_mesh(:,1);
    polar_mesh_trip.cd_mesh(:,2) = polar_mesh_trip.cd_mesh(:,1);
    polar_mesh_trip.cd_mesh(:,3) = polar_mesh_trip.cd_mesh(:,1);
    polar_mesh_trip.cm_mesh(:,2) = polar_mesh_trip.cm_mesh(:,1);
    polar_mesh_trip.cm_mesh(:,3) = polar_mesh_trip.cm_mesh(:,1);
    polar_mesh = polar_mesh_trip; %#ok<NASGU>
    save([case_overview.filename(1:end-4), '_polar_mesh_trip.mat'], 'polar_mesh');
    % Make free transition polar
    case_overview.xtr_top        =   0.99;              % set top    transition position for free
    case_overview.xtr_bot        =   0.99;              % set bottom transition position for free
    polar_mesh_free = make_polar_mesh(case_overview); polar_mesh = polar_mesh_free; %#ok<NASGU>
    save([case_overview.filename(1:end-4), '_polar_mesh_free.mat'], 'polar_mesh')
end
    
%% Make 100% airfoil polar_mesh
case_overview.filename       = 'rotor_integration/airfoil_families/FFA/cylinder.txt';   % Name of airfoil file (can also have folders, absolute or from cd (./))
% Make polar (very simple empirical fit for now)
polar_mesh = make_polar_mesh_for_cylinder(case_overview);
% Save for trip
save([case_overview.filename(1:end-4), '_polar_mesh_trip.mat'], 'polar_mesh');
% Save for free
save([case_overview.filename(1:end-4), '_polar_mesh_free.mat'], 'polar_mesh');


%% Function to make polar meshes
function polar_mesh = make_polar_mesh(case_overview)
    % % Initialize environment
    % Add access to necessary paths
    fs = filesep(); addpath([cd fs 'src']);
    
    % Create and set system context
    SC = system_context;
    SC.N_cores = 1;
    SC.set_context;
    
    % Make simulation protocol object
    SP1 = simulation_protocol('free_transition' , SC); SP1.target_application='RBINKOLIVEIRA_V2';
    SP1.operation = 'alfa_polar_ref_start';                                         % Set the operation mode. 'alfa_polar_ref_start' makes a polar in two sweeps and is usually a good choice
    SP1.operation_parameters = case_overview.polar_range;                           % Make polar from -15 to 20 degrees in 0.5 degree increments, in two sweeps from 0 to 20 and -0.5 to -15
    SP1.xtr_top              = case_overview.xtr_top;                          % set forced transition on top
    SP1.xtr_bot              = case_overview.xtr_bot;                          % set forced transition on bottom
    SP1.N_crit               = case_overview.N_crit;                                % critical amplification factor for Tollmien-Schlichting waves
    SP1.save_bl_parameters   = 'true';
    
    % Make simulation worker object
    SW_1 = simulation_worker('simul_worker_1', SP1, [] , [],  SC);                      % Mariline
    SW_1.app_name            = 'RBINKOLIVEIRA_V2';                                  % For xfoil write SW_1.app_name = 'xfoil'
    SW_1.fig_active          = 0;                                                   % Activate plotting for clean case
    SW_1.parallelized        = 1;
    
    
    % % Now run polars for each demanded Re
    ap_cell = cell(length(case_overview.Re_list), 1);
    
    for n_Re = 1:length(case_overview.Re_list)
        % Set Reynolds
        SP1.Re = case_overview.Re_list(n_Re);
        % Aggregate polars into cell array
        ap_cell{n_Re} = SW_1.run_polar_on_airfoil_file(case_overview.filename);
    end
    
    % % And aggregate everything
    
    % Make reference list of angles of attack
    alpha_range = case_overview.polar_range2(1):case_overview.polar_range2(3):case_overview.polar_range2(2);
    Re_range    = case_overview.Re_list;
    % List of lift, drag and moment coefficients
    cl_mesh  = zeros(length(alpha_range), length(Re_range));
    cd_mesh  = zeros(length(alpha_range), length(Re_range));
    cm_mesh  = zeros(length(alpha_range), length(Re_range));
    
    % Now make sure everything is aggregated together
    for n_Re = 1:length(case_overview.Re_list)
        cl_mesh(:, n_Re) = ap_cell{n_Re}.cl_alpha(alpha_range);
        cd_mesh(:, n_Re) = ap_cell{n_Re}.cd_alpha(alpha_range);
        cm_mesh(:, n_Re) = ap_cell{n_Re}.cm_alpha(alpha_range);
    end
    
    % % Make a structure to keep it all together
    polar_mesh = struct();
    polar_mesh.alpha_range  = alpha_range;
    polar_mesh.Re_range     = Re_range;
    polar_mesh.cl_mesh      = cl_mesh;
    polar_mesh.cd_mesh      = cd_mesh;
    polar_mesh.cm_mesh      = cm_mesh;
end

function polar_mesh = make_polar_mesh_for_cylinder(case_overview)
    % Make reference list of angles of attack
    alpha_range = case_overview.polar_range2(1):case_overview.polar_range2(3):case_overview.polar_range2(2);
    Re_range    = case_overview.Re_list;
    % List of lift, drag and moment coefficients
    cl_mesh  = zeros(length(alpha_range), length(Re_range));
    cd_mesh  = zeros(length(alpha_range), length(Re_range));
    cm_mesh  = zeros(length(alpha_range), length(Re_range));
    % Now fill in
    for n_Re = 1:length(case_overview.Re_list)
        cl_mesh(:, n_Re) = 0;
        cd_mesh(:, n_Re) = interp1([0.5e6, 5e6, 10e6, 100e6], [0.5, 0.8, 0.9, 1.0], case_overview.Re_list(n_Re));
        cm_mesh(:, n_Re) = 0;
    end
    
    % % Make a structure to keep it all together
    polar_mesh = struct();
    polar_mesh.alpha_range  = alpha_range;
    polar_mesh.Re_range     = Re_range;
    polar_mesh.cl_mesh      = cl_mesh;
    polar_mesh.cd_mesh      = cd_mesh;
    polar_mesh.cm_mesh      = cm_mesh;
end




