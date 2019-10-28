%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       OptiFLOW - Parallel MultiObjective Airfoil Optimization System
%           Gael de Oliveira, 2010-2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       File: 
%           Airfoil Flapper / Get an airfoil file, flap it, and export
%           "smoothly" flapped coordinates.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Make a Clean Sheet
clear all; close all; clc

%% Numerical Inputs
x_hinge               = 0.75;                             % Over chord
tc_hinge              = 0.5 ;                             % Relative (as in xfoil definition)
flap_angles           = [-25 -20 -15 -10 -5 0 5 10];      % In degrees
N_points              = 160;                              % Number of points to define airfoil shape
doubleIO              = 0;                                % False for single precision output, True for double!
Re                    = 3e5;

%% File Inputs
airfoil_filename      = 'S1223RTL.p160.dat';
%airfoil_filename      = 'RTL50.clean.p160.dat';
%airfoil_filename      = 'MH60_clean.p160.dat';
airfoil_folder        = '/Users/gael/Desktop/eKite/airfoils/';
airfoil_export_prefix = 'S1223RTL';
%airfoil_export_prefix = 'RTL50';
%airfoil_export_prefix = 'MH60';
airfoil_export_suffix = '.air';


%% System Work necessary to start
% Start by adding the source folders to the matlab path
fs = filesep();      % Folder separator is OS dependent
addpath([cd fs 'src']);
addpath([cd fs 'user_src']);
addpath([cd fs 'gui']);

% Create System Context Object and Set Context
SC = system_context; SC.N_cores = 1; SC.set_context;

%% Instanciate geometric manipulations classes
% Parametrization objects for upper and lower side
p_upper = parametrization( 'cst_upper'  , 8);
p_lower = parametrization( 'cst_lower'  , 8);
 
% Fix number of non-CST-shape parameters
N_dummy_parameters = 3;                                            % Number of parameters that are not shape definition ones for the three flap cases (undeflected, downward and upward)
N_flap_steps_per_side = 1;                                         % Number of flap steps per side 

% Shape Definition Objects using previously defined parametrizations 
SD_0    = shape_definition_cst(p_upper , p_lower , N_dummy_parameters, 'no_flap');     % String is just for identification! NO functional meaning! For now!
SD_flap = shape_definition_cst(p_upper , p_lower , N_dummy_parameters, 'flapped');     % String is just for identification! NO functional meaning! For now!

% Make Shape Fit Object
SF = shape_fit_cst(SD_0, 'cst88_undeflected_flap');                  % 'fitcst88_undeflected_flap' is any arbitrary name you like

% Shape Dynamizer object hooked to SD_flap object
SDD = shape_dynamizer(SD_flap);

%% Load Airfoil
fileToRead = [airfoil_folder ,  airfoil_filename];
x_airfoil = SF.get_parameters_from_file(fileToRead);

%% Regenerate Coordinates of Original Airfoil
% First, for unflapped airfoils
[tx_0, tz_0] = SD_0.generate_coordinates(160, x_airfoil);
% Open a plot
%plot(tx_0, tz_0); axis equal; grid on; hold on; xlabel('x/c'); ylabel('z/c');
% Determine z_hinge
[thickness_at_hinge, camber_at_hinge] = SD_0.thickness_camber_at_tx(x_hinge, x_airfoil);
z_hinge = camber_at_hinge + thickness_at_hinge * (tc_hinge - 0.5);

%% Flap Airfoil
tx_cell     = cell(size(flap_angles));
tz_cell     = cell(size(flap_angles));
legend_cell = cell(size(flap_angles));
figure(1)

% Now run over flap settings
for n_flap = 1:length(flap_angles)
    % First update flap and hinge properties
    SDD.x_hinge    =  x_hinge;
    SDD.z_hinge    =  z_hinge;
    SDD.flap_angle =  flap_angles(n_flap);
    
    % Then generate flapped coordinates
    [tx_cell{n_flap}, tz_cell{n_flap}] = SD_flap.generate_coordinates(N_points, x_airfoil);
    
    % Plot this out 
    plot(tx_cell{n_flap}, tz_cell{n_flap}); axis equal; grid on; hold on; xlabel('x/c'); ylabel('z/c');
    
    % Compose and display output filename
    if SDD.flap_angle < 0 ; flap_sign = 'm'; else flap_sign = ''; end
    airfoil_export_filename = [airfoil_folder, airfoil_export_prefix, '_F' , flap_sign , num2str(abs(SDD.flap_angle)), airfoil_export_suffix];
    disp(['Write to: ' , airfoil_export_filename]);
    
    % Write coordinates to file
    Xt = [tx_cell{n_flap}, tz_cell{n_flap}]; %#ok<NASGU>
    
    % Export to file
    if (doubleIO > 0)
        save( airfoil_export_filename , 'Xt','-ascii', '-double');
    else
        save( airfoil_export_filename , 'Xt','-ascii');
    end
    
    legend_cell{n_flap} = [airfoil_export_prefix, '_F' , flap_sign , num2str(abs(SDD.flap_angle))];
end

% Add legend
legend(legend_cell);
title(['Smooth Flap Airfoil Shapes (Knee=0.1, CST/Kulfan/Oliveira Parametrization)']);
orient landscape
print([airfoil_folder airfoil_export_prefix '_airfoil_plot.pdf'], '-dpdf')


%% Make simulation protocol and worker objects

% Make a Simulation Protocol Object
SP = simulation_protocol('forced_transition' , SC);
SP.target_application='RBINKOLIVEIRA_V2';                  % To use Xfoil as your application set this to 'xfoil'
% Set simulation type and angle step
SP.operation = 'alfa_polar_ref_start';                     % Set the operation mode. 'alfa_polar_ref_start' makes a polar in two sweeps and is usually a good choice 
SP.operation_parameters = [-15 30 0.2 0];                  % Make polar from -15 to 20 degrees in 0.5 degree increments, in two sweeps from 0 to 20 and -0.5 to -15
SP.Re = Re;                                               % Reynolds Number (NaN means inviscid)
SP.xtr_top=.15;                                            % set forced transition on top
SP.xtr_bot=.15;                                            % set forced transition on bottom
SP.N_crit=9;                                               % set amplification to default value (Van Ingen)

% Make simulation worker objects
SW_0    = simulation_worker('no_flap', SP, [] , SD_0   ,  SC);
SW_flap = simulation_worker('flapped', SP, [] , SD_flap,  SC);

%% Airfoil Polars
SDD.flap_angle = -25;
if SDD.flap_angle < 0 ; flap_sign = 'm'; else flap_sign = ''; end
airfoil_export_filename = [airfoil_folder, airfoil_export_prefix, '_F' , flap_sign , num2str(abs(SDD.flap_angle)), airfoil_export_suffix];
ap_flap1 = SW_flap.run_polar_on_airfoil_file(airfoil_export_filename);
current_polar_file = [SW_flap.SC.tmp_subdir  SW_flap.SC.core_subdirs{1} SW_flap.SC.polar_filename];

SDD.flap_angle = -20;
if SDD.flap_angle < 0 ; flap_sign = 'm'; else flap_sign = ''; end
airfoil_export_filename = [airfoil_folder, airfoil_export_prefix, '_F' , flap_sign , num2str(abs(SDD.flap_angle)), airfoil_export_suffix];
ap_flap2 = SW_flap.run_polar_on_airfoil_file(airfoil_export_filename);

SDD.flap_angle = -15;
if SDD.flap_angle < 0 ; flap_sign = 'm'; else flap_sign = ''; end
airfoil_export_filename = [airfoil_folder, airfoil_export_prefix, '_F' , flap_sign , num2str(abs(SDD.flap_angle)), airfoil_export_suffix];
ap_flap3 = SW_flap.run_polar_on_airfoil_file(airfoil_export_filename);

SDD.flap_angle = -10;
if SDD.flap_angle < 0 ; flap_sign = 'm'; else flap_sign = ''; end
airfoil_export_filename = [airfoil_folder, airfoil_export_prefix, '_F' , flap_sign , num2str(abs(SDD.flap_angle)), airfoil_export_suffix];
ap_flap4 = SW_flap.run_polar_on_airfoil_file(airfoil_export_filename);

SDD.flap_angle = - 5;
if SDD.flap_angle < 0 ; flap_sign = 'm'; else flap_sign = ''; end
airfoil_export_filename = [airfoil_folder, airfoil_export_prefix, '_F' , flap_sign , num2str(abs(SDD.flap_angle)), airfoil_export_suffix];
ap_flap5 = SW_flap.run_polar_on_airfoil_file(airfoil_export_filename);

SDD.flap_angle =   0;
if SDD.flap_angle < 0 ; flap_sign = 'm'; else flap_sign = ''; end
airfoil_export_filename = [airfoil_folder, airfoil_export_prefix, '_F' , flap_sign , num2str(abs(SDD.flap_angle)), airfoil_export_suffix];
ap_flap6 = SW_flap.run_polar_on_airfoil_file(airfoil_export_filename);

SDD.flap_angle =   5;
if SDD.flap_angle < 0 ; flap_sign = 'm'; else flap_sign = ''; end
airfoil_export_filename = [airfoil_folder, airfoil_export_prefix, '_F' , flap_sign , num2str(abs(SDD.flap_angle)), airfoil_export_suffix];
ap_flap7 = SW_flap.run_polar_on_airfoil_file(airfoil_export_filename);

SDD.flap_angle =   10;
if SDD.flap_angle < 0 ; flap_sign = 'm'; else flap_sign = ''; end
airfoil_export_filename = [airfoil_folder, airfoil_export_prefix, '_F' , flap_sign , num2str(abs(SDD.flap_angle)), airfoil_export_suffix];
ap_flap8 = SW_flap.run_polar_on_airfoil_file(airfoil_export_filename);


%% Plots
% Make color map
cmap = [ 0         0.4470    0.7410; ...
         0.8500    0.3250    0.0980; ...
         0.9290    0.6940    0.1250; ...
         0.4940    0.1840    0.5560; ...
         0.4660    0.6740    0.1880; ...
         0.3010    0.7450    0.9330; ...
         0.6350    0.0780    0.1840; ...
         0         0.4470    0.7410; ];

% % Polar Plots
figure(2)
ax_1 = ap_flap1.plot;
ap_flap2.plot(ax_1 , cmap(2,:));
ap_flap3.plot(ax_1 , cmap(3,:));
ap_flap4.plot(ax_1 , cmap(4,:));
ap_flap5.plot(ax_1 , cmap(5,:));
ap_flap6.plot(ax_1 , cmap(6,:));
ap_flap7.plot(ax_1 , cmap(7,:));
ap_flap8.plot(ax_1 , cmap(8,:));
legend(ax_1(2), 'off')
%legend(ax_1(1), {'-25', '-20', '-15', '-10', '-5', '0', '5', '10'})
legend(ax_1(1), {'-25', '-20', '', '-15', '', '-10', '', '-5', '', '0', '', '5', '', '10'});
title([airfoil_export_prefix ' - Airfoil Polars (Re=' num2str(SP.Re) ' , xtr=0.15)']);
orient landscape
print([airfoil_folder airfoil_export_prefix '_polar_plot_Re' num2str(SP.Re) '.pdf'], '-dpdf')

% % Flap Plots
figure(3)
aoa = - 5; plot(flap_angles, [ap_flap1.cl_alpha(aoa) , ap_flap2.cl_alpha(aoa), ap_flap3.cl_alpha(aoa), ap_flap4.cl_alpha(aoa), ap_flap5.cl_alpha(aoa), ap_flap6.cl_alpha(aoa), ap_flap7.cl_alpha(aoa), ap_flap8.cl_alpha(aoa)], 'o-'); grid on; ylabel('Cl'); xlabel('Flap Deflection (deg)'); hold on
aoa =   0; plot(flap_angles, [ap_flap1.cl_alpha(aoa) , ap_flap2.cl_alpha(aoa), ap_flap3.cl_alpha(aoa), ap_flap4.cl_alpha(aoa), ap_flap5.cl_alpha(aoa), ap_flap6.cl_alpha(aoa), ap_flap7.cl_alpha(aoa), ap_flap8.cl_alpha(aoa)], 'o-'); grid on; ylabel('Cl'); xlabel('Flap Deflection (deg)'); hold on
aoa =   5; plot(flap_angles, [ap_flap1.cl_alpha(aoa) , ap_flap2.cl_alpha(aoa), ap_flap3.cl_alpha(aoa), ap_flap4.cl_alpha(aoa), ap_flap5.cl_alpha(aoa), ap_flap6.cl_alpha(aoa), ap_flap7.cl_alpha(aoa), ap_flap8.cl_alpha(aoa)], 'o-'); grid on; ylabel('Cl'); xlabel('Flap Deflection (deg)'); hold on
aoa =  10; plot(flap_angles, [ap_flap1.cl_alpha(aoa) , ap_flap2.cl_alpha(aoa), ap_flap3.cl_alpha(aoa), ap_flap4.cl_alpha(aoa), ap_flap5.cl_alpha(aoa), ap_flap6.cl_alpha(aoa), ap_flap7.cl_alpha(aoa), ap_flap8.cl_alpha(aoa)], 'o-'); grid on; ylabel('Cl'); xlabel('Flap Deflection (deg)'); hold on
legend('aoa = -5deg', 'aoa = 0deg', 'aoa = 5deg', 'aoa = 10deg');
title([airfoil_export_prefix  ' - Flap Response (Re=' num2str(SP.Re) ' , xtr=0.15)']);
orient landscape
print([airfoil_folder airfoil_export_prefix '_flap_response_Re' num2str(SP.Re) '.pdf'], '-dpdf')

figure(4)
ap_flap6.plot(); title([airfoil_export_prefix ' - Airfoil Polars (Re=' num2str(SP.Re) ' , xtr=0.15)']);
orient landscape
print([airfoil_folder airfoil_export_prefix '_single_polar_plot_Re' num2str(SP.Re) '.pdf'], '-dpdf')










