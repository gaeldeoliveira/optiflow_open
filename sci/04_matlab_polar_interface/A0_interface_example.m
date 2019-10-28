
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       RFOILSUC Interface (Stripped down from the optimization code)
%
%           TU-Delft, Gael de Oliveira, 2010-2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
%% Create Simulation Objects (Context, Protocol and Worker)
% Start by adding the source folders to the matlab path
fs = filesep();                 % Folder separator is OS dependent
addpath([cd fs 'src']);         % Add optimizer source code folder

% Create System Context Object and Set Context
SC = system_context(); SC.set_context;
% Make a Simulation Protocol Object
SP = simulation_protocol('single_point' , SC);
% Make a Simulation Worker Object
SW = simulation_worker('simul_worker_1', SP, [] , [],  SC);  % Mariline
SW.app_name = '../../bin/rfoilsucplasma';                              % Path to RFOIL executable


%% Set Plasma Actuator Parameters
% PL = plasma_actuator_parameters();
% % Set FIelds to desired values!
% PL.tpp = 0.006;
% PL.lpp = 0.1950;
% PL.cftp = 0;

%% Set Aerodynamic Conditions
% Operation
SP.operation = 'alfa_polar_ref_start';                          % Set the operation mode. 'alfa_polar_ref_start' makes a polar in two sweeps and is usually a good choice 
SP.operation_parameters = [-15 20 0.5 0];                       % Make polar from -15 to 20 degrees in 0.5 degree increments, in two sweeps from 0 to 20 and -0.5 to -15

% Hook Simulation Protocol to Plasma Object so that DBD Actuator is Active
% SP.plasma_actuator_type = 'Single_Plasma';
% SP.plasma_actuator_parameters = PL;

% Operating Point
SP.Re       = 1.6e6;                                            % Reynolds Number (NaN means inviscid)
SP.M        = 0;                                                % Mach number (all validation made for incompressible! NO warranties otherwise!)

% Turbulence and Transition 
SP.N_crit   = 9;                                                % Critical Amplification Factor (for transition)
SP.xtr_top  = 0.99;                                             % x/c of forced transition point for top side (1 = free)
SP.xtr_bot  = 0.99;                                             % x/c of forced transition point for bottom side (1 = free)

% Rotational Effects
SP.cr = 0;                                                      % c/r (for rotational effects) 
SP.cr_factor = 2/3;                                             % c/r correction factor
        
% Numerical Parameters         
SP.save_bl_parameters = 'false';                                % If true, boundary layer parameters are saved (currently only under rfoilsucblind)
SP.ITER     = 300;                                              % Allowable iterations per operating point
SP.repanel = 'false';                                           % Decided wether RFOILs panels himself or uses file discretization! repanel = 'false' or 'true'
SP.EPSV = 0;                                                    % EPS1 convergence criterium (EPSV = 0 keeps everything in default 1e-4 on most versions... double precision allows to go to 1e-7/1e-8 very safely and eventually to 1e-12)
     
%% Make Polar
% Choose airfoil filename
airfoil_file = 'du25.air';                                      % Path to Airfoil
% Run Polar
ap = SW.run_polar_on_airfoil_file(airfoil_file);
% For failed polars, ap=NaN, for sucessful polars, ap is an
% aerodynamic_polar object, able to plot and share data!

%% Plot Polar and Get Data
ap.plot();

% Get Coefficients
alpha_polar = ap.alpha_range();
cl_polar    = ap.cl_alpha(alpha_polar);
cd_polar    = ap.cd_alpha(alpha_polar);
cm_polar    = ap.cm_alpha(alpha_polar);

% Print plot to eps!
print('-dpdf', [airfoil_file '.pdf'])


