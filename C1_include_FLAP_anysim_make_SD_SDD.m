%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       OptiFLOW - Parallel MultiObjective Airfoil Optimization System
%           Gael de Oliveira, 2010-2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       File: 
%           Multiobjective Optimization Example Case for Airfoils with
%           Flaps. 
%           Arbitrary number of polar-configurations per design point.
%
%           Aerodynamic Goal: Glide ration L/D over a range of angles of 
%                             attack and flap angles
%           Structural  Goal: Maximum Airfoil Building Height
%
%           Parallel Execution Enabled
%
%           Include file for: C0_start_example_multisim_parallel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Generate Flapping Angle range vector
N_flap_steps =  2 * N_flap_steps_per_side + 1;          % Always odd to provide for undeflected case in all situations!
% Introduce robustness for case of N_flap_steps_per_side = 0;
if N_flap_steps > 1 
   flap_angle_scalings = linspace(-1,1, N_flap_steps);     
else
   N_flap_steps = 1;
   flap_angle_scalings = [0]; %#ok<NBRAK>
end
   

% Preallocate cell vectors of shape definition and shape dynamizer elements
SD_cell_vector  = cell(size(flap_angle_scalings));
SDD_cell_vector = cell(size(SD_cell_vector));


for n_flap_step = 1:N_flap_steps
    % Create Shape Definition Objects using previously defined parametrizations 
    SD_cell_vector{n_flap_step} = shape_definition_cst(p_upper , p_lower , N_dummy_parameters, ['shape_definition_' num2str(n_flap_step)]);     % String is just for identification! NO functional meaning! For now!
    
    % Instantiate Shape Dynamizer object and hook it to its shape_definition_cst object
    SDD_cell_vector{n_flap_step} = shape_dynamizer(SD_cell_vector{n_flap_step});        
    
    % For each SDD, set a cell vector with a list of fieldnames that are
    % dynamic
    SDD_cell_vector{n_flap_step}.dynamic_variable_list = {'x_hinge' , 'flap_angle'};
     
    % Array with index matching fieldnames to be dynamized
    SDD_cell_vector{n_flap_step}.dynamic_variable_index = [2 , 3];                            % Array with index matching fieldnames with dummy_parameters elements (obtained from SD.breakdown_parameters(parameters))

    % Use a dynamic variable processor to differentiate analyses of a same
    % case
    % Hinge is set at its spot, so first function is identity, while to
    % flap angles are scaled into steps to make a 2d polar case (Cl
    % vs(alpha, beta)
    SDD_cell_vector{n_flap_step}.dynamic_variable_processors = {@(x) x , @(x) x*flap_angle_scalings(n_flap_step)};
end

% Now make an index pointer to the undeflected case, and name the undeflected 
% shape definition SD stbe used by the
% shape fit object!
%(name matters for backward compatibility with GUI, but not for SF)
index_undeflected = N_flap_steps_per_side + 1;
SD = SD_cell_vector{index_undeflected};


%% Now proceed to create simulation worker objects list
% (based on the previously defined simulation protocol SP )

% Preallocate cell vector of simulation workers
SW_cell_vector  = cell(size(SD_cell_vector));

for n_flap_step = 1:N_flap_steps
    SW_cell_vector{n_flap_step} = simulation_worker(['simul_worker_' num2str(n_flap_step)], SP, [] , SD_cell_vector{n_flap_step},  SC);
%  ---- Edit any properties of choice ---- %
    SW_cell_vector{n_flap_step}.app_name = 'RBINKOLIVEIRA_V2';                               % For xfoil write SW_1.app_name = 'xfoil'
    SW_cell_vector{n_flap_step}.parallelized = 1;                                            % Parallelize simulations for cost function computation (to be effective it requires that vectorize is set to true in the genetic algorithm)
    SW_cell_vector{n_flap_step}.multisim_id_value = flap_angle_scalings(n_flap_step);
%  ---- End of property edition ---- %
end
