function sc = cost_function_turbulence_intensity( ~ , experiments_results)
% Extract last design parameter and return it. This is used in the
% probabilistic design runs presented at the Torque 2018 conference
% 
% This is a positive value that we want to maximize. Invalidity is
% manifested by returning a zero.

% Input extraction
% parameters can be ignored

% Extract experiment results from cell capsule
experiment_result   = experiments_results{1};
% Get genotype of experiment
x = experiment_result.x;

%% Extract last coordinate (equal to turbulence intensity)
sc = x(end);

%% Make check to verify that sides do not cross
valid_shape = valid_shape_test(experiment_result);
if not(valid_shape)
   sc = 0; 
end

end