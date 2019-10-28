function sc = cost_function_expected_CLCD_and_CL_at_expected_AoA(parameters , experiments_results)
%function sc = sub_cost_function_CLCD_at_CL(aero_polar, cl_i , w_i)
%   Calculates weighted Cl/Cd of a polar at Cl/Weight combinations from
%   alpha_i and w_i vectors

%% Experiment extraction
% Extract aerodynamic polar from standard input form
experiment_result   = experiments_results{1};
ap = experiment_result.ap; %#ok<NASGU>

%% Run two single objective cost-functions
ld_bar            = cost_function_expected_CLCD_at_expected_AoA(parameters , experiments_results);
cl_bar            = cost_function_expected_CL_at_expected_AoA(  parameters , experiments_results);

%% Extract normalization parameters from parameter structure
ld_bar_ref        = parameters.ld_bar_ref;       % Determined at beginning of process
cl_bar_ref        = parameters.cl_bar_ref;       % Determined at beginning of process

%% Normalize objective function
ld_bar_normalized = ld_bar / ld_bar_ref;
cl_bar_normalized = cl_bar / cl_bar_ref;

%% Make power factors
if ld_bar_normalized > 1; k_ld_bar = 1; else; k_ld_bar = 2; end
if cl_bar_normalized > 1; k_cl_bar = 1; else; k_cl_bar = 2; end

%% Make penalties
ld_bar_penalty    = ((1 - ld_bar_normalized) * 100 )^k_ld_bar;
cl_bar_penalty    = ((1 - cl_bar_normalized) * 100 )^k_cl_bar;

%% Combine parameters
sc = ld_bar_penalty + cl_bar_penalty;

%% Make check to verify that sides do not cross
valid_shape = valid_shape_test(experiment_result);
if not(valid_shape)
   sc = inf;
end


% Make additional check to verify that L/D did not enter into unreasonable
% branch