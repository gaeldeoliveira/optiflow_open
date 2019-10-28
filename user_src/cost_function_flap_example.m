function sc = cost_function_flap_example(~ , experiments_results)
%function sc = cost_function_flap_example(parameters , experiments_results)
%   Is a simple example on how to access data in the flap case (3
%   simulations!)


%% Parameter Parsing
    % No parameters for this case. Replace the tilda (~) with an entry for
    % a parameter structure and set it in the cost function definition of
    % your case if you want to parametrize your cost function to make it more generic!

%% Input extraction
% Extract results structure and aerodynamic polar for each of the three
% flap positions!

% First extract undeflected case
result_undeflected  = experiments_results{1};
ap_undeflected = result_undeflected.ap;

% Second downward deflection case
result_downward  = experiments_results{1};
ap_downward = result_downward.ap;

% Third upward deflection case
result_upward  = experiments_results{1};
ap_upward = result_upward.ap;

%% Process
% First check that polars are valid to avoid unpredictable behaviour and
% destruction of optimization results
if and(and(strcmp(class(ap_undeflected) , 'aerodynamic_polar') , ...
        strcmp(class(ap_downward) , 'aerodynamic_polar')) , ...
        strcmp(class(ap_upward) , 'aerodynamic_polar')); %#ok<STISA>
    % If polar is valid
    sc = ap_upward.cl_max_local_limited - ap_downward.cl_max_local_limited;
else
    % If polar is invalid!
    sc = 0;
end