function outputs = solver_fitini(options)
% outputs = solver_fitini(options)
%
% Solves Fitini, la gentille equation!
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Similarity Profile Solver (ODE)
%           Plasma Development Tool
%
%       October 2014, GNU-GPLv3 or later
%       Gael de Oliveira, Ricardo Pereira
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Retrieve option structure fields
% Integration algorithm parameters
bvp_solver  = options.bvp_solver;                       % bvp4c and bvp5c are the two common choices
N_steps     = options.N_steps;                          % Set number of steps (we are only working in fixed_step mode for now

% Define Plasma Strenght
lambda      = options.lambda;                           % Plasma Equilibrium strenght parameter
t_p         = options.t_p;                              % Plasma Scaled Thickness (to delta)

% Define Boundary Condition Points
eta_A = options.eta_A;
eta_B = options.eta_B;



% % Solution Process

% Define Integration rage
eta_vector  = linspace(eta_A, eta_B, N_steps);

% Generate initial conditions structure
solinit = bvpinit(eta_vector,@odeinit_blasius);

% Solve (and pass lambda to fitini, la gentille équation!
sol = bvp_solver(@(eta, y) odefun_fitini(eta, y, lambda, t_p), @resBC_blasius , solinit);

% Now retrieve solution
y_sol = deval(sol,eta_vector);
y_0   = solinit.y;

% % Fill Output Structure with Useful Things\
outputs.options     = options;
outputs.sol         = sol;
outputs.eta_vector  = eta_vector;
outputs.y_sol       = y_sol;
outputs.y_0         = y_0;
% Generate Speed Handle Function to integrate boundary layer parameters
u            = @(eta) deval(sol, eta, 2);
outputs.u    = u;
dudy         = @(eta) deval(sol, eta, 3);
outputs.dudy = dudy;

% Generate Boundary Layer Parameters
outputs.dstr         = integral( @(eta)          1 - u(eta)     , eta_A, eta_B);
outputs.theta        = integral( @(eta) u(eta).*(1 - u(eta))    , eta_A, eta_B);
outputs.theta_str    = integral( @(eta) u(eta).*(1 - u(eta).^2) , eta_A, eta_B);
outputs.h            = outputs.dstr ./ outputs.theta;
outputs.h_str        = outputs.theta_str ./ outputs.theta;
outputs.cf_coef_x    = sol.y(3,1) * 2;
outputs.cf_coef_t    = sol.y(3,1) * 2 * outputs.theta;

end




