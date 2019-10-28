function [f_val, central_grad_val, forward_grad_val] = fval_with_central_grad(fun_handle, x, step_size)
% [f_val, central_grad_val, forward_grad_val] = fval_with_central_grad(fun_handle, x, step_size)
% Computes the function fun_handle at point x and determines its gradient
% to the components of x with central (and forward) finite differences. The
% size of the step in each direction (forward and backward) is provided as
% argument and outputs are given as:   
%
%   fval             - scalar
%   central_grad_val - array with same size as x
%   forward_grad_val - array with same size as x
%
% For now we assume that outputs of fun_handle are scalar but this not stay
% forcefully so in the future. 
%
% This function could easily be made to handle non-scalar outputs and/or
% enable higher paralelism. This is not a critical issue at the moment
% though - we already have a finer grained paralelism going on - but it may
% become relevant once we really start operating on distributed memory
% architectures.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % First compute value of function at center point (x)
    f_val = fun_handle(x);
    
    % Determine number of dimensions in input
    N_dim_inputs = max(size(x));
    
    % Allocate outputs
    central_grad_val = zeros(size(x));
    forward_grad_val = zeros(size(x));
    
    % Compute gradient
    for n_dim_input = 1:N_dim_inputs
        % (re)Initialize temporary position values
        x_forward                     = x;
        x_backward                    = x;
        % Make forward and backward step positions for current study dimension
        x_forward( n_dim_input)       = x(n_dim_input) + step_size;
        x_backward(n_dim_input)       = x(n_dim_input) - step_size;
        % Compute  function values at forward and backward steps for current study dimension
        f_val_forward                 = fun_handle(x_forward );
        f_val_backward                = fun_handle(x_backward);
        % Now compute current gradient component(s)
        central_grad_val(n_dim_input) = ( f_val_forward - f_val_backward ) / (2*step_size);
        forward_grad_val(n_dim_input) = ( f_val_forward - f_val          ) /    step_size ;
    end
    % Return!
    
end