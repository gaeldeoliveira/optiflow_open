classdef sgd
    % SGD is a static class collecting methods used for stochastic gradient
    % descent minimization on machine learning problems.
    %
    % Pre-processing methods include: 
    %`  -> Partition database index into training, test and validation sets
    %   -> Fold index vector into subsets with prescribed lenght or number 
    %   -> Shuffle index vectors or subsets thereof
    %
    % Minimization runtime methods include: 
    %   -> Central and forward gradient estimation 
    %   -> 
    
    properties
        
    end
    
    methods
        function sgd_object = sgd()
            %Dummy constructor
        end
    end
    
    methods(Static)
        function [full_train_subset, monitoring_subset, validation_subset, cell_of_training_folds, number_of_training_folds] = ... 
            shuffle_partition_and_fold_to_length(reference_array, training_fold_length, monitoring_fraction, validation_fraction, rngstate)
            
            % Get lenght of reference array
            total_length      = length(reference_array);
            
            % Shuffled reference array
            shuffled_array    = sgd.shuffle_array(reference_array, rngstate);
            
            % Make indices of validation set
            validation_index  = 1:floor(validation_fraction*total_length);
            % Make indices of monitoring set
            monitoring_index  = floor(validation_fraction*total_length) + (1:floor(monitoring_fraction*total_length));
            % Make indices of training set (unfolded)
            full_train_index  = (floor(validation_fraction*total_length) + floor(monitoring_fraction*total_length) + 1):total_length;
            
            % Compute number of training elements that exceed needs
            n_excess_train    = round(training_fold_length * ((length(full_train_index) / training_fold_length) - floor(length(full_train_index) / training_fold_length)));
            
            % Add excess elements to monitoring set
            monitoring_index(end+(1:n_excess_train)) = full_train_index(1:n_excess_train);
            % Remove these elements from training set
            full_train_index  = full_train_index((n_excess_train+1):end);
            
            % Now make subsets of shuffled array
            validation_subset = shuffled_array(validation_index);
            monitoring_subset = shuffled_array(monitoring_index);
            full_train_subset = shuffled_array(full_train_index);
            
            % And fold training set!
            [cell_of_training_folds, number_of_training_folds]= sgd.fold_to_lenght(full_train_subset, training_fold_length);
            
            % Done! Return!!!
        end        
        
        function shuffled_array = shuffle_array(reference_array, rngstate)
            % Shuffle an array without repetition: random permutation.
            % Applies to arrays of indices or arbitrary type elements.
            
            % Set random number generator state to taste
            rng(rngstate);
            
            % Shuffle it around
            shuffled_array = reference_array(randperm(length(reference_array)));
            
            % Done! Return!!!
        end
        
        function [cell_of_folds , fold_number]= fold_to_lenght(reference_array, fold_length)
            % Splits reference array into folds without repetition. Last
            % fold is longer when division remainder of ordered index
            % length by fold number is non-zero. Fold lenght specified,
            % fold count determined at runtime. 
            
            % Determine fold lenght
            fold_number                = floor(length(reference_array) / fold_length);
            
            % Allocate cell array of folds
            cell_of_folds              = cell(1,fold_number);
            
            % Make first (fold_number-1) first folds
            for n_fold = 1:(fold_number-1)
                cell_of_folds{n_fold}  = reference_array((fold_length*(n_fold-1)+1):(fold_length*n_fold));
            end
            
            % Make last fold
            cell_of_folds{fold_number} = reference_array((fold_length*(fold_number-1)+1):end);
            
            % Done! Return!!!
        end
        
        function [cell_of_folds, fold_length] = fold_number_of_times(reference_array, fold_number)
            % Splits into folds without repetition. Last fold is longer when
            % division remainder of ordered index length by fold number is
            % non-zero. Fold count specified, fold lenght determined at
            % runtime. 
            
            % Determine fold lenght
            fold_length                = floor(length(reference_array) / fold_number);
            
            % Allocate cell array of folds
            cell_of_folds              = cell(1,fold_number);
            
            % Make first (fold_number-1) first folds
            for n_fold = 1:(fold_number-1)
                cell_of_folds{n_fold}  = reference_array((fold_length*(n_fold-1)+1):(fold_length*n_fold));
            end
            
            % Make last fold
            cell_of_folds{fold_number} = reference_array((fold_length*(fold_number-1)+1):end);
            
            % Done! Return!!!
        end        
        
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
                
    end
end

