
% Split database into folds without repetition. Last fold is longer when
% division remainder of ordered index length by fold number is non-zero.
% Fold count specified, fold lenght determined at runtime.

% % Define fold generator inputs
% Ordered index of database entries
k_ordered   = 1:length(EDB.EC_cell);
% Number of folds
fold_number = 3;


% % Make shuffled index of database entries
% Set random number generator state to taste
rng(rngstate)
% Shuffle it around
k_shuffled = k_ordered(randperm(length(k_ordered)));

% % Now split shuffled index into several folds 
% Determine fold lenght
fold_length = floor(length(k_shuffled) / fold_number);
% Allocate cell array of folds
k_cell      = cell(1,fold_number);
% Make first (fold_number-1) first folds
for n_fold = 1:(fold_number-1)
    k_cell{n_fold} = k_shuffled((fold_length*(n_fold-1)+1):(fold_length*n_fold));
end
% Make last fold
k_cell{fold_number} = k_shuffled((fold_length*(fold_number-1)+1):end);

% % Done! Return!



% Split database into folds without repetition. Last fold is longer when
% division remainder of ordered index length by fold number is non-zero.
% Fold lenght specified, fold count determined at runtime.

% % Define fold generator inputs
% Ordered index of database entries
k_ordered   = 1:length(EDB.EC_cell);
% Number of folds
fold_length = 2;


% % Make shuffled index of database entries
% Set random number generator state to taste
rng(rngstate)
% Shuffle it around
k_shuffled = k_ordered(randperm(length(k_ordered)));

% % Now split shuffled index into several folds 
% Determine fold lenght
fold_number = floor(length(k_shuffled) / fold_length);
% Allocate cell array of folds
k_cell      = cell(1,fold_number);
% Make first (fold_number-1) first folds
for n_fold = 1:(fold_number-1)
    k_cell{n_fold} = k_shuffled((fold_length*(n_fold-1)+1):(fold_length*n_fold));
end
% Make last fold
k_cell{fold_number} = k_shuffled((fold_length*(fold_number-1)+1):end);


% % Make a partition
% Partition inputs
reference_array = 1:36;
training_fold_length = 3;
monitoring_fraction  = 0.20;
validation_fraction  = 0.40;
% rngstate           


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
[cell_of_training_folds, training_fold_number]= sgd.fold_to_lenght(full_train_subset, training_fold_length);


% % Now repeat this with function
% Partition inputs
reference_array = 1:36;
training_fold_length = 3;
monitoring_fraction  = 0.20;
validation_fraction  = 0.40;
% Partition function
[full_train_subset, monitoring_subset, validation_subset, cell_of_training_folds, number_of_training_folds] = sgd.shuffle_partition_and_fold_to_length(reference_array, training_fold_length, monitoring_fraction, validation_fraction, rngstate);






