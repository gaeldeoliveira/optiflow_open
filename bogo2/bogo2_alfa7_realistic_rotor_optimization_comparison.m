%

% Add path to load classes
addpath src/

% 
date_str = 'DO_20180915T200345';
stat_FFA_free          = load([date_str '/stat_FFA_free.mat'   ]);
prob_FFA_free          = load([date_str '/prob_FFA_free.mat'   ]);

plot(stat_FFA_free.BR.mu_vector, stat_FFA_free.BR.a_axi); hold on
plot(prob_FFA_free.BR.mu_vector, prob_FFA_free.BR.a_axi); grid on