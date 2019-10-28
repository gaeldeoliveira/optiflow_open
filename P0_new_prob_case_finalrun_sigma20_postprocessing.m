% Processing of(hopefully, really) last round of probabilistic cases (P0
% series)

% Load case results
load probabilistic/new_results/P0_new_prob_case_finalrun_sigma20_fullcase_20180913T215409.mat

% Sort pareto front


% View pareto front designs
n_design = 0;
n_design = n_design + 1; 
figure(1); 
h_polar = pareto_experiments{n_design , 1}.ap.plot(); 
          pareto_experiments{n_design , 2}.ap.plot(h_polar, 'r'); 
          title(['Number: ' , num2str(n_design)])
figure(2); plot(pareto_experiments{n_design , 1}.coordinates.tx, pareto_experiments{n_design , 1}.coordinates.tz); axis equal; grid on
          title(['fval: ' , num2str(results.fval(n_design, :))])
% 
Design 14 - Sharp stall: high performance design