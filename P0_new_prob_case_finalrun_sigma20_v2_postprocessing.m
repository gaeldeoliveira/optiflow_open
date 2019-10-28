addpath src/
addpath user_src/
addpath dev/


%load(['probabilistic/new_results/P0_new_prob_case_finalrun_sigma20_v2_fullcase_20180914T033727.mat'])
load('P20_v2_with_initial_population_for_reference20180914T035817.mat')

% Make ordering index according to first cost function
[~ , i_fval] = sort(results.fval(:,1));
% Reorder second cost funciton accordingly
f1val_ordered = results.fval(i_fval, 1);
f2val_ordered = results.fval(i_fval, 2);
% Reorder experiments of pareto front
pareto_experiments_ordered = pareto_experiments(i_fval, :);

% Reinterpret initial population experiments
val_obj_list = interpret_population_results(GM.CFG, inital_pop_experiments);
f1val_initial = zeros(size(val_obj_list, 1), 1);
f2val_initial = zeros(size(val_obj_list, 1), 1);
for n = 1:size(val_obj_list, 1)
    f1val_initial(n) = val_obj_list{n}(1);
    f2val_initial(n) = val_obj_list{n}(2);
end

% % Plot pareto front
plot(f1val_ordered, f2val_ordered, '.-'); hold on;
plot(f1val_initial, f2val_initial, 'x' );

% Define Chosen designs
index_chosen  = [5 30];
% Plot beautified pareto front
figure(101)
plot(f2val_initial, -100*f1val_initial, 'o' ); hold on; grid on;
plot([0.6; f2val_ordered], -100*[f1val_ordered(1); f1val_ordered], '.-'); 
plot(f2val_ordered(index_chosen), -100*f1val_ordered(index_chosen), 'k*');
axis([0 0.5 120 200])

% for n = 1:size(GM.CM.names_fitted_bounds,2)
%     text(f2val_initial(n), -100*f1val_initial(n), num2str(n));%GM.CM.names_fitted_bounds{n})
% end

text_index = [1 2 5 8 24 25 27];
text_cell  = {'Commercial' ; ...
              'AH-93W257' ; ...
              'FFA-W3211' ; ...
              'FFA-W3241' ; ...
              'Chap. 6 design' ; ...
              'DU91W2250' ; ...
              'DU00W-212' };
for n_text = 1:length(text_index)
    text(f2val_initial(text_index(n_text)) +0.01, -100*f1val_initial(text_index(n_text)) + 1.5, text_cell{n_text})
end
%text(f2val_ordered(index_chosen(1)) - 0.1, -100*f1val_ordered(index_chosen(1)+2), 'High glide design')
%text(f2val_ordered(index_chosen(2)) - 0.1, -100*f1val_ordered(index_chosen(2)+2), 'Smooth ride design')
text('String','Smooth glide design',...
    'Position',[0.17437461843318 158.74057077935 0]);
text('String','High glide design',...
    'Position',[0.288639826267281 188.494850300688 0]);

xlabel('Cl overshoot in [\alpha_{des},\alpha_{des}+2\sigma_\alpha] range of angles of attack')
ylabel('Expected L/D at \alpha_{des}')
title('Feasible airfoil characteristics (Re=12e6, \sigma_\alpha = 2.0deg)')
legend('Existing designs', 'Feasibility frontier (Pareto front)', 'Chosen designs', 'Location', 'SouthWest')

print -depsc probabilistic/export_figures/integrated_design_compromises.eps
savefig('probabilistic/export_figures/integrated_design_compromises.fig')


% % Now select new designs
% n = 0;
% n = n+1; figure(1); h_p = pareto_experiments_ordered{n, 1}.ap.plot() ; pareto_experiments_ordered{n, 2}.ap.plot(h_p, 'r'); title(num2str(n)); figure(2); plot(pareto_experiments_ordered{n, 1}.coordinates.tx, pareto_experiments_ordered{n, 1}.coordinates.tz); grid on; axis equal; figure(3); plot(f2val_ordered(n), -100*f1val_ordered(n), '*')

% Chosen designs
index_chosen  = [5 30];
figure(201);
subplot(221)
plot(pareto_experiments_ordered{index_chosen(1), 1}.ap.alpha_range, pareto_experiments_ordered{index_chosen(1), 1}.ap.cl_alpha(pareto_experiments_ordered{index_chosen(1), 1}.ap.alpha_range)); hold on
plot(pareto_experiments_ordered{index_chosen(1), 2}.ap.alpha_range, pareto_experiments_ordered{index_chosen(1), 2}.ap.cl_alpha(pareto_experiments_ordered{index_chosen(1), 2}.ap.alpha_range)); grid on
ylabel('Cl'); xlabel('\alpha [deg]')
title(['High glide design Cl_{max} = ' num2str(max(pareto_experiments_ordered{index_chosen(1), 1}.ap.cl_alpha(pareto_experiments_ordered{index_chosen(1), 1}.ap.alpha_range)), 2)])
axis([0 20 0 2.5]);
legend('Free transition', 'Forced transition', 'Location', 'SouthWest')

subplot(222)
plot(pareto_experiments_ordered{index_chosen(2), 1}.ap.alpha_range, pareto_experiments_ordered{index_chosen(2), 1}.ap.cl_alpha(pareto_experiments_ordered{index_chosen(2), 1}.ap.alpha_range)); hold on
plot(pareto_experiments_ordered{index_chosen(2), 2}.ap.alpha_range, pareto_experiments_ordered{index_chosen(2), 2}.ap.cl_alpha(pareto_experiments_ordered{index_chosen(2), 2}.ap.alpha_range)); grid on
ylabel('Cl'); xlabel('\alpha [deg]')
axis([0 20 0 2.5]);
title(['Smooth glide design Cl_{max} = ' num2str(max(pareto_experiments_ordered{index_chosen(2), 1}.ap.cl_alpha(pareto_experiments_ordered{index_chosen(2), 1}.ap.alpha_range)), 2)])
legend('Free transition', 'Forced transition', 'Location', 'SouthWest')

subplot(223)
alpha_range = pareto_experiments_ordered{index_chosen(1), 1}.ap.alpha_range;
plot(alpha_range, pareto_experiments_ordered{index_chosen(1), 1}.ap.cl_alpha(alpha_range) ./ pareto_experiments_ordered{index_chosen(1), 1}.ap.cd_alpha(alpha_range)); hold on
alpha_range = pareto_experiments_ordered{index_chosen(1), 2}.ap.alpha_range;
plot(alpha_range, pareto_experiments_ordered{index_chosen(1), 2}.ap.cl_alpha(alpha_range) ./ pareto_experiments_ordered{index_chosen(1), 2}.ap.cd_alpha(alpha_range)); grid on
ylabel('L/D'); xlabel('\alpha [deg]')
title(['High glide design (L/D)^{clean}_{max} = ' , num2str(max(pareto_experiments_ordered{index_chosen(1), 1}.ap.cl_alpha(alpha_range) ./ pareto_experiments_ordered{index_chosen(1), 1}.ap.cd_alpha(alpha_range)), 3)]) 
axis([0 20 0 200]);
legend('Free transition', 'Forced transition', 'Location', 'SouthWest')

subplot(224)
alpha_range = pareto_experiments_ordered{index_chosen(2), 1}.ap.alpha_range;
plot(alpha_range, pareto_experiments_ordered{index_chosen(2), 1}.ap.cl_alpha(alpha_range) ./ pareto_experiments_ordered{index_chosen(2), 1}.ap.cd_alpha(alpha_range)); hold on
alpha_range = pareto_experiments_ordered{index_chosen(2), 2}.ap.alpha_range;
plot(alpha_range, pareto_experiments_ordered{index_chosen(2), 2}.ap.cl_alpha(alpha_range) ./ pareto_experiments_ordered{index_chosen(2), 2}.ap.cd_alpha(alpha_range)); grid on
ylabel('L/D'); xlabel('\alpha [deg]')
title(['Smooth glide design (L/D)^{clean}_{max} = '  , num2str(max(pareto_experiments_ordered{index_chosen(2), 1}.ap.cl_alpha(alpha_range) ./ pareto_experiments_ordered{index_chosen(2), 1}.ap.cd_alpha(alpha_range)), 3)]) 
axis([0 20 0 200]);
legend('Free transition', 'Forced transition', 'Location', 'SouthWest')

print -depsc probabilistic/export_figures/chosen_airfoil_polars.eps
savefig('probabilistic/export_figures/chosen_airfoil_polars.fig')

figure(301)
subplot(211)
plot(pareto_experiments_ordered{index_chosen(1), 1}.coordinates.tx, pareto_experiments_ordered{index_chosen(1), 1}.coordinates.tz); grid on; axis equal; hold on
plot(pareto_experiments_ordered{index_chosen(2), 1}.coordinates.tx, pareto_experiments_ordered{index_chosen(2), 1}.coordinates.tz); 
axis([-0.05 1.05 -0.15 0.2])
legend('High glide design', 'Smooth glide design', 'Location', 'SouthEast')
print -depsc probabilistic/export_figures/chosen_airfoil_shapes.eps
savefig('probabilistic/export_figures/chosen_airfoil_shapes.fig')

% Now save coordinates
high_glide_coordinates = [pareto_experiments_ordered{index_chosen(1), 1}.coordinates.tx, pareto_experiments_ordered{index_chosen(1), 1}.coordinates.tz];
smooth_glide_coordinates = [pareto_experiments_ordered{index_chosen(2), 1}.coordinates.tx, pareto_experiments_ordered{index_chosen(2), 1}.coordinates.tz];

save('probabilistic/high_glide_design.air'  , 'high_glide_coordinates'  , '-ascii', '-double')
save('probabilistic/smooth_glide_design.air', 'smooth_glide_coordinates', '-ascii', '-double')

% % Now plot

lambda = 0.5

n_design_vector = 1:length(f1val_ordered)

f_composite = @(lambda) lambda * f1val_ordered + (1-lambda) * f2val_ordered;


% 
figure(10001)
plot(f1val_ordered              , f2val_ordered              , '.-'); hold on
plot(f1val_ordered(index_chosen), f2val_ordered(index_chosen), 'k*'); grid on
% Plot for different lambdas
figure(10002)
plot(f1val_ordered              , f2val_ordered              , '.-'); hold on
plot(f1val_ordered , f_composite(0.06), '.-'); 
plot(f1val_ordered , f_composite(0.24), '.-');
plot(f1val_ordered , f_composite(0.48), '.-');
plot(f1val_ordered(index_chosen), f2val_ordered(index_chosen), 'k*'); grid on
axis([-1.9, -1.4, -0.72,0.42])


fcomp_1 = f_composite(0.06); [min_fcomp_1, i_min_fcomp_1] = min(fcomp_1(1:37));
plot(f1val_ordered(i_min_fcomp_1)*[1 1], [-1 1], 'o-.', 'Color', [0.4 0.4 0.4])
plot(f1val_ordered(i_min_fcomp_1), min_fcomp_1, 'ko')
fcomp_2 = f_composite(0.24); [min_fcomp_2, i_min_fcomp_2] = min(fcomp_2(1:37));
plot(f1val_ordered(i_min_fcomp_2)*[1 1], [-1 1], 'o-.', 'Color', [0.4 0.4 0.4])
plot(f1val_ordered(i_min_fcomp_2), min_fcomp_2, 'ko')
fcomp_4 = f_composite(0.48); [min_fcomp_4, i_min_fcomp_4] = min(fcomp_4(1:37));
plot(f1val_ordered(i_min_fcomp_4)*[1 1], [-1 1], 'o-.', 'Color', [0.4 0.4 0.4])
plot(f1val_ordered(i_min_fcomp_4), min_fcomp_4, 'ko')
legend('Pareto front : f_2 = Cl overshoot', ...
    '\lambda f_1 + (1-\lambda) * f_2  with \lambda=0.06', ...
    '\lambda f_1 + (1-\lambda) * f_2  with \lambda=0.24', ...
    '\lambda f_1 + (1-\lambda) * f_2  with \lambda=0.48', ...
    'High and smooth glide designs', ...
    'Compromise levels', ...
    'Location', 'SouthEast');
xlabel('f_1 = - (Expected L / Expected D) / 100')
title('Identification of compromise levels')
text(f1val_ordered(i_min_fcomp_1), -0.3, '\lambda = 0.06')
text(f1val_ordered(i_min_fcomp_2), -0.3, '\lambda = 0.24')
text(f1val_ordered(i_min_fcomp_4), -0.3, '\lambda = 0.48')

print -depsc probabilistic/export_figures/compromise_levels.eps
savefig('probabilistic/export_figures/compromise_levels.fig')




