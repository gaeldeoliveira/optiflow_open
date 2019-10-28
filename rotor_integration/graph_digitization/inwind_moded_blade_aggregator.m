% Load inwind rotor
load('inwind_RWT.mat')

% First copy original blade
inwind_moded_RWT = inwind_RWT;

% Define function for scaling thickness distribution
tc_map_fun = @(tc) 21 + (tc - 24.1423) * (100 - 21) / (100 - 24.1423);

% And apply scaling
inwind_moded_RWT.val_thickness = tc_map_fun(inwind_RWT.val_thickness);

% Save to file
save('inwind_moded_RWT.mat', 'inwind_moded_RWT');

% Then show result
plot(inwind_moded_RWT.r_thickness, inwind_RWT.val_thickness); hold on
plot(inwind_moded_RWT.r_thickness, inwind_moded_RWT.val_thickness); grid on
legend('original', 'scaled')






