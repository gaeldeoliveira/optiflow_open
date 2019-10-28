%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Similarity Profile Solver (ODE)
%           Plasma Development Tool
%
%       October 2014, GNU-GPLv3 or later
%       Gael de Oliveira, Ricardo Pereira
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all


%% Inputs
% Integration algorithm parameters
bvp_solver  = @bvp4c;                           % bvp4c and bvp5c are the two common choices
N_steps     = 2000;                             % Set number of steps (only applicable for fixed_step mode
% fixed_step  = true;                             % Fixed or variable lenght ODE intgration steps (always true for now, choice will be implemented later!)
% Define Boundary Condition Points (range of integration, eta_B ~= infty)
eta_A       = 0;
eta_B       = 20;


% Define Plasma Strenght
lambda_range = [-.4 -.2 -.1 0 .1 .2 .4];               
t_p          = 2.0;


%% Solve Fitini over Cell
% Allocate solution storage cell array
outputs_cell = cell(length(lambda_range), 1);

% Loop over cases 
for n_lambda = 1:size(outputs_cell, 1)
    % (re)Generate Options Structure
    % Integration algorithm parameters
    options.bvp_solver  = bvp_solver;                       % bvp4c and bvp5c are the two common choices
    options.N_steps     = N_steps;                          % Set number of steps (we are only working in fixed_step mode for now
    
    % Define Plasma Strenght
    options.lambda      = lambda_range(n_lambda);             % Plasma Equilibrium strenght parameter
    options.t_p         = t_p;                              % Plasma Scaled Thickness (to delta)
    
    % Define Boundary Condition Points
    options.eta_A = eta_A;
    options.eta_B = eta_B;
    
    % Solve Case!
    outputs_cell{n_lambda} = solver_fitini(options);
end

%% Display Solution
disp('|  lambda   |    t_p    | 2*f_prime |    dstr   |   theta   |    h12    |  h_str  |');
for n_lambda = 1:size(outputs_cell, 1)
    outputs = outputs_cell{n_lambda};
    disp([num2str([outputs.options.lambda outputs.options.t_p outputs.cf_coef_x , outputs.dstr , outputs.theta , outputs.h  , outputs.h_str], '|  %+6.4f  ') , '|']);
end

%% Plot Speed Profiles
% Create figure
figure(1)
hold on ; grid on
% Prepare Colors
color_list = cool(size(outputs_cell, 1));
color_labels = cell(1,size(outputs_cell, 1));

for n_lambda = 1:size(outputs_cell, 1)
    outputs = outputs_cell{n_lambda};
    % Retrieve Color Label
    color_labels{n_lambda} = num2str(outputs.options.lambda);
    % Plot
    plot(outputs.y_sol(2,:) , outputs.eta_vector , 'Color', color_list(n_lambda,:));
end

% Set colormap and colobar with labels
colormap(color_list)
lcolorbar(color_labels,'fontweight','bold');

% Set axis to usual values
axis([-0.05 1.05 0 8])
% Set axis labels
xlabel('u = U/U_e = f_{(\eta)} - Tangential Speed')
ylabel('\eta  - Normal Coordinate')
% Set title
title(['Velocity Profiles  --  t_p = '  , num2str(t_p) , '  --  lambda = colorcoded'])

% Now export figure
set(gcf , 'PaperType', 'A5');
orient landscape;
print ('-dpdf', ['figures/', 't_p' , num2str(t_p, '%+3.2f') , '_velocity_profile.pdf'])

%% Plot Shear Profiles

% Create figure
figure(2)
hold on ; grid on
% Prepare Colors
color_list = cool(size(outputs_cell, 1));
color_labels = cell(1,size(outputs_cell, 1));

for n_lambda = 1:size(outputs_cell, 1)
    outputs = outputs_cell{n_lambda};
    % Retrieve Color Label
    color_labels{n_lambda} = num2str(outputs.options.lambda);
    % Plot
    plot(outputs.y_sol(3,:) , outputs.eta_vector , 'Color', color_list(n_lambda,:));
end

% Set colormap and colobar with labels
colormap(color_list)
lcolorbar(color_labels,'fontweight','bold');

% Set axis to usual values
axis([-0.05 1.05 0 8])
% Set axis labels
xlabel('du/d\eta = f^{,}_{(\eta)} - Scaled Shear Rate')
ylabel('\eta  - Normal Coordinate')
% Set title
title(['Shear Profiles  --  t_p = '  , num2str(t_p) , '  --  lambda = colorcoded'])

% Now export figure
set(gcf , 'PaperType', 'A5');
orient landscape;
print ('-dpdf', ['figures/', 't_p' , num2str(t_p, '%+3.2f') , '_shear_profile.pdf'])





% 
% plot(outputs.eta_vector , outputs.y_sol(1,:) , 'Color', color_list(1,:));
% plot(outputs.eta_vector , outputs.y_sol(2,:) , 'Color', color_list(2,:));
% plot(outputs.eta_vector , outputs.y_sol(3,:) , 'Color', color_list(3,:));
% 
% labels = {num2str(1),num2str(2),num2str(3)};
% lcolorbar(labels,'fontweight','bold');
% 
% grid on



% % Auxiliary Functions (the mat4bvp example script uses nested functions,
% but I am afraid to make variable scope confusions with nested functions,
% so we will stay like this for now!)  

% % Some Results

% lambda = 0
%  2*f_prime |  dstr    |   theta   |    h12   |    h_str    
% 0.66412      1.7208      0.6641      2.5911      1.5726
% not hardcoded anymore!
% 0.66412      1.7208      0.6641      2.5911      1.5726

% lambda = -0.6
% b = 0.5 (t_p, but hardcoded!)
%  2*f_prime |  dstr    |   theta   |    h12   |    h_str    
% -0.33559      3.0869      1.0609      2.9097      1.5738
% not hardcoded anymore!
% -0.33559      3.0869      1.0609      2.9097      1.5738
% t_p = 1.0  (not hardcoded!)
% -0.3377      3.4338      1.0555      3.2533      1.5777

% lambda = 0.6
% b = 0.5 (t_p, but hardcoded!)
%  2*f_prime |  dstr    |   theta   |    h12   |    h_str    
% 4.4369     0.47507      0.2299      2.0664      1.6313
% not hardcoded anymore!
% 4.4369     0.47507      0.2299      2.0664      1.6313

% 
% folder = '/home/gael.deoliveira/Documents/MATLAB/actiairfoil_0.95/projects/new_architecture_nrel/results/airfoils/Polar Plotter/';
% file_list = {'plr.nocq.5e5.nosuc.rough' , 'plr.nocq.5e5.suc0.005.53t69.roug' , 'plr.nocq.5e5.suc0.007.53t69.roug' , 'plr.nocq.5e5.suc0.008.53t69.roug'};
% 
% 
% ap_list = cell(size(file_list));
% legend_list = cell(size(file_list));
% color_list = cool(length(file_list));
% 
% for n_file = 1:length(file_list)
%     % Load results
%     raw_data = read_polar_file([folder file_list{n_file}], 'rfoilsuc');
% 
%     % Make aerodynamic polar object
%     ap = aerodynamic_polar(raw_data);
%        
%     if n_file == 1
%         ap_axes = ap.plot;
%     else
%         ap.plot(ap_axes, color_list(n_file,:));
%     end
%     
%     legend_list{n_file} = [file_list{n_file} , ' L/D_{max}=' , num2str(max(ap.raw_data.cl ./ ap.raw_data.cd), 3)];
% end
% 
% legend(legend_list , 'Location'  , 'NorthOutside')
% 
% set(gcf , 'PaperType', 'A5');
% orient landscape;
% print ('-dpdf', [folder , 'plr.nocq.5e5.rough.pdf'])
% 
% 
% 
% if ischar(color)
%     plot(ap.cd_alpha([]).*10^4 , ap.cl_alpha([]), color);
% else
%     plot(ap.cd_alpha([]).*10^4 , ap.cl_alpha([]), '.-' , 'Color', color);
% end
% 
% colorbar
% colormap
% caxis
% 
% figure; colormap(jet(5))
% labels = {'apples','oranges','grapes','peaches','melons'};
% lcolorbar(labels,'fontweight','bold');
% 
% %% Plot
% figure(1)
% plot(outputs.eta_vector , outputs.y_sol(1,:)  , ...
%      outputs.eta_vector , outputs.y_sol(2,:)  , ...
%      outputs.eta_vector , outputs.y_sol(3,:) );
% 
% hold on
%  
% plot(outputs.eta_vector , outputs.y_0(1,:)  , '--' , ...
%      outputs.eta_vector , outputs.y_0(2,:)  , '--' ,  ...
%      outputs.eta_vector , outputs.y_0(3,:)  , '--');
%  
%  
% legend('g - blasius function'       , ...
%        'dg - velocity profile'      , ...
%        'd2g - shear stress profile' , ...
%        'Initial Guesses');
% 
% grid on
% axis([0 5 -0.2 1.2])    % Manual Axis setting seems fairly universal!
