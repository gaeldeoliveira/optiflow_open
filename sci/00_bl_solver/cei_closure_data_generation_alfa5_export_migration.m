%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Integral Boundary Layer Integrator (ODE)
%           Plasma Development Tool
%
%       August 2014, GNU-GPLv3 or later
%       Gael de Oliveira, Ricardo Pereira
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Energy Interaction Coefficient Closure DataSet Generator
%           following:
%               Modelling the effect of DBD plasma actuators on boundary
%               layer development
%               Internal Report, Gael de Oliveira, Ricardo Pereira
%               September 8, 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
% Clear Stuff
close all; clear all; clc;

% Add some useful paths
addpath 0_closure_relations
addpath 1_definitions
addpath 2_helper_functions

%% Inputs

% % Data Generation Sets
% number of points per dimension
np = 6;
% Set of shape factors (hk)
hk_range = [1 , linspace(1.4, 6.5, np-1)];              % Linear, stop at 6.5 to avoid kink when going above 8 
% Set of Momentum Thickness Reynolds (Re_theta = rt)
rt_range = [400, logspace(2, 5, np-1)*5];               % 400 to 500 to 500.000 (400 has some (mild) imag, but 500 not, so step from there!)
% Set of Momentum Scaled Plasma Force Field Coefficients (t_theta_p)
t_theta_p_range = [0 , logspace(-2, 1, np-1)*3];        % 0 to 0.03 to 30

% % Numerics parameters
N_integration_points = 1000;                            % Integration seems very stable with naive middle point riemann scheme, for 100 points or more (1000 is still a good measure, specially for low shape factors)

% Operation Mode
flag_compute = true;                                    % Compute and Save Dataset to file
flag_load    = false;                                   % Load Dataset from file
dataset_file = 'cei_closure_data_export_migration.mat'; % Dataset filename

%% Object Instanciation
% Construct Swafford Profile Object
SP = swafford_profile();

% Construct Plasma Descriptor Object (used for wy_function only)
PD = plasma_descriptor();

% Construct cei integrand function
cei_integrand = @(yt, b) SP.evaluate_profile(yt) .* PD.wy_function_generic(yt, b) ./ b;

%% Data Allocation
% Allocate case arrays for each study dimension
N_hk        = length(hk_range);
N_rt        = length(rt_range);
N_t_theta_p = length(t_theta_p_range);

% Create gridmesh arrays for each study dimension
% [X,Y,Z] = meshgrid([1 2] ,[4 5], [7 8 9])
 [hk_grid , rt_grid , t_theta_p_grid] = meshgrid(hk_range, rt_range, t_theta_p_range);
%[rt_grid, hk_grid , t_theta_p_grid] = meshgrid(rt_range,hk_range, t_theta_p_range);
% And allocate an array to store computation results
cei_grid = zeros(size(hk_grid));

%% Computation
% Only compute if flag says we should do so!
if (flag_compute)
    
    % Now run through all possible cases! (PARFOR may be a good option here! but race conditions on object definition should be accounted for, so, only parfor the theta loop, and keep the rt and hk loops sequential!)
    % Stick with sequential code for compatiblity (for now!)
    for n_hk = 1:N_hk
        for n_rt = 1:N_rt
            for n_t_theta_p = 1:N_t_theta_p
                % Inner loop, where all computations take place
                % Get current case values
                hk          = hk_grid(n_hk , n_rt, n_t_theta_p);
                rt          = rt_grid(n_hk , n_rt, n_t_theta_p);
                t_theta_p   = t_theta_p_grid(n_hk , n_rt, n_t_theta_p);
                
                % Define Profile
                SP.update_hk_rt_pair(hk, rt);
                
                % Check if we are on a particular case (like hk=1 or 
                if and(not(hk==1), not(t_theta_p == 0))
                    % Define set of Riemann integration points
                    yt_range    = linspace(0, t_theta_p, N_integration_points);                    
                    % Scale Integrand Variable and Perform Naive Riemmann Integration
                    cei_grid(n_hk , n_rt, n_t_theta_p)  = ...
                        riemann_integral(@(yt) cei_integrand(yt, t_theta_p), yt_range);
                elseif hk==1
                    % Analytical Solution! See report (nice trick, I think!)!
                    cei_grid(n_hk , n_rt, n_t_theta_p) = 1;
                elseif t_theta_p==0
                    % Analytical Solution! (Fallback to no plasma case!)
                    cei_grid(n_hk , n_rt, n_t_theta_p) = 0;
                end
                
                
                % Display current results
                disp(['|   hk=' , num2str(hk,'%1.2f') , '   |   Rt=' , num2str(rt, '%04.0f'), ...
                    '   |   t_theta_p=' , num2str(t_theta_p, '%2.2f') , '   |   Cei=' , num2str(cei_grid(n_hk , n_rt, n_t_theta_p), '%2.4f'), '   |']);
            end
        end
    end
    
    % Replace NaNs by zeros (only occur for t_theta_p = 0, bad practice but
    % I am tired and for today it works!)
    % cei_grid(isnan(cei_grid)) = zeros(size(cei_grid(isnan(cei_grid))));

    
    % Save data
    save(dataset_file, 'hk_grid' , 'rt_grid' , 't_theta_p_grid' , 'cei_grid');
    
end

%% Load Data set 
% Only load if flag says we should do so!
if (flag_load)
    % Reload it
    load(dataset_file)
end

%% Export to fortran code (float->symbolic->fortran string->code file)

% % Recast grid into table for human readable fortran code
cei_table = zeros(size(cei_grid, 2), size(cei_grid,1), size(cei_grid, 3));
for i=1:size(cei_table,1)
    for j=1:size(cei_table,2)
        for k=1:size(cei_table,3)
            cei_table(i,j,k) = real(cei_grid(j,i,k));
            if abs(imag(cei_grid(j,i,k))) > 0
                disp('WARNING: Imaginary Components ?(-1)^-3? truncated!')
            end
        end
    end
end

% % Start by making a carriage return with line feed (handy all along)
crlf = char(10); % [char(13),char(10)];

% % Make range vector strings (float->symbolic->fortran string)
%       Convert to symbolic format
HK_RANGE    = sym(hk_range); 
RT_RANGE    = sym(rt_range);
TPPT_RANGE  = sym(t_theta_p_range);
% Generate fortran compatible strings for ranges (concatenate a bit later!)
s1 = fortran(HK_RANGE);
s2 = fortran(RT_RANGE);
s3 = fortran(TPPT_RANGE);

% % Make Range vector Header Strings
% Make header declaration strings for ranges
s1_h = ['      REAL, CLDIMENSION(1,' num2str(length(HK_RANGE  )) ') :: HK_RANGE'];
s2_h = ['      REAL, DIMENSION(1,' num2str(length(RT_RANGE  )) ') :: RT_RANGE'];
s3_h = ['      REAL, DIMENSION(1,' num2str(length(TPPT_RANGE)) ') :: TPPT_RANGE'];

% % Make CEI table strings (float->symbolic->fortran string)
CEI_TABLE = sym(cei_table);  % Convert cei_table to symbolic format
s4 = fortran(CEI_TABLE);     % Export from symbolic format to Fortran code (file output)

% % Make CEI grid header declaration
s4_h = ['      REAL, DIMENSION(' num2str(length(HK_RANGE  )) , ',' ...
                                 num2str(length(RT_RANGE  )) , ',' ...
                                 num2str(length(TPPT_RANGE)) ,     ...
                                 ') :: CEI_TABLE'];
                               
% % Make range size integer strings
s0   = ['      N_HK   = ' , num2str(length(HK_RANGE  )) , crlf, ...
        '      N_RT   = ' , num2str(length(RT_RANGE  )) , crlf, ...
        '      N_TPPT = ' , num2str(length(TPPT_RANGE))];
% % Make range size integer header declaration
s0_h = ['      INTEGER N_HK, N_RT, N_TPPT'];

% Concatenate strings for all data definition!
s_arrays = [s0  , crlf, crlf ...   (size hints)
            s1  , crlf, crlf ...     (hk range)
            s2  , crlf, crlf ...     (rt range)
            s3  , crlf, crlf ...   (tppt range)
            s4  , crlf, crlf]; %     (cei grid)

% Now play the dirty trick, change array contents to replace D's with E's,
% so that constant declaration does not impose double precision!
% (this will only be ok if variable names do not contain D's!)
s_arrays(s_arrays=='D') = 'E';
        
% Concatenate strings for all declaration headers!
s_headers = [s0_h, crlf, ...       (size hints)
             s1_h, crlf, ...         (hk range)
             s2_h, crlf, ...         (rt range)
             s3_h, crlf, ...       (tppt range)
             s4_h, crlf, crlf]; %    (cei grid)


fid = fopen('4_fortan_output/CEI_DATA_migration.f90','w');
if fid>=0
    fprintf(fid, '%s\n', s_headers);
    fprintf(fid, '%s\n', s_arrays);
    fclose(fid)
end


%% Create handle function 
% Use interp3
% Syntax primer: Vq = interp3(X,Y,Z,V,Xq,Yq,Zq) (there may be more
% efficient implementations than this!

cei_function = @(hk, rt, t_theta_p) interp3(hk_grid, rt_grid, t_theta_p_grid, cei_grid, hk, rt, t_theta_p, 'linear');

%% Visualize in Cuts
figure(3); plot(hk_range , cei_function(hk_range, min(rt_range), 1), hk_range , cei_function(hk_range, mean(rt_range), 1) , hk_range , cei_function(hk_range, max(rt_range), 1))
grid on; 
xlabel('hk - Shape Factor') ; ylabel('C_{ei} - Energy Interaction Coefficient');
title('Influence of Momentum Thickness Reynolds @ t_{theta}^p = 1')
legend(['Re_{theta} = ', num2str(min(rt_range)) ], ['Re_{theta} = ' , num2str(mean(rt_range)) ], ['Re_{theta} = ' , num2str(max(rt_range))])

% Print
set(gcf, 'PaperType', 'A5'); orient landscape
print -dpdf ./figures/cei_export_cut_constant_t_theta_p.pdf

figure(4); plot(hk_range , cei_function(hk_range, mean(rt_range), min(t_theta_p_range)), hk_range , cei_function(hk_range, mean(rt_range), mean(t_theta_p_range)) , hk_range , cei_function(hk_range, mean(rt_range), max(t_theta_p_range)))
grid on; 
xlabel('hk - Shape Factor') ; ylabel('C_{ei} - Energy Interaction Coefficient');
title('Influence of Momentum Scaled Plasma Thickness @ Re_{theta}= 1000')
legend(['t_{theta}^p = ', num2str(min(t_theta_p_range))] , ['t_{theta}^p = ' , num2str(mean(t_theta_p_range))] , ['t_{theta}^p =', num2str(max(t_theta_p_range))])

% Print
set(gcf, 'PaperType', 'A5'); orient landscape
print -dpdf ./figures/cei_export_cut_constant_rt.pdf


%% Older Code

% %% Visualize in Slices
% % % Start by plotting a few slices of constant hk, rt, or t_theta_p
% figure(1)
% % slice(hk_grid, rt_grid, t_theta_p_grid, cei_grid, [], min(rt_range), linspace(min(t_theta_p_range) , max(t_theta_p_range), 4))
% % slice(hk_grid, rt_grid, t_theta_p_grid, cei_grid, [], linspace(min(rt_range) , max(rt_range), 3), [])
% slice(hk_grid, rt_grid, t_theta_p_grid, cei_grid, max(hk_range), [mean(rt_range) , max(rt_range)], 0.5 * mean(t_theta_p_range) + 0.5 * min(t_theta_p_range))
% alpha(0.9); view(-58, 24);
% shading interp ; colorbar
% xlabel('hk - Shape Factor') ; ylabel('rt - Reynolds Theta') ; zlabel('t_{theta}^p - Momenutm Scale Thickness')
% title('C_{ei} - Energy Interaction Coefficient -- Slices')
% axis([min(hk_range) , max(hk_range), min(rt_range), max(rt_range) , min(t_theta_p_range), max(t_theta_p_range)])
% 
% % Print
% % set(gcf, 'PaperType', 'A5'); orient landscape
% print -dpng ./figures/cei_slices.png

% %% Visualize in Isosurfaces
% % Lets now attempt some isosurfaces
% figure(2)
% % Make patch on isosurface!
% p1 = patch(isosurface(hk_grid, rt_grid, t_theta_p_grid, cei_grid,0));
% p2 = patch(isosurface(hk_grid, rt_grid, t_theta_p_grid, cei_grid,0.3));
% p3 = patch(isosurface(hk_grid, rt_grid, t_theta_p_grid, cei_grid,0.5));
% % Determine Normals for Lighting
% isonormals(hk_grid, rt_grid, t_theta_p_grid, cei_grid,p1)
% isonormals(hk_grid, rt_grid, t_theta_p_grid, cei_grid,p2)
% isonormals(hk_grid, rt_grid, t_theta_p_grid, cei_grid,p3)
% % Set isosurface color
% set(p1,'FaceColor','blue','EdgeColor','none')
% set(p2,'FaceColor','green','EdgeColor','none')
% set(p3,'FaceColor','red','EdgeColor','none')
% % Set labels, title and axes
% xlabel('hk - Shape Factor') ; ylabel('rt - Reynolds Theta') ; zlabel('t_{theta}^p - Momenutm Scale Thickness')
% title('C_{ei} - Energy Interaction Coefficient -- Isosurfaces')
% axis([min(hk_range) , max(hk_range), min(rt_range), max(rt_range) , min(t_theta_p_range), max(t_theta_p_range)])
% grid on
% % Set view, light and camera
% view(3)
% try 
%     % Must be written in try catch statement to avoid matlab bugs!
%     camlight(0, 0)
% catch
% end
% lighting gouraud
% alpha(0.8)
% legend('C_{ei}=0.0 Isosurface', 'C_{ei}=0.3 Isosurface' , 'C_{ei}=0.5 Isosurface', 'Location', 'NorthEast')
% 
% % Print
% % set(gcf, 'PaperType', 'A5'); orient landscape
% print -dpng ./figures/cei_isosurface.png




% 
% 
% %% Now prepare for data fitting
% 
% 
% [hk_fit, tt_fit] = meshgrid(hk_range, t_theta_p_range);
% 
% rt_fit_500 = 500 * ones(size(hk_fit));
% cei_fit_rt500 = cei_function(hk_fit, rt_fit_500 , tt_fit);
% 
% rt_fit_1000 = 1000 * ones(size(hk_fit));
% cei_fit_rt1000 = cei_function(hk_fit, rt_fit_1000 , tt_fit);
% 
% rt_fit_10000 = 10000 * ones(size(hk_fit));
% cei_fit_rt10000 = cei_function(hk_fit, rt_fit_10000 , tt_fit);
% 
% 
% %% Fit at 500
% 
% 
%    f(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y 
%                     + p12*x*y^2 + p40*x^4 + p31*x^3*y + p22*x^2*y^2
% Coefficients (with 95% confidence bounds):
%        p00 =      0.9983
%        p10 =     -0.8585
%        p01 =      0.1558
%        p20 =      0.2544
%        p11 =    -0.02162
%        p02 =    -0.01172
%        p30 =     -0.0318
%        p21 =   -0.004521
%        p12 =    0.004589
%        p40 =    0.001437
%        p31 =   0.0005442
%        p22 =  -0.0003136
%        
% %% Fit at RT=10000
% 
%        p00 =       1.407
%        p10 =     -0.8978
%        p01 =   -0.007943
%        p20 =      0.1868
%        p11 =     0.05688
%        p02 =    -0.00152
%        p30 =    -0.01613
%        p21 =    -0.01374
%        p12 =   4.121e-05
%        p40 =   0.0004755
%        p31 =    0.000802
%        p22 =   9.712e-05
% 
