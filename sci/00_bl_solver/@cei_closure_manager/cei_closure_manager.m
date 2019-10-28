classdef cei_closure_manager < handle
    %CEI_CLOSURE_MANAGER is a handle class to (re)load
    % Energy Interaction Coefficient closure datasets and interpolate them
    % to create analytical approximations/fits
    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %       Integral Boundary Layer Integrator (ODE)
    %           Plasma Development Tool
    %
    %       August 2014, GNU-GPLv3 or later
    %       Gael de Oliveira, Ricardo Pereira
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %       Energy Interaction Coefficient Closure Dataset Loader
    %           following:
    %               Modelling the effect of DBD plasma actuators on boundary
    %               layer development
    %               Internal Report, Gael de Oliveira, Ricardo Pereira
    %               September 8, 2014
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    properties
        % Dataset Properties
        hk_range                                % Set of shape factors (hk)
        rt_range                                % Set of Momentum Thickness Reynolds (Re_theta = rt)
        t_theta_p_range                         % Set of Momentum Scaled Plasma Force Field Coefficients (t_theta_p)        
        
        % Dataset
        hk_grid                                 % hk dataset
        rt_grid                                 % rt dataset
        t_theta_p_grid                          % t_theta_p dataset
        cei_grid                                % Cei dataset
        
        cei_softening = 1                       % Tunning Factor
               
        % Dataset File
        dataset_file = 'cei_closure_data.mat';  % Dataset filename
        
        % Handles
        PD                                      % Handle to Plasma Descriptor object

    end    
    
    methods
        function CEI = cei_closure_manager(dataset_file)
            % Creator Function: create yourself and load stuff!
            
            % Load dataset
            CEI.load_dataset(dataset_file)
        end
        
        function load_dataset(CEI, dataset_file)
            % Store dataset_file name
            CEI.dataset_file = dataset_file;
            
            % Load data
            dataset = load(CEI.dataset_file, 'hk_grid' , 'rt_grid' , 't_theta_p_grid' , 'cei_grid');
            
            % Store it into the object fields
            CEI.hk_grid         = dataset.hk_grid;
            CEI.rt_grid         = dataset.rt_grid;
            CEI.t_theta_p_grid  = dataset.t_theta_p_grid;
            CEI.cei_grid        = dataset.cei_grid;
            
            % Proceed further into identification of dataset ranges
            CEI.hk_range             = CEI.hk_grid(1,:,1);
            CEI.hk_range             = CEI.hk_range(:)';
            CEI.rt_range             = CEI.rt_grid(:,1,1);
            CEI.rt_range             = CEI.rt_range(:)';
            
            CEI.t_theta_p_range      = zeros(1, size(CEI.t_theta_p_grid(1,1,:), 3));
            CEI.t_theta_p_range(:)   = CEI.t_theta_p_grid(1,1,:);
        end
        
        function cei = cei_function(CEI , hk, rt, t_theta_p)
           % Returns cei value, as interpolated from grid!
           % cei_function = @(hk, rt, t_theta_p) interp3(hk_grid, rt_grid, t_theta_p_grid, cei_grid, hk, rt, t_theta_p);
           
           % Fudge unidimensionally when needed!
           hk = min(hk, max(CEI.hk_range));
           hk = max(hk, min(CEI.hk_range));              
           
           rt = min(rt, max(CEI.rt_range));
           rt = max(rt, min(CEI.rt_range));
           
           t_theta_p = min(t_theta_p, max(CEI.t_theta_p_range));
           t_theta_p = max(t_theta_p, min(CEI.t_theta_p_range));

           % Interpolate against data
%            disp(['|   hk=' , num2str(hk,'%1.2f') , '   |   Rt=' , num2str(rt, '%04.0f'), ...
%                '   |   t_theta_p=' , num2str(t_theta_p, '%2.2f') , '   |']);
           
           cei = interp3(CEI.hk_grid, CEI.rt_grid, CEI.t_theta_p_grid, CEI.cei_grid, hk, rt, t_theta_p) * CEI.cei_softening;
        end
        
        
    end
    
end



% %% Initialization
% % Clear Stuff
% close all; clear all; clc;
% 
% % Add some useful paths
% addpath 0_closure_relations
% addpath 1_definitions
% addpath 2_helper_functions
% 
% %% Inputs
% 
% % % Data Generation Sets
% % Set of shape factors (hk)
% hk_range = linspace(1.4, 8, 40);
% % Set of Momentum Thickness Reynolds (Re_theta = rt)
% rt_range = linspace(500, 5000, 40);
% % Set of Momentum Scaled Plasma Force Field Coefficients (t_theta_p)
% t_theta_p_range = linspace(0, 10, 40);
% 
% % % Numerics parameters
% N_integration_points = 1000;            % Integration seems very stable with naive middle point riemann scheme, for 100 points or more (1000 is still a good measure, specially for low shape factors)
% 
% % Operation Mode
% flag_compute = false;                   % Compute and Save Dataset to file
% flag_load    = true;                    % Load Dataset from file
% dataset_file = 'cei_closure_data.mat';  % Dataset filename
% 
% %% Object Instanciation
% % Construct Swafford Profile Object
% SP = swafford_profile();
% 
% % Construct Plasma Descriptor Object (used for wy_function only)
% PD = plasma_descriptor();
% 
% % Construct cei integrand function
% cei_integrand = @(yt, b) SP.evaluate_profile(yt) .* PD.wy_function_generic(yt, b) ./ b;
% 
% %% Data Allocation
% % Allocate case arrays for each study dimension
% N_hk        = length(hk_range);
% N_rt        = length(rt_range);
% N_t_theta_p = length(t_theta_p_range);
% 
% % Create gridmesh arrays for each study dimension
% % [X,Y,Z] = meshgrid([1 2] ,[4 5], [7 8 9])
% [hk_grid , rt_grid , t_theta_p_grid] = meshgrid(hk_range, rt_range, t_theta_p_range);
% % And allocate an array to store computation results
% cei_grid = zeros(size(hk_grid));
% 
% %% Computation
% % Only compute if flag says we should do so!
% if (flag_compute)
%     
%     % Now run through all possible cases! (PARFOR may be a good option here! but race conditions on object definition should be accounted for, so, only parfor the theta loop, and keep the rt and hk loops sequential!)
%     % Stick with sequential code for compatiblity (for now!)
%     for n_hk = 1:N_hk
%         for n_rt = 1:N_rt
%             for n_t_theta_p = 1:N_t_theta_p
%                 % Inner loop, where all computations take place
%                 % Get current case values
%                 hk          = hk_grid(n_hk , n_rt, n_t_theta_p);
%                 rt          = rt_grid(n_hk , n_rt, n_t_theta_p);
%                 t_theta_p   = t_theta_p_grid(n_hk , n_rt, n_t_theta_p);
%                 
%                 % Define Profile
%                 SP.update_hk_rt_pair(hk, rt);
%                 
%                 % Define set of Riemann integration points
%                 yt_range    = linspace(0, t_theta_p, N_integration_points);
%                 
%                 % Scale Integrand Variable and Perform Naive Riemmann Integration
%                 cei_grid(n_hk , n_rt, n_t_theta_p)  = ...
%                     riemann_integral(@(yt) cei_integrand(yt, t_theta_p), yt_range);
%                 
%                 % Verify that there are no imaginary components, and, if
%                 % there are, declare case as invalid
%                 
%                 
%                 % Display current results
%                 disp(['|   hk=' , num2str(hk,'%1.2f') , '   |   Rt=' , num2str(rt, '%04.0f'), ...
%                     '   |   t_theta_p=' , num2str(t_theta_p, '%2.2f') , '   |   Cei=' , num2str(cei_grid(n_hk , n_rt, n_t_theta_p), '%2.4f'), '   |']);
%             end
%         end
%     end
%     
%     % Replace NaNs by zeros (only occur for t_theta_p = 0, bad practice but
%     % I am tired and for today it works!)
%     cei_grid(isnan(cei_grid)) = zeros(size(cei_grid(isnan(cei_grid))));
% 
%     
%     % Save data
%     save(dataset_file, 'hk_grid' , 'rt_grid' , 't_theta_p_grid' , 'cei_grid');
%     
% end
% 
% %% Load Data set 
% % Only load if flag says we should do so!
% if (flag_load)
%     % Reload it
%     load(dataset_file)
% end
% 
% %% Create handle function 
% % Use interp3
% % Syntax primer: Vq = interp3(X,Y,Z,V,Xq,Yq,Zq) (there may be more
% % efficient implementations than this!
% 
% cei_function = @(hk, rt, t_theta_p) interp3(hk_grid, rt_grid, t_theta_p_grid, cei_grid, hk, rt, t_theta_p);
% 
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
% 
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
% %% Visualize in Cuts
% figure(3); plot(hk_range , cei_function(hk_range, 500, 1), hk_range , cei_function(hk_range, 1000, 1) , hk_range , cei_function(hk_range, 5000, 1))
% grid on; 
% xlabel('hk - Shape Factor') ; ylabel('C_{ei} - Energy Interaction Coefficient');
% title('Influence of Momentum Thickness Reynolds @ t_{theta}^p = 1')
% legend('Re_{theta} =  500' , 'Re_{theta} = 1000' , 'Re_{theta} = 5000')
% 
% % Print
% set(gcf, 'PaperType', 'A5'); orient landscape
% print -dpdf ./figures/cei_cut_constant_t_theta_p.pdf
% 
% figure(4); plot(hk_range , cei_function(hk_range, 1000, 1), hk_range , cei_function(hk_range, 1000, 3.5) , hk_range , cei_function(hk_range, 1000, 6))
% grid on; 
% xlabel('hk - Shape Factor') ; ylabel('C_{ei} - Energy Interaction Coefficient');
% title('Influence of Momentum Scaled Plasma Thickness @ Re_{theta}= 1000')
% legend('t_{theta}^p =  1' , 't_{theta}^p = 3.5' , 't_{theta}^p = 6')
% 
% % Print
% set(gcf, 'PaperType', 'A5'); orient landscape
% print -dpdf ./figures/cei_cut_constant_rt.pdf

