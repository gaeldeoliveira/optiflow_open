classdef vortex_descriptor
    %VORTEX_DESCRIPTOR This handle class stores the deescription of a
    %replicated vortex system and provides methods for evaluating induced
    %velocity fields and ensuring data consistency. %
    %   
    %   Typical Notation
    %       y_v         % Height of vortex pair
    %       z_v         % Distance of vortex to center of vortex pair
    %       gamma_v     % Strenght (circulation) of each individual vortex 
    %       sigma_v     % Core radiua (when applicable)
    
    
    properties
        S                       % Distance between vortex pair centers
        n_cells            =5   % Number of replication cells
        
        induction_function =2   % Type of induction function
                                %   =1  Singular Vortex (dirac forcing/green function)
                                %   =2  Rankine Vortex  (solid core)
                                %   =3  Lamb Vortex     (gaussian core)
    end
    
    methods
        % Constructor Method
        function VD = vortex_descriptor(S, n_cells)
            VD.S = S;
            VD.n_cells = n_cells;
        end
        
        function [w_v, v_v] = induction_vortex_center(VD, z_v, y_v, gamma_v, sigma_v)
            % Function for computing induction of replicated
            % vortex system on position of source (original) vortex
            
            % Set target in the same place as primary source!
            z = z_v;
            y = y_v;
            
            % Compute induction of replicated vortex system on primary
            % source!
            [w_v, v_v] = VD.induction_replicated_duplicated_vortex_system(z , y, z_v, y_v, gamma_v, sigma_v);
        end
        
        function [w_mesh, v_mesh, mag_mesh] = induction_on_mesh(VD, CM, z_v, y_v, gamma_v, sigma_v)
            % Get's velocity field induced by the vortex system (VD) on a mesh of
            % points stored in a crossflow_mesh (CM) object. 
            % Primary vortex descriptors are supplied as arguments: 
            % center (z_v,y_v), strenght (gamma_v) and spreading (sigma_v).
            
            % Extract Mesh Points from Mesh Object
            z_mesh = CM.z_mesh;
            y_mesh = CM.y_mesh;
            
            % Compute convection field
            [w_mesh, v_mesh] = ...
                VD.induction_replicated_duplicated_vortex_system( ...
                                z_mesh, y_mesh, z_v, y_v, gamma_v, sigma_v);
            
            % Compute magnitude of convection field (used for CFL estimation)
            mag_mesh = sqrt(w_mesh.^2 + v_mesh.^2);
    
            % % Then filter out numerical imperfections!
            % Put v to strict 0 (not 1e-16) on first line (y=0, wall)
            v_mesh(1,:) = 0;
            % Put v to strict 0 (not 1e-4) on first and last column (x=+- nS, cell edges)
            w_mesh(:,1) = 0;
            w_mesh(:,end) = 0;
            
            % We would not need to do this if the replication of the vortex 
            % cell system were infinite. Symmetry would be perfect in that 
            % case, but here have some assymetry on the side "central" 
            % symmetry lines! 
            %
            % We have u_mesh in the range of 1e-4, which looks a bit too 
            % big ! What should be noted is that it depends the on number of
            % replications: for example we get u_mesh in the 1e-4 ballpark
            % for n_cells=10 but 1e-7 for n_cells=100!
            % (global symmetry is excellent though u_mesh(:,1) + u_mesh(:,end)=1e-16)
            
            % Restore periodicity explicitly on v_mesh as well!
            v_mesh(:,1  )     = v_mesh(:,CM.j_center);
            v_mesh(:,end)     = v_mesh(:,CM.j_center);
            
        end
        
        function [w, v] = induction_replicated_duplicated_vortex_system(VD, z , y, z_v, y_v, gamma_v, sigma_v)
            % General function for computing induction of replicated
            % vortex system on a list of points (z,y) (can be nD-arrays)
            
            % Replicate Vortices over Cells  (Construct a List)
            [z_v_list, y_v_list, gamma_v_list, sigma_v_list] = ...
                VD.replicate_vortices_sideways(z_v, y_v, gamma_v, sigma_v);
            
            % Duplicate into lower halfplane (Extend the List)
            [z_v_list, y_v_list, gamma_v_list, sigma_v_list] = ...
                vortex_descriptor.duplicate_vortices_to_lower_halfplane(z_v_list, y_v_list, gamma_v_list, sigma_v_list);
            
            % Evaluate Induction of Vortices in List over points z,y (can be n-D arrays!)
            [w, v] = VD.induction_vortex_list(z, y, z_v_list, y_v_list, gamma_v_list, sigma_v_list);
            
        end
        
        function [w, v] = induction_vortex_list(VD,z, y, z_v_list, y_v_list, gamma_v_list, sigma_v_list)
            % x_v, y_v, gamma_v, sigma_v must all be 1d vectors with the same lenght
            % x,y must have same size and number of dimensions, which are arbitrary
            % (except if one of them is scalar, I think? but don't rely on that!)
            
            % Find lenght of list
            N_v = length(z_v_list);
            
            % Now allocate space to gather contributions and return
            w = zeros(size(z));
            v = zeros(size(z));
            
            for n_v = 1:N_v
                % Compute induction of n^th vortex in the list
                [w_n, v_n] = VD.induction_single_vortex(z, y, ...
                    z_v_list(n_v), y_v_list(n_v), ...
                    gamma_v_list(n_v), sigma_v_list(n_v));
                % Add it up!
                w = w + w_n;
                v = v + v_n;
            end
            
            % Done! Return!
        end

        function [w, v] = induction_single_vortex(VD, z, y, z_v, y_v, gamma_v, sigma_v)
            % This function is just a standard interface to call a
            % induction function consistently. Induction functions are
            % static!
            
            if    VD.induction_function == 1
                % Singular (dirac green function) vortex
                [w, v] = vortex_descriptor.induction_singular_vortex(z, y, z_v, y_v, gamma_v);
            elseif VD.induction_function == 2
                % Rankine vortex (solid core)
                [w, v] = vortex_descriptor.induction_rankine_vortex(z, y, z_v, y_v, gamma_v, sigma_v);
            elseif VD.induction_function == 3
                % Lamb vortex (gaussian core)
                [w, v] = vortex_descriptor.induction_lamb_vortex(z, y, z_v, y_v, gamma_v, sigma_v);
            elseif VD.induction_function == 4
                % Lamb vortex (gaussian core)
                [w, v] = vortex_descriptor.induction_lamb_vortex_debugging(z, y, z_v, y_v, gamma_v, sigma_v);
            end
        end
        
        function [z_v_list, y_v_list, gamma_v_list, sigma_v_list] = ...
                replicate_vortices_sideways(VD, z_v, y_v, gamma_v, sigma_v)
            
            % Make list of X coordinates, separate for clockwise and counterclockwise
            % vortices
            x_v_offset_array        = ((-VD.n_cells):1:(VD.n_cells)) * (2*VD.S);
            x_v_counter_array       =  z_v + x_v_offset_array;
            x_v_clockwise_array     = -z_v + x_v_offset_array;
            
            % Make list of Y coordinates, separate for clockwise and counterclockwise
            % vortices
            y_v_counter_array       =  y_v * ones(size(x_v_counter_array));
            y_v_clockwise_array     =  y_v * ones(size(x_v_clockwise_array));
            
            % Make list of strenghts, separate for clockwise and counterclockwise
            % vortices
            gamma_v_counter_array   =  gamma_v * ones(size(x_v_counter_array));
            gamma_v_clockwise_array = -gamma_v * ones(size(x_v_clockwise_array));
            
            % Make list of core radii, separate for clockwise and counterclockwise
            % vortices
            sigma_v_counter_array   =  sigma_v * ones(size(x_v_counter_array));
            sigma_v_clockwise_array =  sigma_v * ones(size(x_v_clockwise_array));
            
            % Concatenate into single list
            z_v_list     = [x_v_counter_array(:)     ;     x_v_clockwise_array(:)];
            y_v_list     = [y_v_counter_array(:)     ;     y_v_clockwise_array(:)];
            gamma_v_list = [gamma_v_counter_array(:) ; gamma_v_clockwise_array(:)];
            sigma_v_list = [sigma_v_counter_array(:) ; sigma_v_clockwise_array(:)];
            
        end
        
    end
    
    methods(Static)
        
        function [u, v] = induction_singular_vortex(z, y, z_v, y_v, gamma_v)
            % Computes Induction due to a single singular vortex at
            % multiple points
            
            % Norm
            r_sq_v_fun = @(x,y,x_v,y_v) (x-x_v).^2 + (y-y_v).^2;
            % x-component of induction
            u_v_fun = @(x,y,x_v,y_v, gamma_v) gamma_v / (2*pi()) * (-(y-y_v)) ./ r_sq_v_fun(x,y,x_v,y_v);
            % y-component of induction
            v_v_fun = @(x,y,x_v,y_v, gamma_v) gamma_v / (2*pi()) * (  x-x_v ) ./ r_sq_v_fun(x,y,x_v,y_v);
            
            
            u   = u_v_fun(z, y, z_v, y_v, gamma_v);
            v   = v_v_fun(z, y, z_v, y_v, gamma_v);
        end

        function [u, v] = induction_rankine_vortex(z, y, z_v, y_v, gamma_v, sigma_v)
            % Computes Induction due to a single rankine vortex at
            % multiple points
            
            %  Use Arrays for Non-Redundant Vectorization
            %  Square of Norm
            r_sq        =  (z-z_v).^2 + (y-y_v).^2 + eps();
            
            % x-component of laplace induction (singular kernel)
            u_laplace =  (-(y-y_v)) ./ r_sq;
            % y-component of laplace induction (singular kernel)
            v_laplace =  (  z-z_v ) ./ r_sq;
            
            % x-component of Rankine induction (solid rotation)
            u_rankine  = (-(y-y_v)) ./ sigma_v^2;
            % y-component of Rankine induction (solid rotation)
            v_rankine  = (  z-z_v ) ./ sigma_v^2;
            
            
            
            u   = gamma_v / (2*pi()) * (u_rankine  .* (r_sq <  sigma_v^2) + ...
                u_laplace  .* (r_sq >= sigma_v^2));
            v   = gamma_v / (2*pi()) * (v_rankine  .* (r_sq <  sigma_v^2) + ...
                v_laplace  .* (r_sq >= sigma_v^2));
        end
        
        function [u, v] = induction_lamb_vortex(z, y, z_v, y_v, gamma_v, sigma_v)
            %  Use Arrays for Non-Redundant Vectorization
            %  Square of Norm
            r_sq         =  (z-z_v).^2 + (y-y_v).^2 + eps();
            
            % x-component of laplace induction (singular kernel)
            u_laplace    =  (-(y-y_v)) ./ r_sq;
            % y-component of laplace induction (singular kernel)
            v_laplace    =  (  z-z_v ) ./ r_sq;
            
            % Lamb damping
            omega_max    = sigma_v;
            damping      =  1 - exp(- abs(pi() * omega_max .* r_sq.^2 ./ gamma_v));
            
            u   = gamma_v / (2*pi()) * (u_laplace) ;%.* damping;
            v   = gamma_v / (2*pi()) * (v_laplace) ;% .* damping;
            
        end
        
        function [u, v] = induction_lamb_vortex_debugging(z, y, z_v, y_v, gamma_v, sigma_v)
            %  Use Arrays for Non-Redundant Vectorization
            %  Square of Norm
            r_sq         =  (z-z_v).^2 + (y-y_v).^2 + eps();
            
            % x-component of laplace induction (singular kernel)
            u_laplace    =  (-(y-y_v)) ./ r_sq;
            % y-component of laplace induction (singular kernel)
            v_laplace    =  (  z-z_v ) ./ r_sq;
            
            % Lamb damping
            omega_max    = sigma_v;
            % Rewrite this to be dimensionally consistent (and results make
            % sense like this! we are out of time... verify this properly
            % ASAP and be transparent about it!) 
            % Original: (type)
            % damping      =  1 - exp(- abs(pi() * omega_max .* r_sq.^2 ./ gamma_v));
            % Corrected typo: (consistent with 2001 wendt paper, but dimensional: [Gamma]=m.s^-1 , [Omega]=s^-1, [r]=2))
            % r_sq is already squared!
            damping      =  1 - exp(- abs(pi() * r_sq .* omega_max ./ gamma_v));
            
            u   = gamma_v / (2*pi()) * (u_laplace) .* damping;%;
            v   = gamma_v / (2*pi()) * (v_laplace) .* damping;% ;
            
        end
        
        function [z_v_list, y_v_list, gamma_v_list, sigma_v_list] = ...
                duplicate_vortices_to_lower_halfplane(z_v_list, y_v_list, gamma_v_list, sigma_v_list)
            
            % Duplicate into lower halfplane (with symmetry, so gamma, inverts!)
            z_v_list     = [z_v_list(:)              ;                z_v_list(:)];
            y_v_list     = [y_v_list(:)              ;              - y_v_list(:)];
            gamma_v_list = [gamma_v_list(:)          ;          - gamma_v_list(:)];
            sigma_v_list = [sigma_v_list(:)          ;            sigma_v_list(:)];
            
        end
    end
    
end

