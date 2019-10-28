classdef field_mixer < handle
    %FIELD_MIXER is a handle class that centralizes context objects for
    %computing the rate of change of the mixing field
    %   Again, we stick with our architectural philosophy: objects store
    %   context but never state. Reaction to state should be independent of
    %   call history!
    
    properties

        VD                      % Vortex_Descriptor object
        CM                      % Crossflow_Mesh object
        SF                      % Shear_Field object
        MF                      % Mixing_Field object
        
        D_u_bar                 % Artificial Viscosity on Forcing Field (set to 0 to avoid very high diffusion near wall, unstabilizing)
        D_u_tilde               % Artificial Viscosity on C_tilde Field (0.002 seems nice for 100*100)
        baldwin_factor = 0      % Turbulent viscosity according to baldwin model (set to 0 for none, set to 1 for activating!)
        
    end
    
    methods
        % Constructor Function
        function FM = field_mixer(VD, CM, SF, MF, D_u_bar, D_u_tilde)
            % Store Inputs
            FM.VD = VD;
            FM.CM = CM;
            FM.SF = SF;
            FM.MF = MF;
            
            FM.D_u_bar   = D_u_bar;
            FM.D_u_tilde = D_u_tilde;            
        end
        
        function [u_tilde, u_bar, du_tilde_dx, dT_dx, mag_mesh , nut_mesh] = ...
                    u_tilde_rate_of_change(FM, state, nue)
            % Returns a u_tilde that is filtered with small fudgings but not
            % updated yet, and, most inportantly, the spatial rate of
            % change of u_tilde!
            
            % Extract data from state
            z_v         = state.z_v;
            y_v         = state.y_v;
            gamma_v     = state.gamma_v;
            sigma_v     = state.sigma_v;
            hk          = state.hk;
            rt          = state.rt;
            ue          = state.ue;
            u_tilde = state.u_tilde;
            ue_over_nue = state.ue / nue;
                                           
            % Get convection field (filtered by vortex_descriptor method itself!)
            [w_mesh, v_mesh, mag_mesh] = FM.VD.induction_on_mesh(FM.CM, z_v, y_v, gamma_v, sigma_v);
            
            % Get u_bar mesh (from Swafford profile)
            u_bar = FM.SF.u_bar_over_mesh(hk, rt, ue, ue_over_nue);
            
            % Filter u_bar, get its gradients and laplacian!
            [u_bar, du_bar_dz, du_bar_dy, lap_u_bar] = ...
                FM.SF.filtered_u_bar_gradients_laplacian(u_bar);
            
            % Filter u_tilde, get its gradients and laplacian!
            [u_tilde, du_tilde_dz, du_tilde_dy, lap_u_tilde] = ...
                FM.MF.filtered_u_tilde_gradients_laplacian(u_tilde);
            
            % Compute contribution to u_tilde due to convection of u_bar
            % (w*ddz + v*ddy)(u_bar)
            u_bar_convection_rate   = w_mesh.*du_bar_dz   + v_mesh.*du_bar_dy;
            % Compute contribution to u_tilde due to convection of u_tilde
            % (w*ddz + v*ddy)(u_tilde)
            u_tilde_convection_rate = w_mesh.*du_tilde_dz + v_mesh.*du_tilde_dy;
            
            % Compute turbulent nut
            if FM.baldwin_factor == 1
                nut_mesh = nut_baldwin_lomax1978(FM, state, nue, u_bar, du_bar_dy, du_tilde_dz, du_tilde_dy, w_mesh, v_mesh);
            elseif FM.baldwin_factor == 2
                nut_mesh = nut_prandtl_brederode(FM, state, nue, u_bar, du_bar_dy, du_tilde_dz, du_tilde_dy, w_mesh, v_mesh);
            else
                nut_mesh = 0;
            end
            % Compute contribution due to diffusion of u_tilde and u_bar
            u_bar_diffusion_rate   = FM.D_u_bar   * lap_u_bar;      % This is deactivated with D_u_bar = 0 , but it should disappear! (There is no physical reason for this in current model/deduction!)
            u_tilde_diffusion_rate = (nue + nut_mesh + FM.D_u_tilde) .* lap_u_tilde;    % FM.D_u_tilde corresponds to a nu for artificial diffusion/dissipation, helping stabilize things!
            
%             % Experiment (rough...)
%             von_karman                 = 0.41;
%             y_over_delta_threshold_bot = 0.5;
%             y_over_delta_threshold_top = 3;  
%             nu_turb_max = 0.001;
%             l_prandtl 		  = von_karman * FM.CM.y_mesh.^2;
% %             y_mesh_mod = FM.CM.y_mesh.^2;
% %             y_mesh_mod(FM.CM.y_range < y_over_delta_threshold_bot*FM.CM.delta,:) = 0;
% %             y_mesh_mod(FM.CM.y_range > y_over_delta_threshold_top*FM.CM.delta,:) = 0;
% %             l_prandtl              = von_karman * y_mesh_mod ;
%             nu_turb                = l_prandtl .* abs(lap_u_tilde);
%             nu_turb(nu_turb>nu_turb_max) = nu_turb_max;
%             u_turb_diffusion_rate  = nu_turb .* lap_u_tilde;
%             u_turb_diffusion_rate(1  , :  ) = 0;
%             u_turb_diffusion_rate(end, :  ) = 0;
            
            %u_tilde_diffusion_rate = u_tilde_diffusion_rate + u_turb_diffusion_rate;
            
            
            % Now elliminate diffusion from edges (this has been moved to
            % MF.filtered_u_tilde_gradients_laplacian and
            % SF.filtered_u_bar_gradients_laplacian)
            % diffusion_rate           = u_bar_diffusion_rate + u_tilde_diffusion_rate;
            % diffusion_rate(1  , :  ) = 0;
            % diffusion_rate(end, :  ) = 0;
            
            
            % Now compute partial time derivative of forcing field (0 for now)
            % ddt(u_bar)
            du_bar_dt = zeros(size(FM.CM.z_mesh));
            
            % Almost finally, make rate of change of u_tilde in lagrangian time (in x)
            du_tilde_dT = - (du_bar_dt + u_bar_convection_rate + u_tilde_convection_rate) + (u_bar_diffusion_rate + u_tilde_diffusion_rate);
            
            % Now make relation between lagrangian time (in x) and space (mesh)
            dT_dx = 1 ./ u_bar;
            % Stabilize first line (goes to inf because 0=u_bar(1,:))
            dT_dx(1,:) = dT_dx(2,:);
            
            % Finally, make rate of change of u_tilde in space (in x)
            du_tilde_dx = dT_dx .* du_tilde_dT;
            
        end
        
        function nut_mesh = nut_baldwin_lomax1978(FM, state, nue, u_bar, du_bar_dy, du_tilde_dz, du_tilde_dy, w_mesh, v_mesh)
            % Returns a nut mesh computed according to the original
            % Baldwin-Lomax 1978 paper!
            %
            % Variant based on intersections (this is the main, preffered
            % and reference variant, faster than the one based on InterX
            % and much more tested by hand!)
            
            % % Define Constants
            A_plus = 26;
            C_cp   = 1.6;
            C_kleb = 0.3;
            C_wk   = 0.25;
            k_small= 0.4;
            K_big  = 0.0168;
            %p_r    = 0.72;     % Unused?
            %p_rt   = 0.9;      % Unused?
            %C_mutm = 14;       % For transition
            
            % % Generate y+ mesh
            % Start by computing re_theta of 
            rt_bar = (state.theta * state.ue) / nue;
            % To compute skin friciton of shear field
            [ cf_bar] = cft_rr( state.hk, rt_bar, 0);
            % And now define the skin friction velocity (ignore spanwise
            % dependency for now, can be added later and shouldn't change
            % things too much!)
            u_tau = state.ue * sqrt(0.5 * cf_bar);
            % Now make y+ mesh
            yplus_mesh = (u_tau / nue) .* FM.CM.y_mesh;
            
            % % Compute mixing lenght for inner layer (Prandtl-Van driest)
            % Now make y+ mesh
            l_prandtl = k_small * FM.CM.y_mesh .* (1 - exp( - yplus_mesh ./ A_plus));
            
            % % Compute mesh of vorticity vector magnitude
            % Get field for debug
            % [w_mesh, v_mesh, mag_mesh] = FM.VD.induction_on_mesh(FM.CM, state.z_v, state.y_v, state.gamma_v, state.sigma_v)
            % [u_tilde, du_tilde_dz, du_tilde_dy, lap_u_tilde] = FM.MF.filtered_u_tilde_gradients_laplacian(state.u_tilde)
            % u_bar = FM.SF.u_bar_over_mesh(state.hk, rt_bar, state.ue, state.ue/nue);
            % [u_bar, du_bar_dz, du_bar_dy, lap_u_bar] = FM.SF.filtered_u_bar_gradients_laplacian(u_bar);
            % Get Gradients of v_tilde and w_tilde fields
            [~           , dw_tilde_dy] = gradient(w_mesh, FM.CM.z_range, FM.CM.y_range);
            [dv_tilde_dz , ~          ] = gradient(v_mesh, FM.CM.z_range, FM.CM.y_range);
            % Compute norm of vorticity vector (gets really large on vortices!)
            omega_norm_mesh = sqrt((dw_tilde_dy - dv_tilde_dz).^2 + du_tilde_dz.^2 + (du_bar_dy + du_tilde_dy).^2);
            % Optionally, smoothen vorticity norm by convolving a kernel on it!
%             K_smooth = 0.045*ones(5);
%             omega_norm_mesh_smooth = conv2(omega_norm_mesh ,K_smooth,'same');
%             omega_norm_mesh = omega_norm_mesh_smooth;
            
            
            
            % Now compute turbulent viscosity (nu_t) for the inner part of
            % the BL (law of the wall region!)
            nut_inner = l_prandtl.^2 .* omega_norm_mesh;
            
            % % Now move on to the outer layer
            % Compute big F(y) function
            f_of_y_mesh = FM.CM.y_mesh .* omega_norm_mesh .* (1 - exp(- yplus_mesh ./ A_plus));
            % Determine max of f for each z stance
            [f_max, i_max] = max(f_of_y_mesh, [], 1);
            % Determine corresponding y position
            y_max_range = FM.CM.y_range(i_max);
            % Optionally smoothen curve a bit (this seems like a good idea, for this, do it twice!)
            y_max_range_smooth  = smoothdata(y_max_range        , 'movmean', round(5 * 101 / length(FM.CM.z_range)));
            y_max_range_smooth2 = smoothdata(y_max_range_smooth , 'movmean', round(5 * 101 / length(FM.CM.z_range)));
            y_max_range_smooth3 = smoothdata(y_max_range_smooth2, 'movmean', round(5 * 101 / length(FM.CM.z_range)));
            y_max_range = y_max_range_smooth3;
            % And replicate it over mesh (from range with z dependency)
            y_max_mesh = zeros(size(FM.CM.y_mesh));
            for j=1:size(y_max_mesh,2)
                y_max_mesh(:,j) = y_max_range(j);
            end
            % Some plotting code for verification (locci are a bit harsh)
            %surf(FM.CM.z_mesh, FM.CM.y_mesh, f_of_y); xlabel('z'); ylabel('y'); zlabel('F of y'); shading flat; hold on
            %plot3(FM.CM.z_range, y_max , f_max, 'ro-')
            % Now compute u_dif (for now assume vortices induce less speed
            % than edge velocity!... make things a bit more stable)
            % u_dif = state.ue;
            % Or do it as it is supposed!
            u_mag = sqrt((u_bar+state.u_tilde).^2 + w_mesh.^2 + v_mesh.^2);
            u_dif = max(u_mag , [], 1);
            % Now compute F_wake components
            F_wake_range_1 = y_max_range .* f_max;
            F_wake_range_2 = C_wk .* y_max_range .* u_dif.^2 ./ f_max;
            % Finally compute complete F_wake
            F_wake_range         = transpose(min([F_wake_range_1(:) , F_wake_range_2(:)], [], 2));
            % Optionally smoothen curve a bit (this seems like a good idea)
            % F_wake_range_smooth = smoothdata(F_wake_range);
            F_wake_range_smooth = smoothdata(F_wake_range, 'movmean', round(5 * 101 / length(FM.CM.z_range)));
            F_wake_range = F_wake_range_smooth;
            % Or switch back to original clauser formula for F_wake_range_1 (first component only theoretically!)
            %F_wake_range = state.hk * state.theta * state.ue * ones(size(FM.CM.z_range))
            % And replicate it over mesh (from range with z dependency)
            F_wake_mesh = zeros(size(FM.CM.y_mesh));
            for j=1:size(F_wake_mesh,2)
                F_wake_mesh(:,j) = F_wake_range(j); 
            end
            %plot(FM.CM.z_range, F_wake_range_1, FM.CM.z_range, F_wake_range_2, FM.CM.z_range, F_wake_range, FM.CM.z_range,F_wake_range_smooth, 'k') ; grid on; hold on; plot(FM.CM.z_range, state.hk * state.theta * state.ue * ones(size(FM.CM.z_range)), 'k--')
            %legend('Component 1 (replace clauser)', 'Component 2 (mostly for wakes)', 'Raw F-wake Baldwin-Lomax', 'Smooth F-wake Baldwin-Lomax', 'Original Clauser')
            % Now define the Klebanoff intermittency factor over mesh
            F_kleb_mesh = (1 + 5.5*(C_kleb*FM.CM.y_mesh./y_max_mesh).^6 ).^(-1);
            % And finally write nut_outer (viscosity in outer layer! nut = mut/rho)
            nut_outer = K_big * C_cp * F_wake_mesh .* F_kleb_mesh;
            
            % Now comes the time to find intersections in y (find yc, the
            % crossover point at each z stance) to make the nut_mesh  
            %tic
            nut_mesh = zeros(size(nut_outer));
            for j=1:size(nut_mesh ,2)
                [y_intersect, nut_intersect] = intersections(FM.CM.y_range, nut_inner(:,j), FM.CM.y_range, nut_outer(:,j)); %#ok<ASGLU>
                flag_keep_inner = FM.CM.y_range<min(y_intersect);
                flag_keep_outer = not(flag_keep_inner);
                i_keep_inner = find(flag_keep_inner);
                i_keep_outer = find(flag_keep_outer);
                nut_mesh(i_keep_inner,j) = nut_inner(i_keep_inner,j);
                nut_mesh(i_keep_outer,j) = nut_outer(i_keep_outer,j);
                %plot(FM.CM.y_range, nut_inner(:,j), FM.CM.y_range, nut_outer(:,j)); hold on; plot(y_intersect,nut_intersect, 'x');plot(FM.CM.y_range,nut_mesh(:,j), 'g--'); legend('Inner', 'Outer', 'Intersections', 'Joined'); grid on;
            end
            %toc
            %figure(1); surf(FM.CM.z_mesh, FM.CM.y_mesh, nut_inner / D_u_tilde); xlabel('z'); ylabel('y'); zlabel('F of y'); shading flat; hold on; view(2); colorbar
            %figure(2); surf(FM.CM.z_mesh, FM.CM.y_mesh, nut_outer / D_u_tilde); xlabel('z'); ylabel('y'); zlabel('F of y'); shading flat; hold on; view(2); colorbar
            %figure(3); surf(FM.CM.z_mesh, FM.CM.y_mesh, nut_mesh  / D_u_tilde); xlabel('z'); ylabel('y'); zlabel('F of y'); shading flat; hold on; view(2); colorbar
            % Some artifacts remain... let's see how it goes!
        end
        
        function nut_mesh = nut_baldwin_lomax1978_interX(FM, state, nue, u_bar, du_bar_dy, du_tilde_dz, du_tilde_dy, w_mesh, v_mesh)
            % Returns a nut mesh computed according to the original
            % Baldwin-Lomax 1978 paper!
            %
            % Variant based on interX (supposedly faster, but actually
            % slower than the one based on intersectoins, (Matlab changed
            % quite a bit since both were written!)
            
            % % Define Constants
            A_plus = 26;
            C_cp   = 1.6;
            C_kleb = 0.3;
            C_wk   = 0.25;
            k_small= 0.4;
            K_big  = 0.0168;
            %p_r    = 0.72;     % Unused?
            %p_rt   = 0.9;      % Unused?
            %C_mutm = 14;       % For transition
            
            % % Generate y+ mesh
            % Start by computing re_theta of 
            rt_bar = (state.theta * state.ue) / nue;
            % To compute skin friciton of shear field
            [ cf_bar] = cft_rr( state.hk, rt_bar, 0);
            % And now define the skin friction velocity (ignore spanwise
            % dependency for now, can be added later and shouldn't change
            % things too much!)
            u_tau = state.ue * sqrt(0.5 * cf_bar);
            % Now make y+ mesh
            yplus_mesh = (u_tau / nue) .* FM.CM.y_mesh;
            
            % % Compute mixing lenght for inner layer (Prandtl-Van driest)
            % Now make y+ mesh
            l_prandtl = k_small * FM.CM.y_mesh .* (1 - exp( - yplus_mesh ./ A_plus));
            
            % % Compute mesh of vorticity vector magnitude
            % Get field for debug
            % [w_mesh, v_mesh, mag_mesh] = FM.VD.induction_on_mesh(FM.CM, state.z_v, state.y_v, state.gamma_v, state.sigma_v)
            % [u_tilde, du_tilde_dz, du_tilde_dy, lap_u_tilde] = FM.MF.filtered_u_tilde_gradients_laplacian(state.u_tilde)
            % u_bar = FM.SF.u_bar_over_mesh(state.hk, rt_bar, state.ue, state.ue/nue);
            % [u_bar, du_bar_dz, du_bar_dy, lap_u_bar] = FM.SF.filtered_u_bar_gradients_laplacian(u_bar);
            % Get Gradients of v_tilde and w_tilde fields
            [~           , dw_tilde_dy] = gradient(w_mesh, FM.CM.z_range, FM.CM.y_range);
            [dv_tilde_dz , ~          ] = gradient(v_mesh, FM.CM.z_range, FM.CM.y_range);
            % Compute norm of vorticity vector (gets really large on vortices!)
            omega_norm_mesh = sqrt((dw_tilde_dy - dv_tilde_dz).^2 + du_tilde_dz.^2 + (du_bar_dy + du_tilde_dy).^2);
            % Optionally, smoothen vorticity norm by convolving a kernel on it!
%             K_smooth = 0.045*ones(5);
%             omega_norm_mesh_smooth = conv2(omega_norm_mesh ,K_smooth,'same');
%             omega_norm_mesh = omega_norm_mesh_smooth;
            
            
            
            % Now compute turbulent viscosity (nu_t) for the inner part of
            % the BL (law of the wall region!)
            nut_inner = l_prandtl.^2 .* omega_norm_mesh;
            
            % % Now move on to the outer layer
            % Compute big F(y) function
            f_of_y_mesh = FM.CM.y_mesh .* omega_norm_mesh .* (1 - exp(- yplus_mesh ./ A_plus));
            % Determine max of f for each z stance
            [f_max, i_max] = max(f_of_y_mesh, [], 1);
            % Determine corresponding y position
            y_max_range = FM.CM.y_range(i_max);
            % Optionally smoothen curve a bit (this seems like a good idea, for this, do it twice!)
            y_max_range_smooth  = smoothdata(y_max_range        , 'movmean', round(5 * 101 / length(FM.CM.z_range)));
            y_max_range_smooth2 = smoothdata(y_max_range_smooth , 'movmean', round(5 * 101 / length(FM.CM.z_range)));
            y_max_range_smooth3 = smoothdata(y_max_range_smooth2, 'movmean', round(5 * 101 / length(FM.CM.z_range)));
            y_max_range = y_max_range_smooth3;
            % And replicate it over mesh (from range with z dependency)
            y_max_mesh = zeros(size(FM.CM.y_mesh));
            for j=1:size(y_max_mesh,2)
                y_max_mesh(:,j) = y_max_range(j);
            end
            % Some plotting code for verification (locci are a bit harsh)
            %surf(FM.CM.z_mesh, FM.CM.y_mesh, f_of_y); xlabel('z'); ylabel('y'); zlabel('F of y'); shading flat; hold on
            %plot3(FM.CM.z_range, y_max , f_max, 'ro-')
            % Now compute u_dif (for now assume vortices induce less speed
            % than edge velocity!... make things a bit more stable)
            % u_dif = state.ue;
            % Or do it as it is supposed!
            u_mag = sqrt((u_bar+state.u_tilde).^2 + w_mesh.^2 + v_mesh.^2);
            u_dif = max(u_mag , [], 1);
            % Now compute F_wake components
            F_wake_range_1 = y_max_range .* f_max;
            F_wake_range_2 = C_wk .* y_max_range .* u_dif.^2 ./ f_max;
            % Finally compute complete F_wake
            F_wake_range         = transpose(min([F_wake_range_1(:) , F_wake_range_2(:)], [], 2));
            % Optionally smoothen curve a bit (this seems like a good idea)
            % F_wake_range_smooth = smoothdata(F_wake_range);
            F_wake_range_smooth = smoothdata(F_wake_range, 'movmean', round(5 * 101 / length(FM.CM.z_range)));
            F_wake_range = F_wake_range_smooth;
            % Or switch back to original clauser formula for F_wake_range_1 (first component only theoretically!)
            %F_wake_range = state.hk * state.theta * state.ue * ones(size(FM.CM.z_range))
            % And replicate it over mesh (from range with z dependency)
            F_wake_mesh = zeros(size(FM.CM.y_mesh));
            for j=1:size(F_wake_mesh,2)
                F_wake_mesh(:,j) = F_wake_range(j); 
            end
            %plot(FM.CM.z_range, F_wake_range_1, FM.CM.z_range, F_wake_range_2, FM.CM.z_range, F_wake_range, FM.CM.z_range,F_wake_range_smooth, 'k') ; grid on; hold on; plot(FM.CM.z_range, state.hk * state.theta * state.ue * ones(size(FM.CM.z_range)), 'k--')
            %legend('Component 1 (replace clauser)', 'Component 2 (mostly for wakes)', 'Raw F-wake Baldwin-Lomax', 'Smooth F-wake Baldwin-Lomax', 'Original Clauser')
            % Now define the Klebanoff intermittency factor over mesh
            F_kleb_mesh = (1 + 5.5*(C_kleb*FM.CM.y_mesh./y_max_mesh).^6 ).^(-1);
            % And finally write nut_outer (viscosity in outer layer! nut = mut/rho)
            nut_outer = K_big * C_cp * F_wake_mesh .* F_kleb_mesh;
            
            % Now comes the time to find intersections in y (find yc, the
            % crossover point at each z stance) to make the nut_mesh  
            %tic
            nut_mesh = zeros(size(nut_outer));
            for j=1:size(nut_mesh ,2)
                %[y_intersect, nut_intersect] = intersections(FM.CM.y_range, nut_inner(:,j), FM.CM.y_range, nut_outer(:,j)); %#ok<ASGLU>
                [P_intersect] = InterX([FM.CM.y_range(:), nut_inner(:,j)]', [FM.CM.y_range(:), nut_outer(:,j)]');
                y_intersect = P_intersect(1,:)'; 
                flag_keep_inner = FM.CM.y_range<min(y_intersect);
                flag_keep_outer = not(flag_keep_inner);
                i_keep_inner = find(flag_keep_inner);
                i_keep_outer = find(flag_keep_outer);
                nut_mesh(i_keep_inner,j) = nut_inner(i_keep_inner,j);
                nut_mesh(i_keep_outer,j) = nut_outer(i_keep_outer,j);
                %plot(FM.CM.y_range, nut_inner(:,j), FM.CM.y_range, nut_outer(:,j)); hold on; plot(y_intersect,nut_intersect, 'x');plot(FM.CM.y_range,nut_mesh(:,j), 'g--'); legend('Inner', 'Outer', 'Intersections', 'Joined'); grid on;
            end
            %toc
            %figure(1); surf(FM.CM.z_mesh, FM.CM.y_mesh, nut_inner / D_u_tilde); xlabel('z'); ylabel('y'); zlabel('F of y'); shading flat; hold on; view(2); colorbar
            %figure(2); surf(FM.CM.z_mesh, FM.CM.y_mesh, nut_outer / D_u_tilde); xlabel('z'); ylabel('y'); zlabel('F of y'); shading flat; hold on; view(2); colorbar
            %figure(3); surf(FM.CM.z_mesh, FM.CM.y_mesh, nut_mesh  / D_u_tilde); xlabel('z'); ylabel('y'); zlabel('F of y'); shading flat; hold on; view(2); colorbar
            % Some artifacts remain... let's see how it goes!
        end
        
        function nut_mesh = nut_prandtl_brederode(FM, state, nue, u_bar, du_bar_dy, du_tilde_dz, du_tilde_dy, w_mesh, v_mesh)
            % Returns a nut mesh computed according to the Prandt model
            % with damping given by Brederode (p. 274) . Shear generalize
            % after Baldwin-Lomax 1978 

            
            % % Define Constants
            rho   = 1.225;
            k_big = 0.41;
            y_over_delta_threshold = 0.02;
            
            
            % % Compute mixing lenght for inner layer (Prandtl-Van driest)
            % Now make lm = K*y mesh (p.274 Brederode)
            l_prandtl = k_big * FM.CM.y_mesh;
            % And fudge it above 0.2*delta
            l_prandtl(FM.CM.y_mesh > y_over_delta_threshold) = 0.082;
            
            % Get Gradients of v_tilde and w_tilde fields
            [~           , dw_tilde_dy] = gradient(w_mesh, FM.CM.z_range, FM.CM.y_range);
            [dv_tilde_dz , ~          ] = gradient(v_mesh, FM.CM.z_range, FM.CM.y_range);
            % Compute norm of vorticity vector (gets really large on vortices!)
            % omega_norm_mesh = sqrt((dw_tilde_dy - dv_tilde_dz).^2 + du_tilde_dz.^2 + (du_bar_dy + du_tilde_dy).^2);
            omega_norm_mesh = sqrt((dw_tilde_dy - dv_tilde_dz).^2 + (du_bar_dy).^2);
            
            % Now compute turbulent viscosity according to brederode
            nut_mesh = rho * l_prandtl.^2 .* omega_norm_mesh;
            
            
            
        end
    end
    
end

