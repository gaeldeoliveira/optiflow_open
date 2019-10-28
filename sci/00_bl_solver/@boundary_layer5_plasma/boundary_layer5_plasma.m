classdef boundary_layer5_plasma < boundary_layer5
    %boundary_layer3  is a handle class definning a boundary_layer model  
    %   similar to that of Rfoil (without rotational effects, and without 
    %   suction or plasma).                                               
    %   Pressure gradient is handled correctly (formally speaking) from   
    %   this version onwards.                                             
    %                                                                     
    %   It also provides for the storage, interpolation (pchip) and       
    %   accurate differentiation of the forcing terms.                    
    %                                                                     
    %   It sustains both the evaluation of the canonic equations RHS      
    %   and the direct system RHS formulated in terms of {theta, h, ctau} 
    %   Finally, the ODE_FUN_WRAPPER method provides an interface to      
    %   Matlab's standard ODE solvers (eg. ode45, ode23)                  
    %                                                                     
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
    %       Boundary Layer ODE RHSs with Plasma Terms
    %           following:
    %               Modelling the effect of DBD plasma actuators on boundary
    %               layer development
    %               Internal Report, Gael de Oliveira, Ricardo Pereira
    %               September 8, 2014
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties        
%         % Viscosity 
%         nue                             % edge viscosity        
%         
%         % Forcing Terms
%         xi_vector                       % x positions definition vector
%         ue_vector                       % ue (edge speed) definition vector
%         msq_vector                      % msq (mach) definition vector
%         ue_over_nue_vector              % ue_over_nue (convective ratio) definition vector
%         
%         % Piecewise Polynomial Interpolants and their Piecewise Polynomial
%         % Derivatives of the Forcing Terms
%         pp_ue                           % edge velocity interpolant (C1 continuity, piecewise C2)
%         pp_due_dxi                      % edge velocity derivative interpolant (C0 continuity, piecewise C1)
%         pp_msq                          % mach number interpolant (C1 continuity, piecewise C2)
%         pp_dmsq_dxi                     % mach number derivative interpolant (C0 continuity, piecewise C1)
%         pp_ue_over_nue                  % convective ratio interpolant (C1 continuity, piecewise C2)
%         pp_due_over_nue_dxi             % convective ratio derivative interpolant (C0 continuity, piecewise C1)        
%         
%         % Model Options (as in Rfoil)
%         scci = 5.60                     % Shear lag constant (add gacon and gbcon in the future!)

        % Handles to Objects
        PD_cell                           % Cell Array of Handles to plasma_descriptor objects
        CEI                               % Handle cei_closure_manager object

    end
    
    methods
        function BL = boundary_layer5_plasma()
            % Constructor Function
        end
                
        function [rhs_theta, rhs_h, rhs_ctau] = transformed_rhs(BL, t, h, ctau, x)
            % % Get canonical rhs
            [rhs_theta_c, rhs_hs_c, rhs_ctau_c] = canonical_rhs(BL, t, h, ctau, x);
            % % Add Plasma terms to it
            [rhs_theta_p, rhs_hs_p, rhs_ctau_p] = add_plasma_to_rhs(BL, t, h, ctau, x, rhs_theta_c, rhs_hs_c, rhs_ctau_c);
            % % Transform rhs
            [rhs_theta, rhs_h, rhs_ctau] = transform_rhs(BL, t, h, ctau, x, rhs_theta_p, rhs_hs_p, rhs_ctau_p);
            % Done! Return!
        end
        
        function [rhs_theta_p, rhs_hs_p, rhs_ctau_p] = add_plasma_to_rhs(BL, t, h, ctau, x, rhs_theta, rhs_hs, rhs_ctau) %#ok<INUSL>
            % Add plasma terms to RHS

            % Get interpolated forcing terms
            ue          = BL.ue_function(x);
            msq         = BL.msq_function(x);
            ue_over_nue = BL.ue_over_nue_function(x);
            
            % Initialize RHS returns as if there were no plasmas
            rhs_theta_p = rhs_theta;
            rhs_hs_p    = rhs_hs;
            rhs_ctau_p = rhs_ctau;
            
            
            % Sum plasma contribution(s) (there is never more than one at a
            % time, in practice!)
            for n_PD = 1:length(BL.PD_cell)

                % Compute Longitudinal plasma weighting function at this point
                wx          = BL.PD_cell{n_PD}.wx_function_CFM(x);
                
                % Retrive Average Plasma Force Field Density Coefficient
                c_phi_x_p   = BL.PD_cell{n_PD}.c_phi_x_p;
                
                % Compute Force Field Momentum Coefficient
                cfm         = wx .* c_phi_x_p;
                
                % Now proceed to write momentum rhs with plasma term
                rhs_theta_p = rhs_theta_p - 0.5 * cfm ./ ue.^2;
                
                % Kinematic Shape Factor
                hk          = hkin( h, msq);                      % Looks ok
                
                % Reynolds Theta
                rt          = re_theta( t, ue_over_nue );         % Looks ok
                
                % Now energy shape factor (H_star)
                hs          = hst( hk, rt, msq);          % Check not possible
                
                % Compute Momentum Scaled Plasma Force Field Thickness
                % (t_theta_p = (T_p / delta_2) = (t_p / theta)
                t_theta_p   = BL.PD_cell{n_PD}.t_p ./ t;
                
                % Now retrieve Energy Interaction Coefficient from closure
                % dataset
                % (consistency argument, but weak, deduction done with h, not
                % hk, but Drela seems to use Hk everywhere! Check when time available!)
                cei = BL.CEI.cei_function(hk, rt, t_theta_p);
                
                % Now compute Force Energy Coefficient
                cfe = cfm .* cei;
                
                % And proceed to write RHS of Energy Shape Factor Equation
                rhs_hs_p    = rhs_hs_p  +  ((0.5 * hs .* cfm)./ ue.^2)  -  (cfe ./ ue.^3);
            
            end
            
        end
                               
        function [y_line] = ode_fun_wrapper(BL, x, y)
            % Wrapper function of BL.transformed_rhs method for matlab
            % standard ode fun interface
            
            t       = y(1);
            dstr    = y(2);
            ctau    = y(3);
            
            [rhs_theta, rhs_dstr, rhs_ctau] = transformed_rhs(BL, t, dstr, ctau, x);
            
            y_line = zeros(length(y),1);
            y_line(1) = rhs_theta;
            y_line(2) = rhs_dstr;
            y_line(3) = rhs_ctau;
        end
        
    end
    
end

