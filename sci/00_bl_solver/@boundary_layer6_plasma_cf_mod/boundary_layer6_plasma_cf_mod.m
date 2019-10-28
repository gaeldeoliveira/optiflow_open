classdef boundary_layer6_plasma_cf_mod < boundary_layer6_plasma
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
        
% Super-Superclass
% %         % Viscosity 
% %         nue                             % edge viscosity        
% %         
% %         % Forcing Terms
% %         xi_vector                       % x positions definition vector
% %         ue_vector                       % ue (edge speed) definition vector
% %         msq_vector                      % msq (mach) definition vector
% %         ue_over_nue_vector              % ue_over_nue (convective ratio) definition vector
% %         
% %         % Piecewise Polynomial Interpolants and their Piecewise Polynomial
% %         % Derivatives of the Forcing Terms
% %         pp_ue                           % edge velocity interpolant (C1 continuity, piecewise C2)
% %         pp_due_dxi                      % edge velocity derivative interpolant (C0 continuity, piecewise C1)
% %         pp_msq                          % mach number interpolant (C1 continuity, piecewise C2)
% %         pp_dmsq_dxi                     % mach number derivative interpolant (C0 continuity, piecewise C1)
% %         pp_ue_over_nue                  % convective ratio interpolant (C1 continuity, piecewise C2)
% %         pp_due_over_nue_dxi             % convective ratio derivative interpolant (C0 continuity, piecewise C1)        
% %         
% %         % Model Options (as in Rfoil)
% %         scci = 5.60                     % Shear lag constant (add gacon and gbcon in the future!)

% Superclass
%         % Handles to Objects
%         PD_cell                           % Cell Array of Handles to plasma_descriptor objects
%         CEI                               % Handle cei_closure_manager object

    end
    
    methods
        function BL = boundary_layer6_plasma_cf_mod()
            % Constructor Function
        end
                
        function [rhs_theta, rhs_hs, rhs_ctau] = canonical_rhs(BL, t, h, ctau, x)
            % Typical Test Values
            % t       = 0.002;
            % dstr    = 0.003;
            % ctau    = 0.06*0.06;         % Max shear stress coefficient
            
            % Get interpolated forcing terms
            ue                  = BL.ue_function(x);
            due_dxi             = BL.due_dxi_function(x);
            msq                 = BL.msq_function(x);
            dmsq_dxi            = BL.dmsq_dxi_function(x); %#ok<*NASGU>
            ue_over_nue         = BL.ue_over_nue_function(x);
            due_over_nue_dxi    = BL.due_over_nue_dxi_function(x);
            
            % % % Start Calculation
            % Shape Factor
%            [ h, h_dstr, h_t] = H(dstr , t);                           % Looks ok
            dstr = h * t;
            % Kinematic Shape Factor
            [ hk, hk_h, hk_msq  ] = hkin( h, msq);                      %#ok<*ASGLU> % Looks ok
            
            % Reynolds Theta
            [rt, rt_ue_over_nue ] = re_theta( t, ue_over_nue );         % Looks ok
            
            % Now calculate Cf
            [ cf, cf_hk, cf_rt, cf_msq  ] = cft_rr( hk, rt, msq);       % Looks ok
            
            % Sum plasma contribution(s) to cf (there is never more than one at a
            % time, in practice!)
            for n_PD = 1:length(BL.PD_cell)

                % Compute Longitudinal plasma weighting function at this point
                wx          = BL.PD_cell{n_PD}.wx_function_CFM(x);
                
                % Retrive Average Plasma Force Field Density Coefficient                
                c_phi_x_p   = BL.PD_cell{n_PD}.c_phi_x_p;
                
                % Compute Force Field Momentum Coefficient
                cfm         = wx .* c_phi_x_p;
            
                % Couette Flow Correction
                % cf = cf + cfm ./ ue.^2;
                % Pipe Flow Correction (deduced!)
                cf = cf + (2/pi()) * cfm ./ ue.^2;
                
            end
            
            % % Momentum Equation RHS Evaluation Draft
            rhs_theta = cf/2 - (h + 2 - msq.^2) .* (t ./ ue) .*  due_dxi;
            
            % Now energy shape factor (H_star)
            [ hs, hs_hk, hs_rt, hs_msq  ] = hst( hk, rt, msq);          % Check not possible
            
            % Slip Speed (necessary for dissipation term)
            [ us, us_h, us_hk, us_hs] = usg(h, hk, hs);                 % Check not possible
            
            % Now proceed to dissipation term
            [ cdi, cdi_us, cdi_cf, cdi_ctau ] = cdissipation( us, cf, ctau);           % based on merchant, requires verification but makes sense!
            
            % Density shape parameter
            [ hc, hc_hk, hc_msq ] = hct( hk, msq );                      % Check not possible, 0 value for mach 0 seems ok!
            
            % % Energy Equation RHS Evaluation Draft (corrected)
            % rhs_hs =  2 * cdi - hs .* (cf / 2) - (2 * hc + hs.*(1 - h)) .* (t ./ ue) .* due_dxi;
            rhs_hs =  (2 * cdi - hs .* (cf / 2) - (2 * hc + hs.*(1 - h)) .* (t ./ ue) .* due_dxi) ./ t;
            
            % Compute Equilibrium Shear Coefficient
            [ ctz, ctz_h , ctz_hk , ctz_hs , ctz_us ] = ctauzero( h , hk , hs , us);
            
            % Head's h1 Shape Parameter
            [ hh, hh_hk ] = hh1cal( hk );
            
            % Boundary Layer Thickness
            [ d , d_h, d_hh, d_t ] = delta_hh( h, hh, t);
            
            
            % Set shear lag constant
            scc = BL.scci;      % (5.60, default!)
            % SCC = SCCI-.950-.950*tanh(.275*hk2-3.5) Choice option May-1998. New 20 june 1996
            
            % % Shear Lag Equation
            rhs_ctau = (ctau ./ d) .* scc .* (sqrt(ctz) - sqrt(ctau));
%            rhs_ctau = (ctau.^2 ./ d) .* scc .* (ctz - ctau);
        end                        
        
    end
    
end

