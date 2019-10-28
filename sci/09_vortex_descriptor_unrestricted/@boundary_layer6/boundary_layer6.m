classdef boundary_layer6 < handle
    %boundary_layer6  is a handle class definning a boundary_layer model  
    %   similar to that of Rfoil (without rotational effects, and without 
    %   suction or plasma).                                                                                            
    %                                                                     
    %   The class provides for the storage, interpolation (pchip) and       
    %   accurate differentiation of the forcing terms.                    
    %                                                                     
    %   It sustains both the evaluation of the canonic equations RHS      
    %   and the direct system RHS formulated in terms of {theta, h, ctau} 
    %
    %   The ODE_FUN_WRAPPER method provides an interface to      
    %   Matlab's standard ODE solvers (eg. ode45, ode23)                  
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    %   Release Notes:
    %
    %       Pressure gradient is handled correctly (formally speaking) from
    %       version 3 onwards.
    %
    %       ODE Integration is integrated from version 6 onwards
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
    
    properties
        % Viscosity 
        u_inf_over_nu_inf               % viscosity scale (from Re)
        
        % Forcing Terms
        xi_vector                       % x positions definition vector
        ue_vector                       % ue (edge speed) definition vector  (adimensionalized to Uinf, that is ue = Ue / Uinf)
        msq_vector                      % msq (mach) definition vector
        ue_over_nue_vector              % ue_over_nue (convective ratio) definition vector
        
        % Piecewise Polynomial Interpolants and their Piecewise Polynomial
        % Derivatives of the Forcing Terms
        pp_ue                           % edge velocity interpolant (C1 continuity, piecewise C2)
        pp_due_dxi                      % edge velocity derivative interpolant (C0 continuity, piecewise C1)
        pp_msq                          % mach number interpolant (C1 continuity, piecewise C2)
        pp_dmsq_dxi                     % mach number derivative interpolant (C0 continuity, piecewise C1)
        pp_ue_over_nue                  % convective ratio interpolant (C1 continuity, piecewise C2)
        pp_due_over_nue_dxi             % convective ratio derivative interpolant (C0 continuity, piecewise C1)        
        
        % Model Options (as in Rfoil)
        scci = 5.60                     % Shear lag constant (add gacon and gbcon in the future!)
        
        % Initial Conditions for Integration
        t0
        h0
        ctau0
        
        % Integration Results
        sol
        
        % Integration algorithm parameters
        ode_solver  = @ode15s;                          % ode45 or ode15s (stiff) are reasonable choices
        fixed_step  = true;                             % Fixed or variable lenght ODE intgration steps
        N_steps     = 1000;                             % Set number of steps (only applicable for fixed_step mode
    end
    
    methods
        function BL = boundary_layer6()
            % Constructor Function
        end
        
        function set_forcing_terms(BL, xi_vector, ue_vector, msq_vector)
            % Store Supplied Vectors
            BL.xi_vector     = xi_vector;
            BL.ue_vector     = ue_vector;
            BL.msq_vector    = msq_vector;
            
            % And Update Interpolants and Direct Dependencies
            BL.update_interpolants()
        end
        
        function update_interpolants(BL)
            % Generate new ue_over_nue vector
            BL.ue_over_nue_vector   = BL.ue_vector .* BL.u_inf_over_nu_inf ;
            
            % Make interpolants
            BL.pp_ue                = pchip(BL.xi_vector, BL.ue_vector);
            BL.pp_due_dxi           = ppdiff(BL.pp_ue,1);
            
            BL.pp_msq               = pchip(BL.xi_vector, BL.msq_vector);
            BL.pp_dmsq_dxi          = ppdiff(BL.pp_msq,1);
            
            BL.pp_ue_over_nue       = pchip(BL.xi_vector, BL.ue_over_nue_vector);
            BL.pp_due_over_nue_dxi  = ppdiff(BL.pp_ue_over_nue,1);
        end
        
        function ue = ue_function(BL, x)
            ue                  = ppval(BL.pp_ue, x);
        end
        
        function due_dxi = due_dxi_function(BL, x)
            due_dxi             = ppval(BL.pp_due_dxi, x);
        end
        
        function msq = msq_function(BL, x)
            msq                 = ppval(BL.pp_msq, x);
        end
                      
        function dmsq_dxi = dmsq_dxi_function(BL, x)
            dmsq_dxi            = ppval(BL.pp_dmsq_dxi, x);
        end
        
        function ue_over_nue = ue_over_nue_function(BL, x)
            ue_over_nue             = ppval(BL.pp_ue_over_nue, x);
        end                
        
        function due_over_nue_dxi = due_over_nue_dxi_function(BL, x)
            due_over_nue_dxi        = ppval(BL.pp_due_over_nue_dxi, x);
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
            rhs_hs =  2 * cdi ./ t  - (hs ./ t) .* (cf / 2) - (2 * hc + hs.*(1 - h)) .* (1 ./ ue) .* due_dxi;
            
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
        
        function [rhs_theta, rhs_h, rhs_ctau] = transformed_rhs(BL, t, h, ctau, x)
            % % Get canonical rhs
            [rhs_theta_c, rhs_hs_c, rhs_ctau_c] = canonical_rhs(BL, t, h, ctau, x);
            % % Transform rhs
            [rhs_theta, rhs_h, rhs_ctau] = transform_rhs(BL, t, h, ctau, x, rhs_theta_c, rhs_hs_c, rhs_ctau_c);
            % Done! Return!
        end
        
        function [rhs_theta, rhs_h, rhs_ctau] = transform_rhs(BL, t, h, ctau, x, rhs_theta, rhs_hs, rhs_ctau)%#ok<INUSL>
            % Get interpolated forcing terms
            ue                  = BL.ue_function(x);
            due_dxi             = BL.due_dxi_function(x);
            msq                 = BL.msq_function(x);
            dmsq_dxi            = BL.dmsq_dxi_function(x);
            ue_over_nue         = BL.ue_over_nue_function(x);
            due_over_nue_dxi    = BL.due_over_nue_dxi_function(x);            
            
            % % Get preliminary variables (inneficient and repeated code, but clearer in terms of namespace)
            % Shape Factor
%            [ h, h_dstr, h_t] = H(dstr , t);                            % Looks ok
            dstr = h * t;

            % Kinematic Shape Factor
            [ hk, hk_h, hk_msq  ] = hkin( h, msq);                      % Looks ok            

            % Reynolds Theta
            [rt, rt_ue_over_nue ] = re_theta( t, ue_over_nue );         % Looks ok
            
            % Now energy shape factor (H_star)
            [ hs, hs_hk, hs_rt, hs_msq  ] = hst( hk, rt, msq);          % Check not possible
                                   
            % % Now start writting transformation terms (new version in boundary_layer3, based on alfa9 revision of pdf)
            gamma               =   hs_hk .* hk_h;                      % Denominator
            lambda_hs           =   1;                                  % hs          specific sensitivity term (kinetic energy shape factor)
            lambda_t            = - ue_over_nue .* hs_rt;               % t           specific sensitivity term (momentum thickness)
            lambda_ue_over_nue  = - t .* hs_rt;                         % ue_over_nue specific sensitivity term (ratio of edge speed to edge viscosity)
            lambda_msq          = - (hs_hk .* hk_msq + hs_msq);         % msq         specific sensitivity term (edge bmach number)                                 
            
            % And transform!
            % ERROR ERROR error % 1/gamma % 
            rhs_h = gamma .* (  lambda_t           .* rhs_theta        + ...
                                lambda_hs          .* rhs_hs           + ...
                                lambda_ue_over_nue .* due_over_nue_dxi + ...
                                lambda_msq         .* dmsq_dxi);
        end
                
        function [y_line] = ode_fun_wrapper(BL, x, y)
            % Wrapper function of BL.transformed_rhs method for matlab
            % standard ode fun interface
            
            t       = y(1);
            h       = y(2);
            ctau    = y(3);
            
            [rhs_theta, rhs_h, rhs_ctau] = transformed_rhs(BL, t, h, ctau, x);
            
            y_line = zeros(length(y),1);
            y_line(1) = rhs_theta;
            y_line(2) = rhs_h;
            y_line(3) = rhs_ctau;
        end
        
        function [d_theta_dx, d_h_dx, d_ctau_dx] = manual_fun_wrapper(BL, t, h, ctau, x)
            % Wrapper function of BL.transformed_rhs method for manual ode
            % integration
            
            [rhs_theta, rhs_h, rhs_ctau] = transformed_rhs(BL, t, h, ctau, x);
            
            d_theta_dx  = rhs_theta;
            d_h_dx      = rhs_h;
            d_ctau_dx   = rhs_ctau;
            
        end
        
        function set_explicit_initial_conditions(BL, t0, h0, ctau0)
            % Simply Store Explicit Initial Conditions
            BL.t0       = t0;               % Adimensionalized Momentum Thickness (t = delta_2 / L)
            BL.h0       = h0;               % Shape Factor
            BL.ctau0    = ctau0;            % Shear Stress Coefficient (includes Turbulent component, for max shear point in BL)
        end
        
        function set_equilibrium_initial_conditions(BL, t0, h0)
            % Store Explicit Initial Conditions
            BL.t0       = t0;               % Adimensionalized Momentum Thickness
            BL.h0       = h0;               % Adimensionalized Shape Factor
            
            % % Guess corresponding ctau (assume equilibrium ctau)
            % From equilibrium ctau (even though even a flat plate is not in perfect equilibrium!)            
            msq0    = BL.msq_vector(1);
            hk      = hkin( h0, msq0);
            rt      = re_theta( t0, BL.u_inf_over_nu_inf );
            cf      = cft( hk, rt, msq0);
            hs      = hst( hk, rt, msq0);
            us      = usg( h0, hk, hs);
            ctau_eq = ctauzero( h0 , hk , hs , us);
            
            % Set Initial Ctau as Equilibrium Ctau            
            BL.ctau0 = ctau_eq;
        end
        
        function y0 = initial_conditions_vector(BL)
            % % Generate initial condition vector for integrator!
            y0    = zeros(1,3);
            
            y0(1) = BL.t0;
            y0(2) = BL.h0;
            y0(3) = BL.ctau0;
        end
        
        function integrate_ODE(BL, x0, xend)
            % Integrate 
            
            % Generate standard format initial condition vector
            y0 = initial_conditions_vector(BL);
            
            % Make suitable handle function reference to wrapper
            odefun = @BL.ode_fun_wrapper;
                        
            % Integrate
            if BL.fixed_step == false
                % Variable Step Size
                BL.sol = BL.ode_solver(odefun,[x0 xend],y0);
            else
                % Fixed Step Size
                BL.sol = BL.ode_solver(odefun,linspace(x0, xend, BL.N_steps),y0);
            end
            
        end
        
    end
    
end

