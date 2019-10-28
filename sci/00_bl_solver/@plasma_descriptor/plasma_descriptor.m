classdef plasma_descriptor < handle
    %PLASMA_DESCRIPTOR is a class used to store and manage the description
    %of a plasma actuator. It stores both the plasma field function, and
    %the plasma properties.
    
    properties        
        % Dimensional Descriptors, provided as INPUTS
        F_T                     % [N/m] Total Actuator Force per Unit Span        
        L_p                     % [m]   Plasma Force Field Lenght
        T_p                     % [m]   Plasma Force Field Thickness
        X0_p                    % [m]   Plasma Force Field Start
        
        % Dimensional Descriptors, DEPENDENT
        phi_x_p                 % [--]  Average Force Field Density
        
        % Adimensional Descriptors, DEPENDENT
        t_p                     % [--]  Adimensional Force Field Thickness        
        l_p                     % [--]  Adimensional Force Field Lenght
        x0_p                    % [--]  Adimensional Force Field Start

        c_phi_x_p               % [--]  Average Force Field Density
        
        % t_bar_p               % [--]  Scaled Force Field Thickness
        % t_theta_p             % [--]  Momentum Scaled Force Field Thickness
                        
        % Adimensionalization Parameters, Free Stream, INDEPENDENT
        L       = 1             % [m]   Lenght Scale (Longitudinal)
        U_inf   = 1             % [m/s] Unperturbed Flow Speed
        rho_inf = 1             % [kg/m3] Unperturbed Flow Density
        
        % Units are provided as an indication, assuming that:
        %       L       is in [m]
        %       U_inf   is in [m/s]
        %       rho_ing is in [kg/m3]
    end
    
    methods
        function PD = plasma_descriptor()
           % Constructor class 
        end
        
        function set_plasma_dimensional_inputs(PD, F_T, L_p, T_p, X0_p)
            % Store Inputs
            PD.F_T  = F_T;
            PD.L_p  = L_p;
            PD.T_p  = T_p;
            PD.X0_p = X0_p;
            
            % Update dependent parameters
            PD.update_dependent_parameters();
        end
        
        function set_adimensionalization_parameters(PD, L, U_inf, rho_inf)
           % Store Inputs
            PD.L        = L;
            PD.U_inf    = U_inf;
            PD.rho_inf  = rho_inf;
            
            % Update dependent parameters
            PD.update_dependent_parameters(); 
        end
        
        function update_dependent_parameters(PD)
            % Make scaled lenght and thickness
            PD.t_p      = PD.T_p  / PD.L;
            PD.l_p      = PD.L_p  / PD.L;
            PD.x0_p     = PD.X0_p / PD.L;
            % Compute Average Force Field
            PD.phi_x_p  = PD.F_T / (PD.L_p * PD.T_p);
            
            % Make Average Force Field Coefficient
            PD.c_phi_x_p = (PD.phi_x_p *  PD.T_p) / (0.5 * PD.rho_inf * PD.U_inf^2);                                   
        end
        
        function wy = wy_function_generic(PD, y, b)
            % wY weighting function, as described in report.
            % Supports vectorized evaluation (conditionality handled as
            % multiplication!)
            wy = 0.5*pi() * sin( 0.5*pi()*(y/b + 1)) .* (y/b<1) .* (y>0);
        end
        
        function wx = wx_function_generic(PD, x, x0, a)
            % wX weighting function, as described in report.
            % Supports vectorized evaluation (conditionality handled as
            % multiplication!)
            wx = 0.5 * pi() * sin(pi() * (x-x0)./a) .* (((x-x0)./a) < 1) .* ((x-x0)>0);                        
        end
        
        function wx = wx_function_CFM(PD, x)
            % wX weighting function with arguments set for CFM evaluation
            % from C_phi_x
            wx = PD.wx_function_generic(x, PD.x0_p, PD.l_p);                        
        end
        
    end
    
end

