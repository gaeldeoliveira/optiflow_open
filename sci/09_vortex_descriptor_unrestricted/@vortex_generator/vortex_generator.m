classdef vortex_generator < handle
    %VORTEX_GENERATOR is a handle class that stores all the information
    % pertaining to the geometric description of a vortex generator. Also
    % implements the stripline and Wendt2001 models for initial
    % circulation and initial peak vorticity.
    
    properties
        d_VG                               % [m    ] distance between VG vanes in pair (not distance between pairs)
        c_VG                               % [m    ] chord of VG  (5cm)
        h_VG                               % [m    ] height of VG
        AoA_VG                             % [deg  ] VG angle of attack (geometric, to free-stream)
        
        sigma_rel = 0.2;                   % [adim.] relative core size (core size over height) 
                                           % (Warning: this is an arbitrary, gut-feeling choice, to get some insight in the short term, but only for that!)

        % Sigma has two meanings:
        %    core size when using rankine vortex induction (stripline model)
        %    omega_max when using lamb    vortex induction (wendt2001 model)
    end
    
    methods
        % Constructor Method
        function VG = vortex_generator(d_VG, c_VG, h_VG, AoA_VG)
            VG.d_VG   = d_VG;
            VG.c_VG   = c_VG;
            VG.h_VG   = h_VG;
            VG.AoA_VG = AoA_VG;
        end
        
        function [y_v0, z_v0, gamma_v0, sigma_v0] = initial_strenght_stripline_model(VG, U_ref)
            % Returns initial vortex filament strenght (gamma_v,
            % circulation per unit lenght) based on 2d flat plate analogy.
            % Also returns an unfadamented guess for sigma_v, the core
            % radius.
            
            % Compute CL based on 2d flat plate expression!
            Cl_VG  = 2*pi()*sin(VG.AoA_VG*pi()/180);
            
            % Return Initial Vortex Descriptors!
            y_v0    = VG.h_VG;                      % [m    ] initial filament height
            z_v0    = VG.d_VG / 2;                  % [m    ] initial filament distance to primary symmetry line (the one that exists with a single vortex pair!) (half of distance between trailing edges of vortex generator pair!)
            gamma_v0 = -0.5*U_ref*VG.c_VG*Cl_VG;    % [m /s ] circulation of VG
            sigma_v0 = VG.sigma_rel*VG.h_VG;        % [m    ] core size of VG
        end
        
        function [y_v0, z_v0, gamma_v0, sigma_v0] = initial_strenght_prandtl_model(VG, U_ref)
            % Returns initial vortex filament strenght (gamma_v,
            % circulation per unit lenght) based on 2d flat plate analogy.
            % Also returns an unfadamented guess for sigma_v, the core
            % radius.
            
            % Compute CL based on 2d flat plate expression!
            Cl_VG  = 2*pi()*sin(VG.AoA_VG*pi()/180);
            
            % Make modified aspect ratio
            AR = (8 * VG.h_VG) / (pi() * VG.c_VG);
            
            % Apply Lifting line correction
            Cl_VG_3d= Cl_VG / (1 + 2/AR);
            
            % Return Initial Vortex Descriptors!
            y_v0    = VG.h_VG;                      % [m    ] initial filament height
            z_v0    = VG.d_VG / 2;                  % [m    ] initial filament distance to primary symmetry line (the one that exists with a single vortex pair!) (half of distance between trailing edges of vortex generator pair!)
            gamma_v0 = -0.5*U_ref*VG.c_VG*Cl_VG_3d; % [m /s ] circulation of VG
            sigma_v0 = VG.sigma_rel*VG.h_VG;        % [m    ] core size of VG
        end
        
        function [y_v0, z_v0, gamma_v0, sigma_v0] = initial_strenght_wendt2001_model(VG, U_ref, delta_ref)
            
            % Implement Wendt2001 model (NASA CR-2001-211144)
            k1 = 1.61;
            k2 = 0.48;
            k3 = 1.41;
            k4 = 1.00;
            xi = 0.29;
            
            % Make modified aspect ratio
            AR = (8 * VG.h_VG) / (pi() * VG.c_VG);
            
            % Get radian angle of attack
            alfa = - VG.AoA_VG * pi() / 180;
            
            % Compute Circulation
            gamma = (k1 * alfa * U_ref * VG.c_VG) ...
                          / (1 + k2 / AR) ...
                                    *  tanh(k3*(VG.h_VG / delta_ref)^k4);
                                
            % Compute beta parameter (for peak vorticity)
            Beta = 1 / (2 * xi.^2 * ( 1-exp(-0.5) )^2 );
            
            % Compute peak vorticity
            omega_max = ( gamma^3 * ( Beta-1 )^2 ) / ...
                            (2 * pi^3 * (alfa)^2 * (VG.c_VG)^2 * (VG.h_VG)^2 * U_ref^2 );
            
            % Return Initial Vortex Descriptors!
            y_v0     = VG.h_VG;                      % [m    ] initial filament height
            z_v0     = VG.d_VG / 2;                  % [m    ] initial filament distance to primary symmetry line (the one that exists with a single vortex pair!) (half of distance between trailing edges of vortex generator pair!)
            gamma_v0 = gamma;                        % [m /s ] circulation of VG
            sigma_v0 = omega_max;                    % [m    ] core size of VG
        end
        
        function [t_v0] = initial_vortex_wendt1995_model(VG, gamma_v0, sigma_v0, nu_v) %#ok<INUSL>
            % Implement the expression for finding the "initial" age of the
            % vortex from its peak vorticity at the trailing edge of the
            % vane tip
            %
            % Admits scalar arguments and/or nD arrays of equal size:
            %       gamma_v0 - circulation at point of interest (vane TE)
            %       sigma_v0 - peak vorticity at point of interest (vane TE)
            %       nu_v     - kinematic viscosity (may differ from molecular nu to incorporate effect of turbulence!)
            %
            
            t_v0 = 1 / (4*pi()) * gamma_v0 ./ (nu_v .* sigma_v0);
            
        end
        
        
        
        
    end
    
end

