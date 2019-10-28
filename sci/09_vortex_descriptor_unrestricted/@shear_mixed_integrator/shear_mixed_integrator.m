classdef shear_mixed_integrator < handle
    %SHEAR_MIXED_INTEGRATOR is a handle class that generates integral 
    % boundary layer closure parameters for the total field by combining 
    % the pure shear field with the mixed field
    
    properties
        CM      % CrossflowMesh object
        SF      % ShearField    object
        MF      % MixedField    object
    end
    
    methods
        % Constructor Method
        function SMI = shear_mixed_integrator(CM, SF, MF)
            % Store Inputs
            SMI.CM = CM;
            SMI.SF = SF;
            SMI.MF = MF;
        end
        
        function [dstr_bar, dstr_tilde, dstr_hat] = dstr_direct(SMI, hk_bar, rt_bar, ue, nue, u_tilde)
            % Get ue_over_nue
            ue_over_nue = ue / nue;
            
            % Get u_bar from pure shear object
            u_bar = SMI.SF.u_bar_over_mesh(hk_bar, rt_bar, ue, ue_over_nue);
            
            % Make u_bar over ue and u_tilde over ue
            u_bar_over_ue   = u_bar   / ue;
            u_tilde_over_ue = u_tilde / ue;
            u_hat_over_ue   = u_bar_over_ue + u_tilde_over_ue;
            
            % Make integrands (according to thicknesses lyx document)
            integrand_bar   = 1 - u_bar_over_ue ;
            integrand_hat   = 1 - u_hat_over_ue ;

            % Integrate!
            dstr_bar_int    = trapz(SMI.CM.y_range,trapz(SMI.CM.z_range,integrand_bar ,2),1) / SMI.CM.z_lenght;
            dstr_hat_int    = trapz(SMI.CM.y_range,trapz(SMI.CM.z_range,integrand_hat ,2),1) / SMI.CM.z_lenght;
            
            % Compute hat value (according to new approach)
            dstr_tilde        = dstr_hat_int - dstr_bar_int;
            
            % Reconstruct bar value from reference
            theta_bar   = rt_bar / ue_over_nue;
            dstr_bar        = hk_bar * theta_bar;
            
            % Reconstruct final hat value using reference bar value (to circumvent numerical integration uncertainties)
            dstr_hat        = dstr_bar + dstr_tilde;
        end
        
        function [theta_bar, theta_tilde, theta_hat] = theta_direct(SMI, hk_bar, rt_bar, ue, nue, u_tilde)
            % Get ue_over_nue
            ue_over_nue = ue / nue;
            
            % Get u_bar from pure shear object
            u_bar = SMI.SF.u_bar_over_mesh(hk_bar, rt_bar, ue, ue_over_nue);
            
            % Make u_bar over ue and u_tilde over ue
            u_bar_over_ue   = u_bar   / ue;
            u_tilde_over_ue = u_tilde / ue;
            u_hat_over_ue   = u_bar_over_ue + u_tilde_over_ue;
            
            % Make integrands (according to new approach)
            integrand_bar   = u_bar_over_ue .* (1 - u_bar_over_ue) ;
            integrand_hat   = u_hat_over_ue .* (1 - u_hat_over_ue) ;

            % Integrate!
            theta_bar_int    = trapz(SMI.CM.y_range,trapz(SMI.CM.z_range,integrand_bar ,2),1) / SMI.CM.z_lenght;
            theta_hat_int    = trapz(SMI.CM.y_range,trapz(SMI.CM.z_range,integrand_hat ,2),1) / SMI.CM.z_lenght;
            
            % Compute hat value (according to thicknesses lyx document)
            theta_tilde        = theta_hat_int - theta_bar_int;
            
            % Reconstruct bar value from reference
            theta_bar   = rt_bar / ue_over_nue;
            
            % Reconstruct final hat value using reference bar value (to circumvent numerical integration uncertainties)
            theta_hat        = theta_bar + theta_tilde;
        end
        
        function [d3_bar, d3_tilde, d3_hat] = delta3_direct(SMI, hk_bar, rt_bar, ue, nue, u_tilde)
            % Get ue_over_nue
            ue_over_nue = ue / nue;
            
            % Get u_bar from pure shear object
            u_bar = SMI.SF.u_bar_over_mesh(hk_bar, rt_bar, ue, ue_over_nue);
            
            % Make u_bar over ue and u_tilde over ue
            u_bar_over_ue   = u_bar   / ue;
            u_tilde_over_ue = u_tilde / ue;
            u_hat_over_ue   = u_bar_over_ue + u_tilde_over_ue;
            
            % Make integrands (according to new approach)
            integrand_bar   = u_bar_over_ue .* (1 - u_bar_over_ue.^2) ;
            integrand_hat   = u_hat_over_ue .* (1 - u_hat_over_ue.^2) ;

            % Integrate!
            d3_bar_int    = trapz(SMI.CM.y_range,trapz(SMI.CM.z_range,integrand_bar ,2),1) / SMI.CM.z_lenght;
            d3_hat_int    = trapz(SMI.CM.y_range,trapz(SMI.CM.z_range,integrand_hat ,2),1) / SMI.CM.z_lenght;
            
            % Compute hat value (according to thicknesses lyx document)
            d3_tilde        = d3_hat_int - d3_bar_int;
            
            % Reconstruct bar value from reference
            theta_bar   = rt_bar / ue_over_nue;                 % Momentum thickness from Re_theta
            msq         = 0;                                    % Force to incompressible!
            hst_bar     = hst( hk_bar, rt_bar, msq);            % Energy Shape Factor from selected closure! (right now IVW2)
            d3_bar      = theta_bar * hst_bar;
            
            % Reconstruct final hat value using reference bar value (to circumvent numerical integration uncertainties)
            d3_hat       = d3_bar + d3_tilde;
        end
        
        function [hk_bar, hk_tilde, hk_hat] = hk_direct(SMI, hk_bar, rt_bar, ue, nue, u_tilde)
            
        end
        
        function [hstr_bar, hstr_tilde, hstr_hat] = hstr_direct(SMI, hk_bar, rt_bar, ue, nue, u_tilde)
            
        end
        
        function [cd_bar, cd_tilde, cd_hat] = cd_direct(SMI, hk_bar, rt_bar, ue, nue, u_tilde)
            
        end
        
        
    end
    
end

% % A note on the integration method!
% Mathworks Example : blogs.mathworks.com/loren/2011/06/13/calculating-the-area-under-a-surface/
%            [X,Y] = meshgrid(x,y);
%            Z = X.^2.*sin(3*(X-Y));
%            trapz(y,trapz(x,Z,2),1);
%
% Recall that:
%    [CM.z_mesh, CM.y_mesh] = meshgrid(CM.z_range, CM.y_range);
%
% To establish the equivalence:
%    x === CM.z_range
%    y === CM.y_range
%    Z === integrand_bar (or tilde!)
%    2 === array direction!
%    1 === array direction!
%
% And write
%    integral_result = trapz(CM.y_range,trapz(CM.z_range,integrand,2),1);

