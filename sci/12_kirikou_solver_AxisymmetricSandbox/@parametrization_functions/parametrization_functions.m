classdef parametrization_functions
    %PARAMETRIZATION_FUNCTIONS is a simple class collecting static wake
    % parametrization functions for the Dogoro and Kirikou solvers
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                                                     %
    %   Part of -   A simple specialized 2d Vorticity Equation Solver for %
    %               Actuator Disk Flows (Kirikou-Dogoro Suite)            %
    %                                                                     %
    %   Date    :   June 2014 to March 2017                               %
    %   Author  :   Gael de Oliveira                                      %
    %                                                                     %
    %   License :   Case by case written agreement limited to specific    %
    %               applications. Distribution to any individual or       %
    %               organization requires explicit written agreement from %
    %               original author.                                      %
    %                                                                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
    end
    
    methods
    end
    
    methods(Static)
        function [h , dh_dbeta] = h_streamline(theta , beta_vector)
            % Returns streamline height function in (h, theta) space
            %    (Vectorized in theta, also returns sensitivities of h to beta)
            
            % Prepare (find degree of polynomial basis)
            P = length(beta_vector)-1;
            
            % Request bernstein polynomial basis vector
            b_vector = parametrization_functions.modified_bernstein_polynomial_basis_vector(theta,P);
            
            % Compute product with coefficients to obtain streamline h at demanded
            % theta
            h =zeros(size(theta(:)));
            for n_beta=1:(P+1)
                h = h + beta_vector(n_beta) * b_vector(:, n_beta);
            end
            
            % h sentitivities to beta coefficients (each row corresponds to a theta
            % value, each column to a beta coefficient)
            dh_dbeta = b_vector;
        end
        
        function b = modified_bernstein_polynomial(theta,P,k)
            % Degree P, order P+1, number k
            tt_o_p = (2/pi) * theta;                            % Two theta over pi
            b = nchoosek(P,k) * tt_o_p.^(k) .* (1-tt_o_p).^(P-k);   % Modified Bernstein Polinomial
        end
        
        function b_vector = modified_bernstein_polynomial_basis_vector(theta,P)
            % Degree P, order P+1, Complete Basis in a Vector (k=0...P)
            
            % Prepare
            b_vector = zeros(length(theta), P+1);                            % Allocate Basis Vector
            tt_o_p = (2/pi) * theta(:);                            % Two theta over pi
            
            % Compute
            for k = 0:P
                b_vector(: , k+1) = nchoosek(P,k) * tt_o_p.^(k) .* (1-tt_o_p).^(P-k);   % Modified Bernstein Polinomial
            end
        end
        
        function [x,y, J11, J12, J21, J22] = S_transformation(h, theta, xs, ys)
            % Transformation from polar shooting coordinate system to cartesian system
            %   Also returns Jabobian:
            %       J11 = dx_dh
            %       J12 = dx_dtheta
            %       J21 = dy_dh
            %       J22 = dy_dtheta
            %
            %   Fully vectorized for h-theta pairs and single origin (xs, ys)
            
            
            
            % Simply write the expression fromt the report!
            x = xs + h.*tan(theta);
            y = ys + h;
            
            % And also compute jacobian (unrolled for vectorization)
            J11 = tan(theta);                   % dx_dh
            J12 = h .* (1+(tan(theta)).^2);     % dx_dtheta
            J21 = ones(size(theta));            % dy_dh
            J22 = zeros(size(theta));           % dy_dtheta
            
        end
        
    end
    
end
