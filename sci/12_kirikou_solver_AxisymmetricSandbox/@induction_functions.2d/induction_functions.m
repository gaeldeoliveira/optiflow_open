classdef induction_functions
    %INDUCTION_FUNCTIONS is a simple class collecting static induction
    % functions for the Dogoro and Kirikou solvers
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
        function [u , v] = constant_strenght_singular_vortex_pair(x, y, x_v, y_v, gamma_v)
            % Two Dimensional Constant Strenght Singular Vortex Pair Induction Function
            %
            % Inputs:
            %    x  ,   y      % Target
            %    x_v,   y_v    % Upper Halfplane Singular Vortex Center
            %    x_v, - y_v    % Lower Halftplane Singular Vortex Center
            %    gamma_v       % Strenght of Upper Singular Vortex
            %  - gamma_v       % Strenght of Lower Singular Vortex
            
            %
            % Outputs:
            %    phi       % Potential (currently zero)
            %    u         % x speed component
            %    v         % y speed component
            %
            %
            % Reference:
            %    Green Function Solution of the Laplace Equation
            %
            % Validation:
            %    For Validation of lifting vortex functions see lifting_vortex_test.m
            %            Katz and Plotkin, Chapter 10, Section 3 (p. 236)
            %
            %
            
            %\ % Lifting Vortex Helper Functions and Variables
            r_sq_v_fun = @(x,y,x_v,y_v) (x-x_v).^2 + (y-y_v).^2;
            u_v_fun = @(x,y,x_v,y_v, gamma_v) gamma_v / (2*pi()) * (-(y-y_v)) ./ r_sq_v_fun(x,y,x_v,y_v);
            v_v_fun = @(x,y,x_v,y_v, gamma_v) gamma_v / (2*pi()) * (  x-x_v ) ./ r_sq_v_fun(x,y,x_v,y_v);
            
            % Now proceed to define induction function for our lifting vortex pair
            % Upper side plus lower side vortices (inverted circulation and y!)
            % (both pull air in/out)!
            u_v = @(x,y) u_v_fun(x,y,x_v,y_v, gamma_v) + u_v_fun(x,y,x_v,-y_v, -gamma_v);
            v_v = @(x,y) v_v_fun(x,y,x_v,y_v, gamma_v) + v_v_fun(x,y,x_v,-y_v, -gamma_v);
            %                       upper^  +gamma_v^                lower^  -gamma_v^
            
            % Now compute
            u = u_v(x, y);
            v = v_v(x, y);
            
            % Cancel out self-induction on center (avoid NaN and sort our issues!)
            %%% Make consistency checks and corrections
            precision_threshold = 1e-15;
            % Condition for point to be considered on center of upper
            condition_on_center_u = r_sq_v_fun(x,y,x_v, y_v) < precision_threshold;
            % Cancel out self-induction on center of upper (avoid NaN and sort our issues!)
            u(condition_on_center_u) = 0;
            v(condition_on_center_u) = 0;
            % Condition for point to be considered on center of lower
            condition_on_center_l = r_sq_v_fun(x,y,x_v,-y_v) < precision_threshold;
            % Cancel out self-induction on center of upper (avoid NaN and sort our issues!)
            u(condition_on_center_l) = 0;
            v(condition_on_center_l) = 0;
            
        end
        
        function [potential , u , v] = constant_strenght_vortex_segment_element(x, y, x1, y1, x2, y2, gamma)
            % Two Dimensional Constant Strenght Singularity Segment Elements
            %
            % Inputs:
            %    x , y     % Target
            %    x1, y1    % Source Start
            %    x2, y2    % Source End
            %
            % Outputs:
            %    phi       % Potential
            %    u         % x speed component
            %    v         % y speed component
            %
            %
            % Reference:
            %    Based on:
            %            Katz and Plotkin, Chapter 10, Section 3 (p. 236)
            %    Represents analytical integral of:
            %             psi = - gamma / 2*pi * int(atan(y/(x-x0), x0, x1, x2);
            %    Generalized by embedded the coordinate transformation so that
            %    arbitrary y, y1 and y2 are allowed
            %
            
            % Squared distances from target to source points
            r1_sq     = (x-x1).^2 + (y-y1).^2 + eps() ;
            r2_sq     = (x-x2).^2 + (y-y2).^2 + eps() ;
            
            % Projected Distances on Element Axis
            l_sq = (x2-x1).^2 + (y2-y1).^2;                 % Square of vortex segment lenght
            l    = sqrt(l_sq);                              % Vortex Segment lenght (always positive)
            
            l1   = ( (x-x1) .* (x2-x1) + (y-y1) .* (y2-y1) ) ./l;   % Distance to from target to start along element axis (can be positive or negative)
            l2   = ( (x-x2) .* (x2-x1) + (y-y2) .* (y2-y1) ) ./l;   % Distance to from target to end   along element axis (can be positive or negative)
            
            % Projected Distances on Normal Axis
            % Simple normal expression for 2d (easier than subtracting colinear
            % component, which is generic for n-D)
            n_versor_x = - (y2-y1)  ./ l;  % x component of normal unit vector
            n_versor_y =   (x2-x1)  ./ l;  % y component of normal unit vector
            
            % Normal distance
            n = (x-x1) .* n_versor_x + (y-y1) .* n_versor_y;
            % Altnernative normal distance
            % n = (x-x2) .* n_versor_x + (y-y2) .* n_versor_y
            
            % Angles from target to source points, as measured from the element axis
            % theta1    = atan(n ./ l1);
            % theta2    = atan(n ./ l2);
            
            theta1    = atan2(n , l1);
            theta2    = atan2(n , l2);
            
            % Potential (generalized for arbitraty y1 and y2, I imagine he thinks it is easy to find!)
            potential = - gamma ./ (2*pi()) .* (...
                l1 .* theta1 ...
                - l2 .* theta2 ...
                + 0.5 .* n  .* log(r1_sq./r2_sq));
            
            % Segment Axis Speed component
            s         =   gamma ./ (2*pi()) .* (theta2 - theta1);
            % Normal Axis speed component
            t         = - gamma ./ (4*pi()) .* log(r1_sq ./ r2_sq);      % Minus sign not present in Katz and Plotkin, but apparently needed!
            
            %%% Make consistency checks and corrections
            precision_threshold = 1e-15;
            % Condition for normal distance to be considered as tending to zero
            condition_on_segment = abs(n)<precision_threshold;
            s(condition_on_segment) = 0;
            if max(size(gamma)) == max(size(condition_on_segment))
                t(condition_on_segment) = - gamma(condition_on_segment) ./ (4*pi()) .* log(l1(condition_on_segment).^2 ./ l2(condition_on_segment).^2);
            else
                t(condition_on_segment) = - gamma ./ (4*pi()) .* log(l1(condition_on_segment).^2 ./ l2(condition_on_segment).^2);
            end
            
            %%% Transform (s,t) speed into (u,v) speed
            % Define l versor
            l_versor_x =  n_versor_y; %  (x2-x1)  ./ l;  % x component of tangential unit vector
            l_versor_y = -n_versor_x; %  (y2-y1)  ./ l;  % y component of tangential unit vector
            
            % Compute speed vector components (combine l,n system basis with s,t coefficients)
            u = s .* l_versor_x + t .* n_versor_x;
            v = s .* l_versor_y + t .* n_versor_y;
            
        end
        
        function [potential , u , v, u_dx, u_dy, v_dx, v_dy] = constant_strenght_vortex_segment_element_derivatives(x, y, x1, y1, x2, y2, gamma)
            % Two Dimensional Constant Strenght Singularity Elements
            %
            % Inputs:
            %    x , y     % Target
            %    x1, y1    % Source Start
            %    x2, y2    % Source End
            %
            % Outputs:
            %    phi       % Potential
            %    u         % x speed component
            %    v         % y speed component
            %    u_dx      % x derivative of x speed component
            %    u_dy      % y derivative of x speed component
            %    v_dx      % x derivative of y speed component
            %    v_dy      % y derivative of y speed component
            %
            %    Derivatives are fudged to 0 on segment (useful for my purposes but not
            %    strictly correct!)
            %
            %
            % Reference:
            %    Based on:
            %            Katz and Plotkin, Chapter 10, Section 3 (p. 236)
            %    Represents analytical integral of:
            %             psi = - gamma / 2*pi * int(atan(y/(x-x0), x0, x1, x2);
            %    Generalized by embedded the coordinate transformation so that
            %    arbitrary y, y1 and y2 are allowed
            %
            
            % Squared distances from target to source points
            r1_sq      = (x-x1).^2 + (y-y1).^2 + eps() ;
            r2_sq      = (x-x2).^2 + (y-y2).^2 + eps() ;
            % Make derivatives to x,y
            r1_sq_dx   = 2*(x-x1)              + eps() ;
            r1_sq_dy   =              2*(y-y1) + eps() ;
            
            r2_sq_dx   = 2*(x-x2)              + eps() ;
            r2_sq_dy   =              2*(y-y2) + eps() ;
            
            % Projected Distances on Element Axis
            l_sq = (x2-x1).^2 + (y2-y1).^2;                 % Square of vortex segment lenght
            l    = sqrt(l_sq);                              % Vortex Segment lenght (always positive)
            
            l1    = ( (x-x1) .* (x2-x1) + (y-y1) .* (y2-y1) ) ./l;   % Distance to from target to start along element axis (can be positive or negative)
            l2    = ( (x-x2) .* (x2-x1) + (y-y2) .* (y2-y1) ) ./l;   % Distance to from target to end   along element axis (can be positive or negative)
            % Make derivatives to x,y
            l1_dx = ( (1   ) .* (x2-x1)                     ) ./l;
            l1_dy = (                     (1   ) .* (y2-y1) ) ./l;
            l2_dx = ( (1   ) .* (x2-x1)                     ) ./l;
            l2_dy = (                     (1   ) .* (y2-y1) ) ./l;
            
            
            
            
            
            
            % Projected Distances on Normal Axis
            % Simple normal expression for 2d (easier than subtracting colinear
            % component, which is generic for n-D)
            n_versor_x = - (y2-y1)  ./ l;  % x component of normal unit vector
            n_versor_y =   (x2-x1)  ./ l;  % y component of normal unit vector
            
            % Normal distance
            n    = (x-x1) .* n_versor_x + (y-y1) .* n_versor_y;
            % Derivatives to x-y
            n_dx = (1   ) .* n_versor_x;
            n_dy =                        (1   ) .* n_versor_y;
            
            % Altnernative normal distance
            % n = (x-x2) .* n_versor_x + (y-y2) .* n_versor_y
            
            % Angles from target to source points, as measured from the element axis
            % theta1    = atan(n ./ l1);
            % theta2    = atan(n ./ l2);
            
            theta1    = atan2(n , l1);
            theta2    = atan2(n , l2);
            
            % Derivatives (via symbolic toolbox)
            theta1_dl1   = -n ./ ((n.^2      + l1.^2));
            theta1_dn    =  1 ./ ((n.^2 ./l1 + l1  ));
            
            theta2_dl2   = -n ./ ((n.^2      + l2.^2));
            theta2_dn    =  1 ./ ((n.^2 ./l2 + l2  ));
            
            % Combine with chain rule
            theta1_dx   = theta1_dl1 .* l1_dx + theta1_dn .* n_dx;
            theta1_dy   = theta1_dl1 .* l1_dy + theta1_dn .* n_dy;
            
            theta2_dx   = theta2_dl2 .* l2_dx + theta2_dn .* n_dx;
            theta2_dy   = theta2_dl2 .* l2_dy + theta2_dn .* n_dy;
            
            
            
            
            % Potential (generalized for arbitraty y1 and y2, I imagine he thinks it is easy to find!)
            potential = - gamma ./ (2*pi()) .* (...
                l1 .* theta1 ...
                - l2 .* theta2 ...
                + 0.5 .* n  .* log(r1_sq./r2_sq));
            
            % Segment Axis Speed component
            s         =   gamma ./ (2*pi()) .* (theta2 - theta1);
            % Normal Axis speed component
            t         = - gamma ./ (4*pi()) .* log(r1_sq ./ r2_sq);      % Minus sign not present in Katz and Plotkin, but apparently needed!
            
            % And s (Segment Axis) speeds derivatives to theta angles
            s_dtheta1 =   gamma ./ (2*pi()) .* (       - 1     );
            s_dtheta2 =   gamma ./ (2*pi()) .* (1              );
            % and x,y
            s_dx      = s_dtheta1 * theta1_dx + s_dtheta2 * theta2_dx;
            s_dy      = s_dtheta1 * theta1_dy + s_dtheta2 * theta2_dy;
            
            % Now t (Normal direction) speed derivatives to (axis)projected distances
            t_dr1_sq  = - gamma ./ (4*pi()) ./     r1_sq;
            t_dr2_sq  =   gamma ./ (4*pi()) ./     r2_sq;
            % and chain rule to get x-y derivatives
            t_dx      =   t_dr1_sq .* r1_sq_dx  + t_dr2_sq  .* r2_sq_dx;
            t_dy      =   t_dr1_sq .* r1_sq_dy  + t_dr2_sq  .* r2_sq_dy;
            
            
            %%% Make consistency checks and corrections
            precision_threshold = 1e-15;
            % Condition for normal distance to be considered as tending to zero
            condition_on_segment = abs(n)<precision_threshold;
            s(condition_on_segment) = 0;
            if max(size(gamma)) == max(size(condition_on_segment))
                t(condition_on_segment) = - gamma(condition_on_segment) ./ (4*pi()) .* log(l1(condition_on_segment).^2 ./ l2(condition_on_segment).^2);
                
                % Fudge derivatives to 0 on segment (informally incorrect!)
                s_dx(condition_on_segment) = 0;
                s_dy(condition_on_segment) = 0;
                t_dx(condition_on_segment) = 0;
                t_dy(condition_on_segment) = 0;
            else
                t(condition_on_segment) = - gamma ./ (4*pi()) .* log(l1(condition_on_segment).^2 ./ l2(condition_on_segment).^2);
                
                % Fudge derivatives to 0 on segment (informally incorrect!)
                s_dx(condition_on_segment) = 0;
                s_dy(condition_on_segment) = 0;
                t_dx(condition_on_segment) = 0;
                t_dy(condition_on_segment) = 0;
            end
            
            %%% Transform (s,t) speed into (u,v) speed
            % Define l versor
            l_versor_x =  n_versor_y; %  (x2-x1)  ./ l;  % x component of tangential unit vector
            l_versor_y = -n_versor_x; %  (y2-y1)  ./ l;  % y component of tangential unit vector
            
            % Compute speed vector components (combine l,n system basis with s,t coefficients)
            u   = s .* l_versor_x + t .* n_versor_x;
            v   = s .* l_versor_y + t .* n_versor_y;
            % And do the same for the s-t speed derivatives, to get the u-v speeds
            % field jacobian, but in two steps!
            u_ds =      l_versor_x;
            u_dt =                        n_versor_x;
            v_ds =      l_versor_y                  ;
            v_dt =                        n_versor_y;
            % Now combine with chain rule
            u_dx = u_ds .* s_dx + u_dt .* t_dx;
            u_dy = u_ds .* s_dy + u_dt .* t_dy;
            v_dx = v_ds .* s_dx + v_dt .* t_dx;
            v_dy = v_ds .* s_dy + v_dt .* t_dy;
            % Done!!!
            
        end
        
        function [potential , u , v] = constant_strenght_vortex_segment_element_gpu(x, y, x1, y1, x2, y2, gamma)
            % Two Dimensional Constant Strenght Singularity Elements
            %
            % Inputs:
            %    x , y     % Target
            %    x1, y1    % Source Start
            %    x2, y2    % Source End
            %
            % Outputs:
            %    phi       % Potential
            %    u         % x speed component
            %    v         % y speed component
            %
            %
            % Reference:
            %    Based on:
            %            Katz and Plotkin, Chapter 10, Section 3 (p. 236)
            %    Represents analytical integral of:
            %             psi = - gamma / 2*pi * int(atan(y/(x-x0), x0, x1, x2);
            %    Generalized by embedded the coordinate transformation so that
            %    arbitrary y, y1 and y2 are allowed
            %
            
            x       = gpuArray(double(x    ));
            y       = gpuArray(double(y    ));
            x1      = gpuArray(double(x1   ));
            y1      = gpuArray(double(y1   ));
            x2      = gpuArray(double(x2   ));
            y2      = gpuArray(double(y2   ));
            gamma   = gpuArray(double(gamma));
            
            % Squared distances from target to source points
            r1_sq     = (x-x1).^2 + (y-y1).^2;
            r2_sq     = (x-x2).^2 + (y-y2).^2;
            
            % Projected Distances on Element Axis
            l_sq = (x2-x1).^2 + (y2-y1).^2;                 % Square of vortex segment lenght
            l    = sqrt(l_sq);                              % Vortex Segment lenght (always positive)
            
            l1   = ( (x-x1) .* (x2-x1) + (y-y1) .* (y2-y1) ) ./l;   % Distance to from target to start along element axis (can be positive or negative)
            l2   = ( (x-x2) .* (x2-x1) + (y-y2) .* (y2-y1) ) ./l;   % Distance to from target to end   along element axis (can be positive or negative)
            
            % Projected Distances on Normal Axis
            % Simple normal expression for 2d (easier than subtracting colinear
            % component, which is generic for n-D)
            n_versor_x = - (y2-y1)  ./ l;  % x component of normal unit vector
            n_versor_y =   (x2-x1)  ./ l;  % y component of normal unit vector
            
            % Normal distance
            n = (x-x1) .* n_versor_x + (y-y1) .* n_versor_y;
            % Altnernative normal distance
            % n = (x-x2) .* n_versor_x + (y-y2) .* n_versor_y
            
            % Angles from target to source points, as measured from the element axis
            % theta1    = atan(n ./ l1);
            % theta2    = atan(n ./ l2);
            
            theta1    = atan2(n , l1);
            theta2    = atan2(n , l2);
            
            % Potential (generalized for arbitraty y1 and y2, I imagine he thinks it is easy to find!)
            potential = - gamma / (2*pi()) * (...
                l1 .* theta1 ...
                - l2 .* theta2 ...
                + 0.5 .* n  .* log(r1_sq./r2_sq));
            
            % Segment Axis Speed component
            s         = gamma / (2*pi()) .* (theta2 - theta1);
            % Normal Axis speed componen
            t         = gamma / (4*pi()) .* log(r1_sq ./ r2_sq);
            
            %%% Transform (s,t) speed into (u,v) speed
            % Define l versor
            l_versor_x =  n_versor_y; %  (x2-x1)  ./ l;  % x component of tangential unit vector
            l_versor_y = -n_versor_x; %  (y2-y1)  ./ l;  % y component of tangential unit vector
            
            % Compute speed vector components (combine l,n system basis with s,t coefficients)
            u = s .* l_versor_x + t .* n_versor_x;
            v = s .* l_versor_y + t .* n_versor_y;
        end
        
        function [potential , u , v] = constant_strenght_vortex_segment_element_restricted(x, y, x1, x2, gamma)
            % Two Dimensional Constant Strenght Singularity Elements
            %
            % Inputs:
            %    x , y     % Target
            %    x1, y1    % Source Start
            %    x2, y2    % Source End
            %
            % Outputs:
            %    phi       % Potential
            %    u         % x speed component
            %    v         % y speed component
            %
            %
            % Reference:
            %    Based on:
            %            Katz and Plotkin, Chapter 10, Section 3 (p. 236)
            %    Represents analytical integral of:
            %             psi = - gamma / 2*pi * int(atan(y/(x-x0), x0, x1, x2);
            %    Original Expression : y1 and y2 restricted to 0; used for validation
            %    purposes!
            
            %%% (Katz and Plotkin expression, for y1=y2=0)
            % Angles from target to source points
            % (ATAN2 is not recomended here!)
            theta1    = atan2(y , (x-x1));
            theta2    = atan2(y , (x-x2));
            
            
            % Squared distances from target to source points
            r1_sq     = (x-x1).^2 + (y).^2;
            r2_sq     = (x-x2).^2 + (y).^2;
            
            
            % Potential
            potential = - gamma / (2*pi()) * (...
                (x-x1).* theta1 ...
                - (x-x2).* theta2 ...
                + 0.5 *  y    .* log(r1_sq./r2_sq));
            
            % x speed component
            u         = gamma / (2*pi()) .* (theta2 - theta1);
            % y speed component (singularity of r2sq can be dealt by inverting and
            % minus on log)
            v         =   gamma / (4*pi()) .* log(r1_sq ./ r2_sq);
             
        end
        
    end
    
end
