classdef disk_loading < handle
    % Disk_Loading is a class that describes the loading distribution of an
    %   actuator disk
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
        % % New Model
        % % Define Average Loading
        f1_average             % Average Loading Density
        
        % % Define Loading Shape Coefficient (default: 1 0 1 1 1 1)
        ab0  =   1.0;          % Average Value (Relative, typically 1)
        ab2  =   0.0;          % Tip Offset from Average Value (Relative)
        ab4  =   1.0;          % 4th order Curtosis Parameter (1.0 = order four  and above)
        ab6  =   1.0;          % 6th order Curtosis Parameter (1.0 = order six   and above)
        ab8  =   1.0;          % 8th order Curtosis Parameter (1.0 = order eight and above)
        ab10 =   1.0;          % 8th order Curtosis Parameter (1.0 = order eight and above)
        
        % % Define Generator (B) Functions (for all n but 0)
        bn = @(eta , n) ( eta.^n - 1/(n+1) ) / (1 - 1/(n+1));
        
        % Build a flag to indicate whether we are in the disk
        f_w_disk_flag = @(eta) and(eta > -1, eta < 1);
        
        
    end
    
    methods
        function DL = disk_loading(f1_average)
            % Store inputs
            DL.f1_average = f1_average;
        end
        
        function f_w = f_w_function(DL, y_over_r)
            % Returns loading density at some y_over_r position

            % Particular Generator (B) Functions 
            b0  = @(eta) eta.^0;
            b2  = @(eta) DL.bn(eta, 2 );
            b4  = @(eta) DL.bn(eta, 4 );
            b6  = @(eta) DL.bn(eta, 6 );
            b8  = @(eta) DL.bn(eta, 8 );
            b10 = @(eta) DL.bn(eta, 10);
        
            % % Define Adimensional Loading Function (recursive linear combination)
            C10 = @(eta)            b10(eta);
            C8  = @(eta) (1-DL.ab10) * b8( eta) + DL.ab10 * C10(eta);
            C6  = @(eta) (1-DL.ab8 ) * b6( eta) + DL.ab8  * C8( eta);
            C4  = @(eta) (1-DL.ab6 ) * b4( eta) + DL.ab6  * C6( eta);
            C2  = @(eta) (1-DL.ab4 ) * b2( eta) + DL.ab4  * C4( eta);
            C0  = @(eta)    DL.ab0   * b0( eta) + DL.ab2  * C2( eta);

            % Compute loading and fudge to zero outside actuator
            f_w = C0(y_over_r) .* DL.f_w_disk_flag(y_over_r) * DL.f1_average; 
        end
    end
    
end

