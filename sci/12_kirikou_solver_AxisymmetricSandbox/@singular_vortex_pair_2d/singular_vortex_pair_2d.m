classdef singular_vortex_pair_2d < handle
    %singular_vortex_pair_2d
    %   Holds information and helpful functions to compute effect a pair of
    %   lifting vortices (symmetry on x axis)
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
        x_v         % x-stance of lifting vortex pair (symmetry on x axis)
        r_v         % distance of lifting vortex pair to x axis (symmetry on x axis)
        gamma_v     % strenght of each vortex in the pair
    end
    
    methods
        function VP = singular_vortex_pair_2d(x_v, r_v, gamma_v)
           % Constructor function 
           VP.x_v = x_v;
           VP.r_v = r_v;
           VP.gamma_v = gamma_v;
        end
        
        function [pot, u , v] = induced_speed_on_many_points(VP, x, y)
            [u , v] = induction_functions.constant_strenght_singular_vortex_pair(x, y, VP.x_v, VP.r_v, VP.gamma_v);
            % Make a dummy potential (temporary)(zero array, with right size!)
            pot = u * 0 + v * 0;
        end
    end
    
end

