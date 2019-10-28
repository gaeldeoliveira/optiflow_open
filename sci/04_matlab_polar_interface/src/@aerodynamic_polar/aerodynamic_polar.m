classdef aerodynamic_polar
    %POLAR is a non-handle class storing a polar. In addition to storing
    %raw polar data it provides preprocessed interpolation routines
    
    properties
        % The data as it comes from R/Xfoil
        name = '';
        
        cl_max_local                % Cl of greatest local cl maximum
        alpha_cl_max_local          % Alpha of greatest local cl maximum
        
        alpha_stall_limit_range = 4 % Max range above end of monotonous region for detection of limited stall
        cl_max_local_limited        % Cl of greatest local cl maximum for alphas below  {end of monotonous region + alpha_stall_limit_range}
        alpha_cl_max_local_limited  % Alpha of greatest local cl maximum for alphas below  {end of monotonous region + alpha_stall_limit_range}
        
        invalid = false;
    end
    
    properties (SetAccess= private)
        pp_refs
        interpolation_method = 'linear'
        raw_data = struct();
        pp_alpha_cl
        cl_vector_cl_range = [];
    end
    
    methods
        function ap = aerodynamic_polar(raw_data, varargin)
            
            % The aerodynamic polar is a container for aerodynamic polar
            % data providing some handy interpolation and plotting
            % routines.
            %
            % The constructor function syntax is:
            % AP = AERODYNAMIC_POLAR(RAW_DATA,...)
            % AP = AERODYNAMIC_POLAR(RAW_DATA, NAME , ...)
            % AP = AERODYNAMIC_POLAR(RAW_DATA, NAME , INTERPOLATION_METHOD)
            %
            % RAW_DATA is a structure containing at least three array
            % fields:
            %   .alpha
            %   .cl
            %   .cd
            % Other fields will be interpolated in function of alpha and
            % are acessible through the field_alpha() method.
            %
            % NAME is a string specifying airfoil/polar name
            %
            % INTERPOLATION_METHOD can be any string supported in interp1:
            %  'linear' , 'pchip' , 'spline'
            % The default is 'linear', 'pchip' is also a good option, as it
            % is shape preserving.
            % Spline should be used with more caution as it is not shape
            % preserving. This property could be used to find a maximum
            % more easily
            
            % Assign name if supplied
            if ~isempty(varargin)
                ap.name = varargin{1};
                if length(varargin) == 2
                    ap.interpolation_method = varargin{2};
                end
            end
            
            % Assign Raw Data structure
            ap.raw_data = raw_data;
            % Make Piecewise Polynomials for interpolation
            ap = ap.make_pps();
            % Detect largest monotonous cl vs alpha region and make PP
            ap = ap.detect_monotonous_subregion();
        end
        
        
        function ap = make_pps(ap)
            % makes a piecewise interpolating polynomial for each field of
            % the raw data
            
            data_fields = fields(ap.raw_data);
            x = ap.raw_data.alpha;
            for n_field = 1:length(data_fields)
                % Extract data from structure
                y = ap.raw_data.(data_fields{n_field});
                
                % Filter Out NaNs
                y_non_nans = ~isnan(y);
                y_valid = y(y_non_nans);
                x_valid = x(y_non_nans);
                
                
                % Make interpolating polynomial
                pp = interp1(x_valid, y_valid, ap.interpolation_method , 'pp');
                % Store Interpolating polynomial
                ap.pp_refs.(data_fields{n_field}) = pp;
            end
            
        end
        
        function alpha = alpha_cl(ap , cl)
            if and(cl <= max(ap.cl_vector_cl_range) , cl >= min(ap.cl_vector_cl_range))
                alpha = ppval(ap.pp_alpha_cl , cl);
            else
                alpha = NaN;
            end
        end
        
        function cl = cl_alpha(ap , alpha)
            cl = ap.field_alpha('cl' , alpha);
        end
        
        function cd = cd_alpha(ap , alpha)
            cd = ap.field_alpha('cd' , alpha);
        end
        
        function cm = cm_alpha(ap , alpha)
            cm = ap.field_alpha('cm' , alpha);
        end
        
        function field = field_alpha(ap, field, alpha)
            if isempty(alpha)
                % Original:
                %   alpha = ap.raw_data.alpha(1):(1/5):ap.raw_data.alpha(end);
                % Changed to min():1/5:max() ? Probably more
                % robust for multi-streak polars (with BL and start point 
                % reinit for robustness)               
                alpha = min(ap.raw_data.alpha):(1/5):max(ap.raw_data.alpha);
            end
            if and(alpha <= max(ap.raw_data.alpha), alpha >= min(ap.raw_data.alpha))
                field = ppval(ap.pp_refs.(field)  ,  alpha);
            else
                field = Inf;
            end
        end
        
        function alpha_range = alpha_range(ap)
            alpha_range = min(ap.raw_data.alpha):(1/5):max(ap.raw_data.alpha);
        end
        
        function alpha_data = alpha_data(ap)
            alpha_data = ap.raw_data.alpha;
        end
        
        function [axes_vector] = plot(ap, varargin)
            
            % First Plot the Cl vs Cd
            
            if isempty(varargin)
                % Create subplot only if necessary
                axes_vector(1) = subplot('position',[0.1 0.1 0.40 0.8]);
                color ='b';
            else
                axes_vector = varargin{1};
                axes(axes_vector(1)); %#ok<MAXES>
                hold on
                color  = varargin{2};
            end
            if ischar(color)
                plot(ap.cd_alpha([]).*10^4 , ap.cl_alpha([]), color);
            else
                plot(ap.cd_alpha([]).*10^4 , ap.cl_alpha([]), 'Color', color);                
            end
            
            xlabel('C_d x 10^4')
            ylabel('C_l')
            grid on
            hold off
            
            a = axis;
            a(1) = 0;
            a(2) = 300;
            axis(a);
            
            
            % Second Plot the Cl vs Alpha
            
            if isempty(varargin)
                % Create subplot only if necessary
                axes_vector(2) = subplot('position',[0.60 0.1 0.35 0.8]);
            else
                axes(axes_vector(2)); %#ok<MAXES>
                hold on
            end
            
            if ischar(color)
                % Plot with color code
                plot(ap.raw_data.alpha(1):(1/5):ap.raw_data.alpha(end) , ap.cl_alpha([]), color , ap.alpha_cl_max_local_limited  , ap.cl_max_local_limited , 'rx');
            else
                % Plot with color RGB vector
                plot(ap.raw_data.alpha(1):(1/5):ap.raw_data.alpha(end) , ap.cl_alpha([]), 'Color' , color )
                hold on
                plot(ap.alpha_cl_max_local_limited  , ap.cl_max_local_limited , 'rx');
            end
            
            cl_cd_max = max(ap.raw_data.cl ./ ap.raw_data.cd);
            
            xlabel('Angle of Attack (Degrees)');
            ylabel('C_l');
            legend(['Cl - max(C_l/C_d) = ' num2str(cl_cd_max)]  , 'Dectected Stall Point', 'Location', 'SouthEast');
            grid on
            hold off
            
            
            
        end
        
        function ap = detect_monotonous_subregion(ap)
            
            % Initialize variables to start for loop
            n_data = 2;
            % Determine direction at beginning of polar
            direction = sign(ap.raw_data.alpha(n_data) - ap.raw_data.alpha(n_data-1));
            
            if direction > 0
                % If we are growing then change to
                %       growing occurs at start
                vector_n_change_to_growing = 1;
                %       decreasing is yet to come
                vector_n_change_to_decreasing = [];
                
                % And Cl
                vector_cl_change_to_growing = ap.raw_data.cl(1);
                vector_cl_change_to_decreasing = [];
            else
                % If we are decreasing then change to:
                %       growing is yet to come
                vector_n_change_to_growing = [];
                %       decreasing occurs at start (from nothing to decreasing)
                vector_n_change_to_decreasing = 1;
                
                % And Cl
                vector_cl_change_to_growing = [];
                vector_cl_change_to_decreasing = ap.raw_data.cl(1);
            end
            
            
            % Run loop to determine vector of direction changes and their
            % corresponding cl ranges
            
            for n_data = 3:length(ap.raw_data.alpha)
                % Store previous direction
                direction_nminus1 = direction;
                % Determine New Direction
                direction = sign(ap.raw_data.cl(n_data) - ap.raw_data.cl(n_data-1));
                
                % If there is a direction change store by appending to changes vector
                % Notice that we store the point n-1 as this is the one
                % that corresponds to the extreme (point of change)
                if direction > direction_nminus1
                    vector_n_change_to_growing = [vector_n_change_to_growing (n_data-1)]; %#ok<AGROW>
                    vector_cl_change_to_growing = [vector_cl_change_to_growing ap.raw_data.cl(n_data-1)]; %#ok<AGROW>
                    %        disp('Change to Growing')
                    
                    % Handle direction = 0 case
                    direction = 1;
                end
                
                if direction < direction_nminus1
                    vector_n_change_to_decreasing = [vector_n_change_to_decreasing (n_data-1)]; %#ok<AGROW>
                    vector_cl_change_to_decreasing = [vector_cl_change_to_decreasing ap.raw_data.cl(n_data-1)]; %#ok<AGROW>
                    %       disp('Change to Decreasing')
                    
                    % Handle direction = 0 case
                    direction = -1;
                end
                
            end
            
            % Pad Vectors
            % In some cases, the length of the two vectors will be different. When the
            % total the number of changes is even, the shortest vector should be filled
            % with the last point of the polar
            if length(vector_cl_change_to_growing) > length(vector_cl_change_to_decreasing)
                vector_n_change_to_decreasing(end+1)        = length(ap.raw_data.cl);
                vector_cl_change_to_decreasing(end+1)       = ap.raw_data.cl(end);
            end
            
            if length(vector_cl_change_to_growing) < length(vector_cl_change_to_decreasing)
                vector_n_change_to_growing(end+1)        = length(ap.raw_data.cl);
                vector_cl_change_to_growing(end+1)       = ap.raw_data.cl(end);
            end
            
            % Determine size of monotonous cl intervals
            vector_cl_ranges = vector_cl_change_to_decreasing - vector_cl_change_to_growing;
            
            % Find greatest Cl range
            [~, index_cl_range_largest] = max(vector_cl_ranges);
            
            % Extract Indices of boundaries of cl range
            n_cl_range_start = vector_n_change_to_growing(index_cl_range_largest);
            n_cl_range_end   = vector_n_change_to_decreasing(index_cl_range_largest);
            
            % Extract and store Cl and Alpha vectors for interpolation
            ap.cl_vector_cl_range   =  ap.raw_data.cl(n_cl_range_start:n_cl_range_end);
            alpha_vector_cl_range   =  ap.raw_data.alpha(n_cl_range_start:n_cl_range_end);
            
            % Make and store interpolation Polynomial
            x_valid = ap.cl_vector_cl_range;
            y_valid = alpha_vector_cl_range;
            
            ap.pp_alpha_cl = interp1(x_valid, y_valid, ap.interpolation_method , 'pp');
            
            % plot(x_valid, y_valid)
            % xlabel('Cl')
            % ylabel('Cd')
            
            
            % Another Strategy for determining cl stall                                    
            
            % Find AOA's of local maxima
            vector_alpha_change_to_decreasing = ...
                transpose(ap.raw_data.alpha(vector_n_change_to_decreasing));
            
            % Find AOA of end of monotonous interval
            alpha_end_monotonous = max(alpha_vector_cl_range);
            
            % Look for stall angle in the alpha_stall_limit_range degrees that follow
            % the end of the longest monotonous interval
            % Look for angles of cl change to decreasing close enough to the end of the
            % monotonous interval to be candidate stall angles
            
            index_possible_stall = vector_alpha_change_to_decreasing < ...
                (alpha_end_monotonous+ap.alpha_stall_limit_range);
            
            vector_alpha_possible_stall = ...
                vector_alpha_change_to_decreasing(index_possible_stall);
            vector_cl_possible_stall = ...
                vector_cl_change_to_decreasing(index_possible_stall);
            
            % Now find the highest local(which can also be absolute) maximum in the
            % selected candidates
            
            [ap.cl_max_local_limited index_cl_max_local] = ...
                max(vector_cl_possible_stall);
            
            ap.alpha_cl_max_local_limited = ...
                    vector_alpha_possible_stall(index_cl_max_local);
            
%             % % % This if clause might not be necessary (probably doing nothing thanks
%             % % % to padding... but leave for now, as unexpected failures of monotonous
%             % % % detection could occur)
%             if isscalar(ap.cl_max_local_limited)
%                 % A stall angle was found
% 
%             else
%                 % No local maximum was found, set maximum angle of attack as
%                 % angle of maximum cl overall ( I think this case can't happen, unless
%                 % monotonous interval detection failed...)
%                 [ap.cl_max_local_limited , n_cl_max_overall] = max(ap.raw_data.cl);
%                 ap.alpha_cl_max_local_limited = ap.raw_data.alpha(n_cl_max_overall);
%             end                                    
        end
        
    end
end


