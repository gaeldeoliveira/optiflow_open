classdef case_comparator
    %CASE_COMPARATOR is a static class meant to store static methods for
    %comparing airfoil polars
    % (IMPORTANT CORRECTION: sqrt were not taken in the right place in
    % editions before and up to TORQUE2018 abstract. Corrected 21/02/2018 )
    
    properties
        
    end
    
    methods
        function comparator = case_comparator()
            % Dummy Constructor!
        end
    end
    
    methods(Static)
        
        function accuracy_metric = compare_cl_data(sim_results , exp_case)
            if strcmp(class(sim_results.ap), 'aerodynamic_polar') %#ok<STISA>
                N_quad_points = 100; % TODO: add handling of polars without points on one of polars (shoots out NaN at the moment)
                % Now find ranges of alpha for which values are defined in both cases
                % (simulated and experimental)
                alpha_range_sim = sim_results.ap.alpha_range();
                alpha_range_exp = exp_case.al_alcl_polar_std;
                
                % Now find overlap between simulated and experimental ranges
                alpha_range_overlap_min = max(min(alpha_range_sim), min(alpha_range_exp)) + eps();
                alpha_range_overlap_max = min(max(alpha_range_sim), max(alpha_range_exp)) - eps();
                
                % Make range of overlapping angles of attack for reinterpolation and integration
                alpha_range_overlap = linspace(alpha_range_overlap_min, alpha_range_overlap_max, N_quad_points);
                
                % Interpolate simulated and experimental cl values
                cl_sim = sim_results.ap.cl_alpha(alpha_range_overlap);
                cl_exp = exp_case.cl_alpha(alpha_range_overlap);
                
                % Compute difference in cl, between simulated and experimental case
                delta_cl = cl_sim - cl_exp;
                
                % Now integrate squared differences with trapezoidal rule (Z = trapz(X,Y))
                RMS_delta_cl = sqrt(trapz(alpha_range_overlap , delta_cl.^2));
                
                % Finally, normalize with length of integration region
                % (IMPORTANT CORRECTION: sqrt added after TORQUE2018 abstract, this corresponds)
                % normalized_RMS_delta_cl = RMS_delta_cl / (alpha_range_overlap_max - alpha_range_overlap_min);
                normalized_RMS_delta_cl = RMS_delta_cl / sqrt(alpha_range_overlap_max - alpha_range_overlap_min);
                
                % And keep value
                accuracy_metric = normalized_RMS_delta_cl;
                
                % Handle missing points as a linear penalty, in case simulation data (and hence overlap) does not span all experimental data range.
                if or(alpha_range_overlap_min > min(alpha_range_exp), alpha_range_overlap_max < max(alpha_range_exp))
                    accuracy_metric = accuracy_metric * (alpha_range_overlap_max - alpha_range_overlap_min) / (max(alpha_range_exp) - min(alpha_range_exp));
                end
            else
                % If aerodynamic polar was invalid, set accuracy index to
                % NaN
                accuracy_metric = nan;
            end
        end
        
        function accuracy_metric = compare_cm_data(sim_results , exp_case)
            if strcmp(class(sim_results.ap), 'aerodynamic_polar') %#ok<STISA>
                N_quad_points = 100; % TODO: add handling of polars without points on one of polars (shoots out NaN at the moment)
                % Now find ranges of alpha for which values are defined in both cases
                % (simulated and experimental)
                alpha_range_sim = sim_results.ap.alpha_range();
                alpha_range_exp = exp_case.al_alcm_polar_std;
                
                % Now find overlap between simulated and experimental ranges
                alpha_range_overlap_min = max(min(alpha_range_sim), min(alpha_range_exp)) + eps();
                alpha_range_overlap_max = min(max(alpha_range_sim), max(alpha_range_exp)) - eps();
                
                % Make range of overlapping angles of attack for reinterpolation and integration
                alpha_range_overlap = linspace(alpha_range_overlap_min, alpha_range_overlap_max, N_quad_points);
                
                % Interpolate simulated and experimental cl values
                cm_sim = sim_results.ap.cm_alpha(alpha_range_overlap);
                cm_exp = exp_case.cm_alpha(alpha_range_overlap);
                
                % Compute difference in cl, between simulated and experimental case
                delta_cm = cm_sim - cm_exp;
                
                % Now integrate squared differences with trapezoidal rule (Z = trapz(X,Y))
                RMS_delta_cm = sqrt(trapz(alpha_range_overlap , delta_cm.^2));
                
                % Finally, normalize with length of integration region
                % (IMPORTANT CORRECTION: sqrt added after TORQUE2018 abstract, this corresponds)
                % normalized_RMS_delta_cm = RMS_delta_cm / (alpha_range_overlap_max - alpha_range_overlap_min);
                normalized_RMS_delta_cm = RMS_delta_cm / sqrt(alpha_range_overlap_max - alpha_range_overlap_min);
                
                % And keep value
                accuracy_metric = normalized_RMS_delta_cm;
                
                % Handle missing points as a linear penalty, in case simulation data (and hence overlap) does not span all experimental data range.
                if or(alpha_range_overlap_min > min(alpha_range_exp), alpha_range_overlap_max < max(alpha_range_exp))
                    accuracy_metric = accuracy_metric * (alpha_range_overlap_max - alpha_range_overlap_min) / (max(alpha_range_exp) - min(alpha_range_exp));
                end
            else
                % If aerodynamic polar was invalid, set accuracy index to
                % NaN
                accuracy_metric = nan;
            end
        end
        
        function accuracy_metric = compare_cd_data(sim_results , exp_case)
            if strcmp(class(sim_results.ap), 'aerodynamic_polar') %#ok<STISA>
                N_quad_points = 100; % TODO: add handling of polars without points on one of polars (shoots out NaN at the moment)
                % Now find ranges of alpha for which values are defined in both cases
                % (simulated and experimental)
                alpha_range_sim = sim_results.ap.alpha_range();
                alpha_range_exp = exp_case.al_alclcd_polar_std;
                
                % Now find overlap between simulated and experimental ranges
                alpha_range_overlap_min = max(min(alpha_range_sim), min(alpha_range_exp)) + eps();
                alpha_range_overlap_max = min(max(alpha_range_sim), max(alpha_range_exp)) - eps();
                
                % Make range of overlapping angles of attack for reinterpolation and integration
                alpha_range_overlap = linspace(alpha_range_overlap_min, alpha_range_overlap_max, N_quad_points);
                
                % Interpolate simulated and experimental cl values
                cd_sim = sim_results.ap.cd_alpha(alpha_range_overlap);
                cd_exp = exp_case.cd_alpha(alpha_range_overlap);
                
                % Compute difference in cl, between simulated and experimental case
                delta_cd = cd_sim - cd_exp;
                
                % Now integrate squared differences with trapezoidal rule (Z = trapz(X,Y))
                RMS_delta_cd = sqrt(trapz(alpha_range_overlap , delta_cd.^2));
                
                % Finally, normalize with length of integration region
                % (IMPORTANT CORRECTION: sqrt added after TORQUE2018 abstract, this corresponds)
                % normalized_RMS_delta_cd = RMS_delta_cd / (alpha_range_overlap_max - alpha_range_overlap_min);
                normalized_RMS_delta_cd = RMS_delta_cd / sqrt(alpha_range_overlap_max - alpha_range_overlap_min);
                
                % And keep value
                accuracy_metric = normalized_RMS_delta_cd;
                
                % Handle missing points as a linear penalty, in case simulation data (and hence overlap) does not span all experimental data range.
                if or(alpha_range_overlap_min > min(alpha_range_exp), alpha_range_overlap_max < max(alpha_range_exp))
                    accuracy_metric = accuracy_metric * (alpha_range_overlap_max - alpha_range_overlap_min) / (max(alpha_range_exp) - min(alpha_range_exp));
                end
            else
                % If aerodynamic polar was invalid, set accuracy index to
                % NaN
                accuracy_metric = nan;
            end
        end
        
        function global_accuracy_metric = group_case(accuracy_index_array)
            % Find index of valid points (not NaN, which means invalid polar)
            not_nan = not(isnan(accuracy_index_array));
            % Count number of valid points
            n_valid_points = sum(not_nan);
            % Count total number of points
            n_total_points = length(accuracy_index_array);
            % RMS sum and normalize
            % (IMPORTANT CORRECTION: sqrt moved after TORQUE2018 abstract, new version is correct RMS logic)
            % global_accuracy_metric = sqrt(sum(accuracy_index_array(not_nan).^2)) / n_valid_points;
            global_accuracy_metric = sqrt(sum(accuracy_index_array(not_nan).^2) / n_valid_points);
            % Penalize (no effect if all points are valid)
            global_accuracy_metric = global_accuracy_metric * (n_total_points / n_valid_points);
        end
        
        function normalized_accuracy_profile = normalize_accuracy_profile(accuracy_profile, accuracy_profile_0)
            % Normalize accuracy_profile with respect to reference accuracy_profile_0
            
            % Get field names of accuracy profile
            accuracy_profile_fields = fieldnames(accuracy_profile);
            
            % Allocate normalized accuracy profile
            normalized_accuracy_profile = simulation_accuracy_profile();
            
            for n_accuracy_profile_field = 1:length(accuracy_profile_fields)
                normalized_accuracy_profile.(accuracy_profile_fields{n_accuracy_profile_field}) = ...
                    accuracy_profile.(accuracy_profile_fields{n_accuracy_profile_field}) ./ accuracy_profile_0.(accuracy_profile_fields{n_accuracy_profile_field});
            end
        end
        
        function normalized_accuracy_profile = normalize_accuracy_profile_on_minibatch(accuracy_profile, accuracy_profile_0, index_minibatch)
            % Normalize accuracy_profile with respect to reference accuracy_profile_0
            % Later on do this as a method of each type of
            % simulation_accuracy_profile() with 
            
            % Allocate normalized accuracy profile
            normalized_accuracy_profile = simulation_accuracy_profile();
            
            normalized_accuracy_profile.cl_global_accuracy  = accuracy_profile.cl_global_accuracy ./ accuracy_profile_0.cl_global_accuracy; % Accuracy of Cl computation on airfoil_polar cases
            normalized_accuracy_profile.cm_global_accuracy  = accuracy_profile.cm_global_accuracy ./ accuracy_profile_0.cm_global_accuracy; % Accuracy of Cm computation on airfoil_polar cases
            normalized_accuracy_profile.cd_global_accuracy  = accuracy_profile.cd_global_accuracy ./ accuracy_profile_0.cd_global_accuracy; % Accuracy of Cd computation on airfoil_polar cases
        
            normalized_accuracy_profile.cl_accuracy_metric_array =  accuracy_profile.cl_accuracy_metric_array ./ accuracy_profile_0.cl_accuracy_metric_array(index_minibatch);  % Accuracy of Cl computation for each airfoil_polar case (index_minibatch reordering)
            normalized_accuracy_profile.cm_accuracy_metric_array =  accuracy_profile.cm_accuracy_metric_array ./ accuracy_profile_0.cm_accuracy_metric_array(index_minibatch);  % Accuracy of Cm computation for each airfoil_polar case (index_minibatch reordering)
            normalized_accuracy_profile.cd_accuracy_metric_array =  accuracy_profile.cd_accuracy_metric_array ./ accuracy_profile_0.cd_accuracy_metric_array(index_minibatch);  % Accuracy of Cd computation for each airfoil_polar case (index_minibatch reordering)
        end
        
        function scalar_accuracy_metric = make_scalar_accuracy_metric_from_normalized_accuracy_profile(normalized_accuracy_profile)
            % Simple RMS sum for now
            scalar_accuracy_metric = sqrt(normalized_accuracy_profile.cl_global_accuracy.^2 + ...
                                          normalized_accuracy_profile.cm_global_accuracy.^2 + ...
                                          normalized_accuracy_profile.cd_global_accuracy.^2);
        end
    end
end

