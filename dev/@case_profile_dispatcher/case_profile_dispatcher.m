classdef case_profile_dispatcher < handle
    %CASE_PROFILE_DISPATCHER receives profile evaluation requests and dispatches
    % them to the appropriate evaluators. It returns accuracy metrics.
    %
    % Evaluators are wrapped in functions and set their own local context
    % to enable distributed computing in the future (with batch/job/task
    % scheduling). 
    
    properties
        EDB     % Experimental case database (handle)
        CPR     % Cluster profile (value class)
    end
    
    methods
        function CPD = case_profile_dispatcher(EDB, CPR)
            %CASE_PROFILE_DISPATCHER Constructor receives and experimental
            % case database (EPR) and a cluster profile (CPR). 
            
            % Store handle to EDB object
            CPD.EDB = EDB;
            % Store handle to CPR object
            CPD.CPR = CPR;
        end
        
        function [accuracy_profile, APE] = dispatch_simulation_profile(CPD, SPR)
            % Receives a simulation profile (SPR) and returns global
            % accuracy metrics as an accuracy profile.
            % Optional outputs include (APE and BLE?)objects for more
            % detailed inspection of the results of executed simulations.
            
            % Find airfoil polar cases
            index_airfoil_cases = CPD.EDB.find_cases_by_IDkind('airfoil_polar');
            % And get them into a cell array
            EC_cell_airfoil_cases = CPD.EDB.EC_cell(index_airfoil_cases);
            
            % Get database parametrization orders
            u_order = CPD.EDB.SD.parametrization_handles.upper.order;
            l_order = CPD.EDB.SD.parametrization_handles.lower.order;
            
            % Dispatch airfoi polar cases
            % [cl_global_accuracy , cm_global_accuracy , cd_global_accuracy, APE] = APE_function_wrapper(CPD.CPR, SPR, EC_cell_airfoil_cases);
            [cl_global_accuracy , cm_global_accuracy , cd_global_accuracy, APE] = APE_function_wrapper_custom_p(CPD.CPR, SPR, EC_cell_airfoil_cases, u_order, l_order);
            
            % Replicate above lines for other types of simulation (e.g.
            % bl_run)
            
            % Create a simulation_accuracy_profile
            accuracy_profile = simulation_accuracy_profile();
            % Feed results into it!
            accuracy_profile.cl_global_accuracy       = cl_global_accuracy;
            accuracy_profile.cm_global_accuracy       = cm_global_accuracy;
            accuracy_profile.cd_global_accuracy       = cd_global_accuracy;
            % And also by airfoil!
            accuracy_profile.cl_accuracy_metric_array = APE.cl_accuracy_metric_array;
            accuracy_profile.cm_accuracy_metric_array = APE.cm_accuracy_metric_array;
            accuracy_profile.cd_accuracy_metric_array = APE.cd_accuracy_metric_array;
        end
        
        function [accuracy_profile, APE] = dispatch_simulation_profile_on_minibatch(CPD, SPR, index_minibatch)
            % Receives a simulation profile (SPR) and returns global
            % accuracy metrics as an accuracy profile.
            % Optional outputs include (APE and BLE?)objects for more
            % detailed inspection of the results of executed simulations.
            
            % Find airfoil polar cases that belong to minibatch
            % index_airfoil_cases = CPD.EDB.find_cases_by_IDkind('airfoil_polar');
            index_airfoil_cases = CPD.EDB.find_cases_by_IDkind_on_dbsubset('airfoil_polar', index_minibatch);
            
            % And get them into a cell array
            EC_cell_airfoil_cases = CPD.EDB.EC_cell(index_airfoil_cases);
            
            % Get database parametrization orders
            u_order = CPD.EDB.SD.parametrization_handles.upper.order;
            l_order = CPD.EDB.SD.parametrization_handles.lower.order;
            
            % Dispatch airfoi polar cases
            %[cl_global_accuracy , cm_global_accuracy , cd_global_accuracy, APE] = APE_function_wrapper(         CPD.CPR, SPR, EC_cell_airfoil_cases);
            [cl_global_accuracy , cm_global_accuracy , cd_global_accuracy, APE] = APE_function_wrapper_custom_p(CPD.CPR, SPR, EC_cell_airfoil_cases, u_order, l_order);
            
            % Replicate above lines for other types of simulation (e.g.
            % bl_run)
            
            % Create a simulation_accuracy_profile
            accuracy_profile = simulation_accuracy_profile();
            % Feed results into it!
            accuracy_profile.cl_global_accuracy       = cl_global_accuracy;
            accuracy_profile.cm_global_accuracy       = cm_global_accuracy;
            accuracy_profile.cd_global_accuracy       = cd_global_accuracy;
            % And also by airfoil!
            accuracy_profile.cl_accuracy_metric_array = APE.cl_accuracy_metric_array;
            accuracy_profile.cm_accuracy_metric_array = APE.cm_accuracy_metric_array;
            accuracy_profile.cd_accuracy_metric_array = APE.cd_accuracy_metric_array;
        end
        
    end
end

