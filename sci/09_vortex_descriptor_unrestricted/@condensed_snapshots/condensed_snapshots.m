classdef condensed_snapshots < handle
    %CONDENSED_SNAPSHOTS collects data from interaction field snapshots
    %into a large 3d array (z,y,n_snapshot)
    
    properties
        u_tilde_over_ue  = [];  % Stored snapshots 
        y_v_over_d0_data = [];  % Stored vortex core heights
        z_v_over_S_data  = [];  % Stored vortex core spanwise positions
        
        y_over_d_mesh    = [];  % Spatial mesh for reconstructing images
        z_over_s_mesh    = [];  % Spatial mesh for reconstructing images
    end
    
    methods
        function CSN = condensed_snapshots()
            % Dummy contructor method
        end
        
        function add_snapshots(CSN, snapshots)
            % Adds snapshots from a new source to main arrays after
            % checking their validity
            for n_snapshot = 2:size(snapshots.u_tilde_data, 3)
                % Test snapshot validity
                snapshot_too_large = max(max(abs(snapshots.u_tilde_data(:,:,n_snapshot)))) / snapshots.ue_data(n_snapshot) > 1;
                snapshot_has_inf   = not(sum(find(isinf(snapshots.u_tilde_data(:,:,n_snapshot)))) == 0);
                snapshot_has_nan   = not(sum(find(isnan(snapshots.u_tilde_data(:,:,n_snapshot)))) == 0);
                snapshot_undefined = or(sum(sum(snapshots.u_tilde_data(:,:,n_snapshot))) == 0, snapshots.ue_data(n_snapshot) == 0);
                % If it is valid (not invalid! inclusive or)
                if not(or(snapshot_too_large , or(snapshot_has_inf, or(snapshot_has_nan, snapshot_undefined))))
                    % Add it
                    disp(['Adding snapshot Nr: ' , num2str(n_snapshot) , '    from case:' , snapshots.vg_case.friendly_name])
                    if isempty(CSN.u_tilde_over_ue)
                        CSN.u_tilde_over_ue(:,:,    1) = snapshots.u_tilde_data(:,:,n_snapshot) / snapshots.ue_data(n_snapshot);
                        CSN.y_v_over_d0_data(       1) = snapshots.y_v_over_d0_data(n_snapshot) ;
                        CSN.z_v_over_S_data(        1) = snapshots.z_v_over_S_data( n_snapshot) ;
                    else
                        CSN.u_tilde_over_ue(:,:,end+1) = snapshots.u_tilde_data(:,:,n_snapshot) / snapshots.ue_data(n_snapshot);
                        CSN.y_v_over_d0_data(   end+1) = snapshots.y_v_over_d0_data(n_snapshot) ;
                        CSN.z_v_over_S_data(    end+1) = snapshots.z_v_over_S_data( n_snapshot) ;
                    end
                else
                    disp(['Reject snapshot Nr: ' , num2str(n_snapshot) , '    from case:' , snapshots.vg_case.friendly_name])
                end
            end
            
            % Get mesh from any structure
            CSN.y_over_d_mesh = snapshots.y_mesh / snapshots.delta0;
            CSN.z_over_s_mesh = snapshots.z_mesh / snapshots.S     ;
        end
    end
end

