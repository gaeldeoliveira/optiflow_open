classdef experimental_polar_conditions
    %EXPERIMENTAL_POLAR_CONDITIONS is a simple class for holding data. For
    %now, we just use it as structure template (non-handle, without methods)
    
    properties
        M                           % Mach number
        Re                          % Reynolds Number
        N_crit                      % Critical Amplification Factor of Tollmien-Schlichted waves (Van Ingen e^N method)
        xtr_top                     % x/c of forced transition point for top side (1 = free)
        xtr_bot                     % x/c of forced transition point for bottom side (1 = free)
        INGE = false                % boolean for toggling (when true) non-linear growth of Tollmien Schlichting waves (van Ingen method)
    end
    
    methods
        function EPC = experimental_polar_conditions()
            % Dummy constructor! This class exists mostly for documentation
            % purposes!
        end
    end
end

