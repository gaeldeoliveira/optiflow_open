function [cf1_epRM , cf2_epRM, x_epRM, i_epRM, cf_scalar_epR] = select_compromise_in_pareto_front(cf1_pareto, cf2_pareto, x_pareto, epR, varargin)
    % Extract offset if applicable
    if ~isempty(varargin)
        i_offset = varargin{1};
    else
        i_offset = 0;
    end
    % Make scalar function
    cf_scalar_epR = (1 - epR) * cf1_pareto + epR* (cf2_pareto);
    % Find its maximum
    [~ , i_epRM] = max(cf_scalar_epR);
    % Make point pairs and coordinates of maximum
    cf1_epRM     = cf1_pareto(i_epRM + i_offset   ); 
    cf2_epRM     = cf2_pareto(i_epRM + i_offset   );
    x_epRM       = x_pareto(  i_epRM + i_offset, :);
end
