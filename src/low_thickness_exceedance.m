function [pn75, pn90, bh_75, bh_90]= low_thickness_exceedance(parameters, experiment_result)
%LOW_THICKNESS_PENALIZATION

    % Extract thresholds if applicable
    if isfield(parameters, 'bh_75_threshold')
        % Typical thresholds for 21% thick airfoils 
        bh_75_threshold = parameters.bh_75_threshold;
    else
        % Typical thresholds for 21% thick airfoils 
        % bh_75_threshold = 0.0591; % 0.75 * bh_75 of   S809
        bh_75_threshold = 0.0320;   % 0.6  * bh_75 of Riso21
    end
    if isfield(parameters, 'bh_90_threshold')
        % Typical thresholds for 21% thick airfoils 
        bh_90_threshold = parameters.bh_90_threshold;
    else
        % Typical thresholds for 21% thick airfoils 
        % bh_90_threshold = 0.0160; % 0.8 * bh_90 of  S809
        bh_90_threshold = 0.0117; %   0.8 * bh_90 of FFA21
    end

    % Extract coordinates
    tx = experiment_result.coordinates.tx_analysis;
    tz = experiment_result.coordinates.tz_analysis;

    % Find leading edge
    [~   , i_le] = min(tx);
    % Top
    tx_top       = flipud(tx(1:i_le));
    tz_top       = flipud(tz(1:i_le));
    % Bot
    tx_bot       = tx(i_le:end);
    tz_bot       = tz(i_le:end);
    % Interpolate at 75% chord
    tz_top_75    = interp1(tx_top, tz_top, 0.75);
    tz_bot_75    = interp1(tx_bot, tz_bot, 0.75);
    % Interpolate at 90% chord
    tz_top_90    = interp1(tx_top, tz_top, 0.90);
    tz_bot_90    = interp1(tx_bot, tz_bot, 0.90);
    % Find building heigh at both chords
    bh_75        = tz_top_75 - tz_bot_75;
    bh_90        = tz_top_90 - tz_bot_90;
    % Compute threshold exceedances
    bh_75_exceedance = (bh_75 - bh_75_threshold) / bh_75_threshold;
    bh_90_exceedance = (bh_90 - bh_90_threshold) / bh_90_threshold;
    
    % Compute penalty for 75%
    if bh_75_exceedance > 0
        pn75 = 0;
    else
        pn75 = abs(bh_75_exceedance);
        disp(['Pn75 =' , num2str(pn75)]);
    end
    % Compute penalty for 90%
    if bh_90_exceedance > 0
        pn90 = 0;
    else
        pn90 = abs(bh_90_exceedance);
        disp(['Pn90 =' , num2str(pn90)]);
    end

end
    


