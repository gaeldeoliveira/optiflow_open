function penalty_factor = penalty_factor_building_height_violation(tx_no_check, min_thickness, coordinates)
%Computes a penalty factor for violation of a minimum building height
%constraint. This is a multiplicative factor:
%   * 1 - means full constraint violation
%   * 0 - means full constraint violation
% Value is fudged for security, ensuring return is always between 1 and 0.

% Make criteria for side crossing or very thin
%tx_no_check    = 0.025;
%min_thickness  = 0.003;

% Extract Coordinates
tx = coordinates.tx;
tz = coordinates.tz;

% Find Leading Edge
[~ , i_LE] = min(tx);
% And make sorting indexes
i_upper = i_LE:-1:1;
i_lower = i_LE:1:length(tx);
% Now sort upper and lower sides, first for x-component
tx_upper = tx(i_upper);
tx_lower = tx(i_lower);
% Now sort upper and lower sides, then for z-component
tz_upper = tz(i_upper);
tz_lower = tz(i_lower);
% Verify
% plot(tx, tz); hold on; axis equal; grid on;
% plot(tx_upper, tz_upper, 'o');
% plot(tx_lower, tz_lower, 'x');

% Now filter out points before and after limit, first for x-component
tx_upper_for_check = tx_upper(and(tx_upper > tx_no_check, tx_upper < (1-tx_no_check))); %#ok<NASGU>
tx_lower_for_check = tx_lower(and(tx_lower > tx_no_check, tx_lower < (1-tx_no_check))); %#ok<NASGU>
% Now filter out points before and after limit, then for z-component
tz_upper_for_check = tz_upper(and(tx_upper > tx_no_check, tx_upper < (1-tx_no_check)));
tz_lower_for_check = tz_lower(and(tx_lower > tx_no_check, tx_lower < (1-tx_no_check)));

% New verification
% plot(tx, tz); hold on; axis equal; grid on;
% plot(tx_upper_for_check, tz_upper_for_check, 'o');
% plot(tx_lower_for_check, tz_lower_for_check, 'x');

% Compute local building height!
tz_height = tz_upper_for_check - tz_lower_for_check;

% And determine penalty factor
relative_violation = (tz_height - min_thickness) ./ min_thickness;
penalty_factor = 1 - mean(abs(relative_violation(relative_violation < 0)));

% Fudge into 0..1 range, should be there, but let's just ensure everything
% is going fine:
penalty_factor = max([0 , sum(penalty_factor)]);
penalty_factor = min([1 , sum(penalty_factor)]);

end

