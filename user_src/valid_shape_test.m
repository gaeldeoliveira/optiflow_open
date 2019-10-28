function valid_shape = valid_shape_test(experiment_result)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
tx = experiment_result.coordinates.tx_analysis;
tz = experiment_result.coordinates.tz_analysis;

% Find leading edge
[~ , ile] = min(tx);

% Extract sided tx coordinates
tx_top = tx(ile:(-1):1  );
tx_bot = tx(ile:(+1):end);
% Extract sided tz coordinates
tz_top = tz(ile:(-1):1  );
tz_bot = tz(ile:(+1):end);

% Make interpolants
tz_top_fun = @(tx) interp1(tx_top / tx_top(end), tz_top, tx);
tz_bot_fun = @(tx) interp1(tx_bot / tx_bot(end), tz_bot, tx);

% Make tespoint vector (interpolant needed because two sides may not have
% same discretization in all cases! use top discretization as reference!)
tx_test = tx_top;

% And proceed to test crossing condition
dtz_test = tz_top_fun(tx_test) - tz_bot_fun(tx_test);

% Check if it is met (true = valid shape, false = invalid shape)
% valid_shape = isempty(find(dtz_test < 0));
valid_shape = isempty(find(dtz_test < 0, 1));

end








