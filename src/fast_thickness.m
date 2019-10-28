function [th, bh] = fast_thickness(results)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
tx = results.coordinates.tx_analysis;
tz = results.coordinates.tz_analysis;

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

% Compute thickness
th = max(tz_top) - min(tz_bot);
% And maximum building height
bh = max(tz_top_fun(tx_top) - tz_bot_fun(tx_top));
end








