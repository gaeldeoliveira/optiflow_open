function az = nan_to_zero(a)
% Fudges all NaNs in a vector to zero and returns that vector but without
% the NaNs
    a(isnan(a)) = 0; 
    az = a;
end