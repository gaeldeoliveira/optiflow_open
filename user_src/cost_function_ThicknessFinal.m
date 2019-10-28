function sc = cost_function_ThicknessFinal( ~ , experiments_results)
% Calculated weighted Cl/Cd of a polar at AOA/Weight combinations from
% alpha_i and w_i vectors

% Input extraction
% parameters can be ignored

% Extract aerodynamic polar from standard experiment input form
experiment_result   = experiments_results{1};
% parameters = experiment_result.x;
coordinates = experiment_result.coordinates;


NumPoints = size(coordinates.tx,1);

Ztop = coordinates.tz(1:ceil(NumPoints/2));
Zbottom = coordinates.tz(end:-1:ceil(NumPoints/2));

Thickness = Ztop - Zbottom;


[MaxThickness ind]= max(Thickness);
Location = coordinates.tx(ind);

    spaces10 = '          ';
    Ttext = num2str(round(MaxThickness*100));
    LocText = num2str(round(Location*100));

   % Display Information on Selection
%    disp(['     Tmax = ' num2str(round(MaxThickness*100)) '%    @ ' num2str(round(Location*100)) '% of the chord'])
   disp(['     Tmax = ' Ttext '%' spaces10(1:end-1-length(Ttext)) '@  ' LocText '%' spaces10(1:4-length(LocText))  'of the chord'])
 

   sc = MaxThickness;
end