function [ hc, hc_hk, hc_msq ] = hct( hk, msq )
%HCT Density Shape Parameter (from Whitfield)
%   Detailed explanation goes here

hc     = msq .* (0.064 ./ (hk-0.8) + 0.251);
hc_hk  = msq .* (-.064 ./ (hk-0.8).^2     );
hc_msq =         0.064 ./ (hk-0.8) + 0.251 ;

end

% Original Fortran Code
%       SUBROUTINE HCT( HK, MSQ, HC, HC_HK, HC_MSQ )
%       REAL MSQ
% !
% !---- density shape parameter    (from Whitfield)
%       HC     = MSQ * (0.064/(HK-0.8) + 0.251)
%       HC_HK  = MSQ * (-.064/(HK-0.8)**2     )
%       HC_MSQ =        0.064/(HK-0.8) + 0.251
% !
%       RETURN
%       END