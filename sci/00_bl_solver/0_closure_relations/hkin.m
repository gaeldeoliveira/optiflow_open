function [ hk, hk_h, hk_msq  ] = hkin( h, msq)
%HKIN Calculate kinematic shape parameter (assuming air)
% !     (from Whitfield )

hk     = (h - 0.29*msq)     ./  (1.0 + 0.113*msq);
hk_h   =  1.0               ./  (1.0 + 0.113*msq);
hk_msq = (-.29 - 0.113*hk)  ./  (1.0 + 0.113*msq);

end


% %% Original Fortran Code
%       SUBROUTINE HKIN( H, MSQ, HK, HK_H, HK_MSQ )
%       REAL MSQ
% !
% !---- calculate kinematic shape parameter (assuming air)
% !     (from Whitfield )
%       HK     = (H - 0.29*MSQ)/(1.0 + 0.113*MSQ)
%       HK_H   =  1.0          /(1.0 + 0.113*MSQ)
%       HK_MSQ = (-.29 - 0.113*HK) / (1.0 + 0.113*MSQ)
% !
%       RETURN
%       END



