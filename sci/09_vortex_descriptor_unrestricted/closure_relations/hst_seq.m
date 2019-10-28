function [ hs, hs_hk, hs_rt, hs_msq  ] = hst_seq( hk, rt, msq)
%HST Turbulen HS correlation
%       new correlation  29 Nov 91
%       (from  arctan(y+) + Schlichting  profiles)
% Note: HSTIVW is currently not used by Rfoilsuc

%      DATA HSMIN, DHSINF / 1.500, 0.04  /
hsmin  = 1.5;
dhsinf = 0.04;

%       IF(RT.GT.400.0) THEN
%        h0    = 3.0 + 400.0/RT
%        h0_RT =     - 400.0/RT**2
%       ELSE
%        h0    = 4.0
%        h0_RT = 0.0
%       ENDIF

if rt > 400 
    h0      = 3 + 400 ./ rt;
    h0_rt   =   - 400 ./ rt.^2;         % This part had been translated wrong, and was responsible for transformation problems!
else
    h0      = 4.0;
    h0_rt   = 0.0;
end

if hk < h0
%       IF(HK.LT.h0) THEN
    % !----- attached branch

    hr    = ( h0 - hk) ./ (h0-1.0);
    hr_hk =      - 1.0 ./ (h0-1.0);
    hr_rt = (1.0 - hr) ./(h0-1.0) .* h0_rt;
    
    hs    = (2.0 - hsmin-4.0 ./rt) .*hr.^2  * 1.5./(hk+0.5)  +  hsmin+4.0./rt;
    hs_hk =-(2.0 - hsmin-4.0 ./rt) .*hr.^2  * 1.5./(hk+0.5)^2 + (2.0-hsmin-4.0./rt) .* hr *2.0 * 1.5 ./(hk+0.5) * hr_hk;
    hs_rt = (2.0 - hsmin-4.0 ./rt) .*hr * 2 * 1.5./(hk+0.5) * hr_rt + (hr.^2 * 1.5./(hk+0.5) - 1.0)*4.0./rt.^2;
else
    % !----- separated branch
    grt  = log(rt);
    hdif = hk - h0;
    rtmp = hk - h0 + 4.0 ./ grt;
    
    htmp    = 0.007 * grt./rtmp.^2  + dhsinf ./ hk;
    htmp_hk = -.014 * grt./rtmp.^3  - dhsinf ./hk.^2;
    htmp_rt = -.014 * grt./rtmp.^3 .* (-h0_rt - 4.0./(grt.^2 .* rt)) + 0.007 ./ (rtmp.^2 * rt);
    
    hs    = hdif.^2  .* htmp + hsmin + 4.0./rt;
    hs_hk = hdif*2.0 .* htmp + hdif.^2 .* htmp_hk;
    hs_rt = hdif.^2  .* htmp_rt      - 4.0/rt.^2 + hdif * 2.0 .* htmp .* (-h0_rt);
end

%!---- Whitfield's minor additional compressibility correction
fm     =  1.0 + 0.014 * msq;
hs     = ( hs + 0.028 * msq ) ./ fm;
hs_hk  = ( hs_hk            ) ./ fm;
hs_rt  = ( hs_rt            ) ./ fm;
hs_msq = 0.028 ./ fm  -  0.014 *hs ./ fm;


end


% Original Fortran Code
%     SUBROUTINE HST( HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ )
%       IMPLICIT REAL (A-H,M,O-Z)
% !
% !---- Turbulent HS correlation
% !
% !      DATA HSMIN, DHSINF / 1.500, 0.015 /
%       DATA HSMIN, DHSINF / 1.500, 0.04  /
% !
%       IF(RT.GT.400.0) THEN
%        h0    = 3.0 + 400.0/RT
%        h0_RT =     - 400.0/RT**2
%       ELSE
%        h0    = 4.0
%        h0_RT = 0.0
%       ENDIF
% !
%       IF(HK.LT.h0) THEN
% !----- new correlation  29 Nov 91
% !-     (from  arctan(y+) + Schlichting  profiles)
%        hr    = ( h0 - HK)/(h0-1.0)
%        HR_HK =      - 1.0/(h0-1.0)
%        HR_RT = (1.0 - hr)/(h0-1.0) * h0_RT
%        HS    = (2.0-HSMIN-4.0/RT)*HR**2  * 1.5/(HK+0.5)  +  HSMIN+4.0/RT
%        HS_HK =-(2.0-HSMIN-4.0/RT)*HR**2  * 1.5/(HK+0.5)**2 &
%              + (2.0-HSMIN-4.0/RT)*hr*2.0 * 1.5/(HK+0.5) * HR_HK
%        HS_RT = (2.0-HSMIN-4.0/RT)*HR*2.0 * 1.5/(HK+0.5) * HR_RT &
%              + (HR**2 * 1.5/(HK+0.5) - 1.0)*4.0/RT**2
% !
%       ELSE
% !
% !----- separated branch
%        GRT = ALOG(RT)
%        hdif = HK - h0
%        rtmp = HK - h0 + 4.0/GRT
%        HTMP    = 0.007*GRT/RTMP**2 + DHSINF/HK
%        HTMP_HK = -.014*GRT/RTMP**3 - DHSINF/HK**2
%        HTMP_RT = -.014*GRT/RTMP**3 * (-h0_RT - 4.0/GRT**2 / RT) &
%                + 0.007    /RTMP**2 / RT
%        HS    = HDIF**2 * HTMP + HSMIN + 4.0/RT
%        HS_HK = HDIF*2.0* HTMP &
%              + HDIF**2 * HTMP_HK
%        HS_RT = HDIF**2 * HTMP_RT      - 4.0/RT**2 &
%              + hdif*2.0* HTMP * (-h0_RT)
% !
%       ENDIF
% !
% !---- Whitfield's minor additional compressibility correction
%       FM = 1.0 + 0.014*MSQ
%       HS     = ( HS + 0.028*MSQ ) / FM
%       HS_HK  = ( HS_HK          ) / FM
%       HS_RT  = ( HS_RT          ) / FM
%       HS_MSQ = 0.028/FM  -  0.014*HS/FM
% !
%       RETURN
%       END
% 
