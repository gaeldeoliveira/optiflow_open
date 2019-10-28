function [ cf, cf_hk, cf_rt, cf_msq  ] = cft( hk, rt, msq)
%CFT Turbulent skin friction function  ( Cf )    (Default Swafford in original)
%       rfoilsuc default is cft_rr (default changed by Gael)

    cf     = zeros(size(hk));
    cf_hk  = zeros(size(hk));
    cf_rt  = zeros(size(hk));
    cf_msq  = zeros(size(hk));
    for n_hk = 1:length(hk)
        [cf(n_hk), cf_hk(n_hk), cf_rt(n_hk), cf_msq(n_hk)] = cft_seq( hk(n_hk), rt, msq);
    end

end

% Original Fortran code
%       SUBROUTINE CFT( HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ )
%       REAL MSQ
% !
% !---- Turbulent skin friction function  ( Cf )    (Default Swafford)
%       FC = SQRT(1.0 + 0.2*MSQ)
%       GRT = ALOG(RT/FC)
%       GRT = AMAX1(GRT,3.0)
%       GEX = -1.74 - 0.31*HK
%       ARG = 4.0 - HK/0.875
%       ARG = AMIN1( 10.0, ARG )
%       ARG = AMAX1(-10.0, ARG )
%       
%       
%       CFO =  0.3*EXP(-1.33*HK) * (GRT/2.3026)**GEX
%       CF     = ( CFO  +  1.1E-4*(TANH(ARG)-1.0) ) / FC
%       CF_HK  = (-1.33*CFO - 0.31*ALOG(GRT/2.3026)*CFO &
%                - 1.1E-4/COSH(ARG)**2 / 0.875    ) / FC
%       CF_RT  = GEX*CFO/(FC*GRT) / RT
%       CF_MSQ = GEX*CFO/(FC*GRT) * (-0.1/FC**2)  -  0.1*CF/FC**2
% !
%       RETURN
%       END
