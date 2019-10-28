function [ cf, cf_hk, cf_rt, cf_msq  ] = cft_seq( hk, rt, msq)
%CFT Turbulent skin friction function  ( Cf )    (Default Swafford in original)
%       rfoilsuc default is cft_rr (default changed by Gael)

fc  = sqrt(1.0 + 0.2 * msq);
grt = log(rt./fc);
grt = max(grt,3.0);
gex = -1.74 - 0.31 * hk;
arg = 4.0 - hk / 0.875;
arg = min( 10.0, arg );
arg = max(-10.0, arg );

cf0     =  0.3 * exp(-1.33 * hk) .* (grt/2.3026) .^gex;
cf      = ( cf0                                    + 1.1E-4*(tanh(arg)-1.0) )       / fc;
cf_hk   = (-1.33*cf0 - 0.31*log(grt/2.3026) .* cf0 - 1.1E-4./cosh(arg).^2 / 0.875 ) / fc;
cf_rt   = gex .* cf0 ./ (fc .* grt .* rt);
cf_msq  = gex .* cf0 ./ (fc*grt) .* (-0.1./fc.^2)  -  0.1*cf ./ fc.^2;

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
