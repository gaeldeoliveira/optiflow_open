function [ cf, cf_hk, cf_rt, cf_msq  ] = cft_rr( hk, rt, msq)
%CFT_RR     Turbulent skin friction function  ( Cf )    (Higher Cf at low Hk)
%       Default on new rfoilsuc version (default choice changed by gael)

fc = sqrt(1.0 + 0.2*msq);
grt = log(rt./fc);
grt = max(grt,3.0);
gex = -1.75 - 0.25 * hk;
arg = 4.0 - hk/0.5;
arg = min( 10.0, arg );
arg = max(-10.0, arg );

cf0     =  0.325*exp(-1.37*hk) .* (grt/2.3026).^gex;
cf      = ( cf0 +  1.0E-4*(tanh(arg)-1.0) ) ./ fc;

cf_hk  = (-1.37*cf0 - 0.25 * log(grt/2.3026) .* cf0 - 1.0E-4./cosh(arg).^2 / 0.5 ) ./ fc;

cf_rt  = gex*cf0/(fc*grt) ./ rt;
cf_msq = gex*cf0/(fc*grt) .* (-0.1./fc.^2)  -  0.1*cf./fc.^2;

end

% Original Fortran Code
%       SUBROUTINE CFT_RR( HK, RT, MSQ, CF, CF_HK, CF_RT, CF_MSQ )
%       REAL MSQ
% !
% !---- Turbulent skin friction function  ( Cf )    (Higher Cf at low Hk)
%       FC = SQRT(1.0 + 0.2*MSQ)
%       GRT = ALOG(RT/FC)
%       GRT = AMAX1(GRT,3.0)
%       GEX = -1.75 - 0.25*HK
%       ARG = 4.0 - HK/0.5
%       ARG = AMIN1( 10.0, ARG )
%       ARG = AMAX1(-10.0, ARG )
%       
%       CFO =  0.325*EXP(-1.37*HK) * (GRT/2.3026)**GEX
%       CF     = ( CFO  +  1.0E-4*(TANH(ARG)-1.0) ) / FC
%       CF_HK  = (-1.37*CFO - 0.25*ALOG(GRT/2.3026)*CFO &
%                - 1.0E-4/COSH(ARG)**2 / 0.5    ) / FC
%       CF_RT  = GEX*CFO/(FC*GRT) / RT
%       CF_MSQ = GEX*CFO/(FC*GRT) * (-0.1/FC**2)  -  0.1*CF/FC**2
% !
%       RETURN
%       END