function [ cq, cq_h , cq_hk , cq_hs , cq_us, cq_sv ] = cqt( h , hk , hs , us, sv )
%CQT is the square root of the total shear stress coefficient CTT, composed of:
%         A term due to suction CTS
%         A term due to the rest CTZ
%      CQ = equilibrium wake layer shear coefficient (Ctau)EQ ** 1/2

%      Start by calling auxilliary functions for each contribution to shear stress
%       CALL CTAUZERO(H , HK , HS , US, CTZ, CTZ_H , CTZ_HK , CTZ_HS , CTZ_US)
%       CALL CTAUSUC(H , HK , HS , US, SV, CTS, CTS_H , CTS_HK , CTS_HS , CTS_US, CTS_SV)
[ ctz, ctz_h , ctz_hk , ctz_hs , ctz_us ] = ctauzero( h , hk , hs , us);
[ cts, cts_h , cts_hk , cts_hs , cts_us, cts_sv ] = ctausuc( h , hk , hs , us, sv);

%      Now compose Maximum shear stress and take square root
ctt     = cts + ctz;
cq      = sqrt(ctt);

%      Now proceed proceed to derivatives with handy identity derived in Gael's master thesis
cq_h  = (0.5 ./ cq) .* (cts_h  + ctz_h );
cq_hk = (0.5 ./ cq) .* (cts_hk + ctz_hk);
cq_hs = (0.5 ./ cq) .* (cts_hs + ctz_hs);
cq_us = (0.5 ./ cq) .* (cts_us + ctz_us);
cq_sv = (0.5 ./ cq) .* (cts_sv);

end

% %% Original Fortran Code
%       SUBROUTINE CQT(H , HK , HS , US, SV, CQ, CQ_H , CQ_HK , CQ_HS , CQ_US, CQ_SV)
% !      Gael de Oliveira
% !      CQT is the square root of the total shear stress coefficient CTT, composed of:
% !         A term due to suction CTS
% !         A term due to the rest CTZ
% !      CQ = equilibrium wake layer shear coefficient (Ctau)EQ ** 1/2
% 
% !      Start by calling auxilliary functions for each contribution to shear stress
%       CALL CTAUZERO(H , HK , HS , US, CTZ, CTZ_H , CTZ_HK , CTZ_HS , CTZ_US)
%       CALL CTAUSUC(H , HK , HS , US, SV, CTS, CTS_H , CTS_HK , CTS_HS , CTS_US, CTS_SV)
% 
% !      Now compose Maximum shear stress and take square root
%       CTT = CTS + CTZ
%       CQ = SQRT(CTT)
% 
% !      Now proceed proceed to derivatives with handy identity we derived in report
%       CQ_H  = (0.5 / CQ) * (CTS_H + CTZ_H)
%       CQ_HK = (0.5 / CQ) * (CTS_HK + CTZ_HK)
%       CQ_HS = (0.5 / CQ) * (CTS_HS + CTZ_HS)
%       CQ_US = (0.5 / CQ) * (CTS_US + CTZ_US)
%       CQ_SV = (0.5 / CQ) * (CTS_SV)
% 
%       RETURN
%       END
