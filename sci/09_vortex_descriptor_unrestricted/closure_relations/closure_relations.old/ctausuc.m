function [ cts, cts_h , cts_hk , cts_hs , cts_us, cts_sv ] = ctausuc( h , hk , hs , us, sv)
%CTAUSUC    CTS is contribution of suction to the equilibrium shear stress (Gael)

%       CTS is contribution of suction to the equilibrium shear stress (Gael)
cts = sv .* 0.5 .* (hs + hs .* (1 - hk) ./ h - 1) ./ (1 - us);

%       And its derivatives
cts_sv  =           0.5  .* (hs + hs .* (1 - hk) ./ h - 1.0) ./ (1 - us);
cts_hk  = - sv.^2 * 0.5  .*        1 ./ h                    ./ (1 - us);
cts_h   =      sv * 0.5 .* hs  .* (hk - 1) .* (1./h.^2)      ./ (1 - us);
cts_hs  =      sv * 0.5 .*  (1 + (1 - hk) ./ h)              ./ (1 - us);
cts_us  =                                                cts ./ (1 - us);


end

% %% Original Fortran Code
%       SUBROUTINE CTAUSUC(H , HK , HS , US, SV, CTS, CTS_H , CTS_HK , CTS_HS , CTS_US, CTS_SV)
% !      Gael de Oliveira
% !      CTS is contribution of suction to the equilibrium shear stress
%       CTS = SV * 0.5 * (HS + HS * (1.0 - HK) / H - 1.0) / (1.0-US)
% 
% !      And its derivatives
%       CTS_SV = 0.5 * (HS + HS * (1.0 - HK) / H - 1.0) / (1.0-US)
%       CTS_HK = - SV * 0.5 * SV * (1.0/H) / (1.0-US)
%       CTS_H = SV * 0.5 * HS * (HK - 1.0) * (1.0/H**2) / (1.0-US)
%       CTS_HS = SV * 0.5 * (1.0 + (1.0 - HK) / H) / (1.0 - US)
%       CTS_US = CTS / (1.0-US)
% 
%       RETURN
%       END

