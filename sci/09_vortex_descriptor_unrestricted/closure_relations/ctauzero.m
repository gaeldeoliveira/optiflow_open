function [ ctz, ctz_h , ctz_hk , ctz_hs , ctz_us ] = ctauzero( h , hk , hs , us)
%CTAUZERO Equilibrium shear coefficient 
%   (used in dissipation coefficient calculation, as forcing of shear lag equation!)
%   (comes from G-beta locus)

%     ---- Set G-beta locus according to IVW. (R.v.Rooij)
gacon  = 6.75;
gbcon  = 0.83;
ctcon  = (0.5/(gbcon.*gacon.^2));

%     ---- CTZ is the equilibrium shear stress without the contribution of suction
hkb = hk - 1.0;
usb = 1.0 - us;
ctz = ctcon .* hs .* hkb.^3  ./ (usb .* h .* hk.^2);

%     ---- And its derivatives
ctz_hs  =  ctcon .*       hkb.^3 ./ (usb .* h .* hk.^2);
ctz_us  =  ctcon .* hs .* hkb.^3 ./ (usb .* h .* hk.^2) ./ usb;
ctz_hk  =  ctcon .* hs .* hkb.^2 ./ (usb .* h .* hk.^2) * 3.0 - ctcon .* hs .* hkb.^3 ./ (usb*h*hk.^3) * 2.0;
ctz_h   = -ctcon .* hs .* hkb.^3 ./ (usb .* h .* hk.^2) ./ h;


end



% %% Original Fortran Code
%       SUBROUTINE CTAUZERO(H , HK , HS , US, CTZ, CTZ_H , CTZ_HK , CTZ_HS , CTZ_US)
% !     Gael de Oliveira (with code from R.v.Rooij and Drela)
% 
% 
% !     ---- Set G-beta locus according to IVW. (R.v.Rooij)
%       GACON  = 6.75
%       GBCON  = 0.83
%       CTCON = (0.5/(GBCON*GACON**2))
% 
% 
% !     CTZ is the equilibrium shear stress without the contribution of suction
%       HKB = HK - 1.0
%       USB = 1.0 - US
%       CTZ     = CTCON*HS*HKB**3 / (USB*H*HK**2)
% 
% !     And its derivatives
%       CTZ_HS = CTCON  *  HKB**3 / (USB*H*HK**2)
%       CTZ_US = CTCON*HS*HKB**3 / (USB*H*HK**2) / USB
%       CTZ_HK = CTCON*HS*HKB**2 / (USB*H*HK**2) * 3.0 - CTCON*HS*HKB**3 / (USB*H*HK**3) * 2.0
%       CTZ_H  =-CTCON*HS*HKB**3 / (USB*H*HK**2) / H
% 
%       RETURN
%       END