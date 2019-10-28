function [ hh, hh_hk ] = hh1cal( hk )
%HH1CAL Head's shape parameter (H1 allows boundary layer thickness estimation)
%   Many possible versions here, the default one, used here is "Modified Green, after Drela"

hh      = 3.15 + 1.72 ./ (hk - 1);
hh_hk   =      - 1.72 ./ (hk - 1).^2;

end


% %% Stripped down Fortran Code 
% % (to follow execution path given goto, label and return keywords)
%       SUBROUTINE HH1CAL(HK,HH,HH_HK)
%       GOTO 200
% 
% !---- Modified Green, after Drela
%   200 HH      = 3.15 + 1.72/(HK - 1.)
%       HH_HK   =      - 1.72/( (HK-1)*(HK-1) )
% !
%       RETURN
% !
%       END


% %% Original Fortran Code
%       SUBROUTINE HH1CAL(HK,HH,HH_HK)
%       GOTO 200
% !
% !---- Head's shape parameter, Melnik approximation of Le Balleur
%       IF (HK .LT. 4.0) THEN
%         HH    = (0.5*HK + 1.0)*HK / (HK-1.0)
%         HH_HK = (0.5*HK*HK-HK-1.0) /( (HK-1.0)*(HK-1.0) )
%       ELSE
%         HH    = 1.75 + 5.52273*HK/(HK+5.818181)
%         HH_HK = 32.13224/( (HK+5.818181)**2. )
%       ENDIF
% !
%       RETURN
% !
% !---- Modified Le Balleur (after Houwink)
%   100 IF (HK .LT. 2.732) THEN
%         HH    = (0.5*HK + 1.0)*HK / (HK-1.0)
%         HH_HK = (0.5*HK*HK-HK-1.0) /( (HK-1.0)*(HK-1.0) )
%       ELSE
%         HS    = 0.5*(HK-2.732) + 2.732
%         IF (HS.LT.4.0) THEN
%           HH    = (0.5*HS + 1.0)*HS / (HS-1.0)
%           HH_HK = 0.5*(0.5*HS*HS-HS-1.0) /( (HS-1.0)*(HS-1.0) )
%         ELSE
%           HH    = 1.75 + 5.52273*HS/(HS+5.818181)
%           HH_HK = 16.06612/( (HS+5.818181)**2. )
%         ENDIF
%       ENDIF
% !
%       RETURN
% !
% !---- Modified Green, after Drela
%   200 HH      = 3.15 + 1.72/(HK - 1.)
%       HH_HK   =      - 1.72/( (HK-1)*(HK-1) )
% !
%       RETURN
% !
% !---- "Head's shape" parameter by HEAD
%   300 HH     = 1.535/((HK-.7)**2.715) + 3.3
%       HH_HK  = -4.1675/((HK-.7)**3.715)
% !
%       RETURN
% !
% !---- "Head's shape" parameter, Lock approximation
%   400 IF (HK .LT. 4.0) THEN
%         HH    = 2.+ 1.5*(1.12/(HK-1.))**1.093+ .5*((HK-1.)/1.12)**1.093
%         HH_HK = -1.8557/((HK-1.)**2.093) + .48278*(HK-1.)**.093
%       ELSE
%         HH    = 4. + (HK-4.)/3
%         HH_HK = - 1.0/3.0
%       ENDIF
% !
%       RETURN
% !
%       END
% 
