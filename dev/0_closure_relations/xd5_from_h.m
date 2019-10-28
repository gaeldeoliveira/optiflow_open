function [XD5, XD5_H] = xd5_from_h(H, THMIN, THMAX)

      % SUBROUTINE XD5_FROM_H(H, THMIN, THMAX, XD5, XD5_H)
      % IMPLICIT REAL (A-Z) 
      % ! Map H intervention region into unit interval
      XD5   = (H - THMIN) ./ (THMAX - THMIN);
      XD5_H =  1.0        ./ (THMAX - THMIN);
      % ! Fudge for values of H above intervention region
      if (H > THMAX)
             XD5   = 1.0;
             XD5_H = 0.0;
      end
      % ! Fudge for values of H below intervention region
      if (H < THMIN)
             XD5   = 0.0;
             XD5_H = 0.0;
      end
      % RETURN
      % END
      
end
      
      
      
% For later reference
%       SUBROUTINE XD5_FROM_H(H, THMIN, THMAX, XD5, XD5_H)
%       IMPLICIT REAL (A-Z) 
%       ! Map H intervention region into unit interval
%       XD5   = (H - THMIN) / (THMAX - THMIN)
%       XD5_H =  1.0        / (THMAX - THMIN)
%       ! Fudge for values of H above intervention region
%       IF(H .GT. THMAX) THEN
%              XD5   = 1.0
%              XD5_H = 0.0
%       ENDIF
%       ! Fudge for values of H below intervention region
%       IF(H .LT. THMIN) THEN
%              XD5   = 0.0
%              XD5_H = 0.0
%       ENDIF
%       RETURN
%       END