function [HS, HS_HK, HS_RT, HS_MSQ] = hst_ml_seq( HK, RT, MSQ, SPR)

     % % Extract data from SPR  
     TA1H = SPR.TA1H;
     TA2H = SPR.TA2H;
     TA3H = SPR.TA3H;
     TA4H = SPR.TA4H;
     TA5H = SPR.TA5H;
     TA6H = SPR.TA6H;
     THMIN= SPR.THMIN;
     THMAX= SPR.THMAX;
     
     % % Now move on to fortran clone

      %SUBROUTINE HST( HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ )
      %IMPLICIT REAL (A-H,M,O-Z)
      %INCLUDE 'PLASMA_OPTIONS.f90'
      
      % CALL HST2( HK, RT, MSQ, HS0, HS0_HK, HS0_RT, HS0_MSQ )
      [ HS0, HS0_HK, HS0_RT, HS0_MSQ] = hst_seq( HK, RT, MSQ);
      
%!      HS     = HS0
%!      HS_HK  = HS0_HK
%!      HS_RT  = HS0_RT
%!      HS_MSQ = HS0_MSQ
      
      % CALL XD5_FROM_H(HK, THMIN, THMAX, XD5, XD5_H)
      [XD5, XD5_H] = xd5_from_h(HK, THMIN, THMAX);
      
      % CALL BST_D5( XD5, TA1H, TA2H, TA3H, TA4H, TA5H, TA6H, BD5, BD5_XD5)
      [BD5, BD5_XD5] = bst_d5(XD5, TA1H, TA2H, TA3H, TA4H, TA5H, TA6H);
      
      % CALL SD5_FROM_H(BD5, BD5_XD5, XD5_H, SD5, SD5_H)
      [SD5, SD5_H] = sd5_from_h(BD5, BD5_XD5, XD5_H);
      
      HS     = SD5   * HS0                  ;
      HS_HK  = SD5_H * HS0 + SD5 * HS0_HK   ;
      HS_RT  =               SD5 * HS0_RT   ;
      HS_MSQ =               SD5 * HS0_MSQ  ;
      
%!      WRITE(*,1000) TA1C, TA2C, TA3C
%!      WRITE(*,1010) TA4C, TA5C, TA6C
%!      WRITE(*,1020) HK, THMIN, THMAX
%!      WRITE(*,1030) XD5, BD5, BD5_XD5
%!      WRITE(*,1040) XD5_H, SD5, SD5_H
      
%! 1000 FORMAT('Hstar TA1C=', F12.4, '   TA2C=', F12.4, '    TA3C=', F12.4)
%! 1010 FORMAT(' TA4C=', F12.4, '   TA5C=', F12.4, '    TA6C=', F12.4)
%! 1020 FORMAT('    H=', F12.4, '  THMIN=', F12.4, '   THMAX=', F12.4)
%! 1030 FORMAT('  XD5=', F12.4, '     D5=', F12.4, '  D5_XD5=', F12.4)
%! 1040 FORMAT('XD5_H=', F12.4, '    SD5=', F12.4, '   SD5_H=', F12.4)
      
%      RETURN
%      END

end
      
      
% For later reference
%       SUBROUTINE HST( HK, RT, MSQ, HS, HS_HK, HS_RT, HS_MSQ )
%       IMPLICIT REAL (A-H,M,O-Z)
%       INCLUDE 'PLASMA_OPTIONS.f90'
%       
%       CALL HST2( HK, RT, MSQ, HS0, HS0_HK, HS0_RT, HS0_MSQ )
%       
% !      HS     = HS0
% !      HS_HK  = HS0_HK
% !      HS_RT  = HS0_RT
% !      HS_MSQ = HS0_MSQ
%       
%       CALL XD5_FROM_H(HK, THMIN, THMAX, XD5, XD5_H)
%       CALL BST_D5( XD5, TA1H, TA2H, TA3H, TA4H, TA5H, TA6H, BD5, BD5_XD5)
%       CALL SD5_FROM_H(BD5, BD5_XD5, XD5_H, SD5, SD5_H)
%       
%       HS     = SD5   * HS0
%       HS_HK  = SD5_H * HS0 + SD5 * HS0_HK
%       HS_RT  =               SD5 * HS0_RT
%       HS_MSQ =               SD5 * HS0_MSQ
%       
% !      WRITE(*,1000) TA1C, TA2C, TA3C
% !      WRITE(*,1010) TA4C, TA5C, TA6C
% !      WRITE(*,1020) HK, THMIN, THMAX
% !      WRITE(*,1030) XD5, BD5, BD5_XD5
% !      WRITE(*,1040) XD5_H, SD5, SD5_H
%       
% ! 1000 FORMAT('Hstar TA1C=', F12.4, '   TA2C=', F12.4, '    TA3C=', F12.4)
% ! 1010 FORMAT(' TA4C=', F12.4, '   TA5C=', F12.4, '    TA6C=', F12.4)
% ! 1020 FORMAT('    H=', F12.4, '  THMIN=', F12.4, '   THMAX=', F12.4)
% ! 1030 FORMAT('  XD5=', F12.4, '     D5=', F12.4, '  D5_XD5=', F12.4)
% ! 1040 FORMAT('XD5_H=', F12.4, '    SD5=', F12.4, '   SD5_H=', F12.4)
%       
%       RETURN
%       END