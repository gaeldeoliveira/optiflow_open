function [CF, CF_HK, CF_RT , CF_MSQ] = cft_ml_seq( HK, RT, MSQ, SPR)
     
     % % Extract data from SPR  
     TA1C = SPR.TA1C;
     TA2C = SPR.TA2C;
     TA3C = SPR.TA3C;
     TA4C = SPR.TA4C;
     TA5C = SPR.TA5C;
     TA6C = SPR.TA6C;
     THMIN= SPR.THMIN;
     THMAX= SPR.THMAX;
     TDCF = SPR.TDCF ;
     
     % % Now move on to fortran clone

      % SUBROUTINE CFTGMERCHANT( HK, RT, MSQ, US, SV, CF, CF_HK, CF_RT , CF_MSQ, CF_US, CF_SV)
      % IMPLICIT REAL (A-H,M,O-Z)
      % INCLUDE 'PLASMA_OPTIONS.f90'
      
      % CALL CFTGMERCHANT2( HK, RT, MSQ, US, SV, CF0, CF0_HK, CF0_RT , CF0_MSQ, CF0_US, CF0_SV)
      [ CF0, CF0_HK, CF0_RT, CF0_MSQ  ] = cft_seq( HK, RT, MSQ);
      
% !      CF     = CF0
% !      CF_HK  = CF0_HK
% !      CF_RT  = CF0_RT
% !      CF_MSQ = CF0_MSQ
% !      CF_US  = CF0_US
% !      CF_SV  = CF0_SV
      
      % CALL XD5_FROM_H(HK, THMIN, THMAX, XD5, XD5_H)
      [XD5, XD5_H] = xd5_from_h(HK, THMIN, THMAX);
      
      % CALL BST_D5( XD5, TA1C, TA2C, TA3C, TA4C, TA5C, TA6C, BD5, BD5_XD5)
      [BD5, BD5_XD5] = bst_d5(XD5, TA1C, TA2C, TA3C, TA4C, TA5C, TA6C);
      
      % CALL SD5_FROM_H(BD5, BD5_XD5, XD5_H, SD5, SD5_H)
      [SD5, SD5_H] = sd5_from_h(BD5, BD5_XD5, XD5_H);
      
      
      CF     = SD5   * (CF0 + TDCF) - TDCF;
      CF_HK  = SD5_H * (CF0 + TDCF) + SD5  * CF0_HK;
      CF_RT  =                        SD5  * CF0_RT;
      CF_MSQ =                        SD5  * CF0_MSQ;
      % CF_US  =                        SD5  * CF0_US
      % CF_SV  =                        SD5  * CF0_SV
      
% !      WRITE(*,1000) TA1C, TA2C, TA3C
% !      WRITE(*,1010) TA4C, TA5C, TA6C
% !      WRITE(*,1020) HK, THMIN, THMAX
% !      WRITE(*,1030) XD5, BD5, BD5_XD5
% !      WRITE(*,1040) XD5_H, SD5, SD5_H
      
% ! 1000 FORMAT(' TA1C=', F12.4, '   TA2C=', F12.4, '    TA3C=', F12.4)
% ! 1010 FORMAT(' TA4C=', F12.4, '   TA5C=', F12.4, '    TA6C=', F12.4)
% ! 1020 FORMAT('    H=', F12.4, '  THMIN=', F12.4, '   THMAX=', F12.4)
% ! 1030 FORMAT('  XD5=', F12.4, '     D5=', F12.4, '  D5_XD5=', F12.4)
% ! 1040 FORMAT('XD5_H=', F12.4, '    SD5=', F12.4, '   SD5_H=', F12.4)
 
      % RETURN
      % END
      
end

% For later reference 
%       SUBROUTINE CFTGMERCHANT( HK, RT, MSQ, US, SV, CF, CF_HK, CF_RT , CF_MSQ, CF_US, CF_SV)
%       IMPLICIT REAL (A-H,M,O-Z)
%       INCLUDE 'PLASMA_OPTIONS.f90'
%       
%       CALL CFTGMERCHANT2( HK, RT, MSQ, US, SV, CF0, CF0_HK, CF0_RT , CF0_MSQ, CF0_US, CF0_SV)
%       
% !      CF     = CF0
% !      CF_HK  = CF0_HK
% !      CF_RT  = CF0_RT
% !      CF_MSQ = CF0_MSQ
% !      CF_US  = CF0_US
% !      CF_SV  = CF0_SV
%       
%       CALL XD5_FROM_H(HK, THMIN, THMAX, XD5, XD5_H)
%       CALL BST_D5( XD5, TA1C, TA2C, TA3C, TA4C, TA5C, TA6C, BD5, BD5_XD5)
%       CALL SD5_FROM_H(BD5, BD5_XD5, XD5_H, SD5, SD5_H)
%       
%       CF     = SD5   * (CF0 + TDCF) - TDCF
%       CF_HK  = SD5_H * (CF0 + TDCF) + SD5  * CF0_HK
%       CF_RT  =                        SD5  * CF0_RT
%       CF_MSQ =                        SD5  * CF0_MSQ
%       CF_US  =                        SD5  * CF0_US
%       CF_SV  =                        SD5  * CF0_SV
%       
% !      WRITE(*,1000) TA1C, TA2C, TA3C
% !      WRITE(*,1010) TA4C, TA5C, TA6C
% !      WRITE(*,1020) HK, THMIN, THMAX
% !      WRITE(*,1030) XD5, BD5, BD5_XD5
% !      WRITE(*,1040) XD5_H, SD5, SD5_H
%       
% ! 1000 FORMAT(' TA1C=', F12.4, '   TA2C=', F12.4, '    TA3C=', F12.4)
% ! 1010 FORMAT(' TA4C=', F12.4, '   TA5C=', F12.4, '    TA6C=', F12.4)
% ! 1020 FORMAT('    H=', F12.4, '  THMIN=', F12.4, '   THMAX=', F12.4)
% ! 1030 FORMAT('  XD5=', F12.4, '     D5=', F12.4, '  D5_XD5=', F12.4)
% ! 1040 FORMAT('XD5_H=', F12.4, '    SD5=', F12.4, '   SD5_H=', F12.4)
%  
%       RETURN
%       END