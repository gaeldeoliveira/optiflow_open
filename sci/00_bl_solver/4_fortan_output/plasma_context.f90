! -----------------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------------
!
!     Gael de Oliveira : Closure functions for flow control devices (Xfoil compatible)
!
! -----------------------------------------------------------------------------------------
! -----------------------------------------------------------------------------------------

      SUBROUTINE USG(H2 , HK2 , HS2 , US2, US2_H2 , US2_HK2 , US2_HS2)
      IMPLICIT REAL (A-Z)
!     ---- Gael de Oliveira 	(Rededuction of Drela for generic B)
      US2     = 0.5*HS2*( 3.0 - 4.0*(HK2-1.0)/H2   )/3.0
      US2_HS2 = 0.5  *  ( 3.0 - 4.0*(HK2-1.0)/H2   )/3.0
      US2_HK2 = 0.5*HS2*(     - 4.0          /H2   )/3.0
      US2_H2  = 0.5*HS2*(       4.0*(HK2-1.0)/H2**2)/3.0
      
      !     ---- Xfoil  	(updated, from Drela)      
      !      US2     = 0.5*HS2*( 1.0 - (HK2-1.0)/(GBCON*H2) )
      !      US2_HS2 = 0.5  *  ( 1.0 - (HK2-1.0)/(GBCON*H2) )
      !      US2_HK2 = 0.5*HS2*(     -  1.0     /(GBCON*H2) )
      !      US2_H2  = 0.5*HS2*        (HK2-1.0)/(GBCON*H2**2)
      
      RETURN
      END

      SUBROUTINE CFTCMERCHANT( US , SV, CFTC, CFTC_US, CFTC_SV)
      IMPLICIT REAL (A-Z)
!      CFTC is value of correction to CF for accounting for transpiration according to merchant
!      Gael de Oliveira
      CFTC = - 2.0 * SV * US
      CFTC_US = - 2.0 * SV
      CFTC_SV = - 2.0 * US
!      CFTC_US = 0.
!      CFTC_SV = 0.
      RETURN
      END


      SUBROUTINE CTAUSUC(H , HK , HS , US, SV, CTS, CTS_H , CTS_HK , CTS_HS , CTS_US, CTS_SV)
!      Gael de Oliveira
!      CTS is contribution of suction to the equilibrium shear stress
      CTS = SV * 0.5 * (HS + HS * (1.0 - HK) / H - 1.0) / (1.0-US)

!      And its derivatives
      CTS_SV = 0.5 * (HS + HS * (1.0 - HK) / H - 1.0) / (1.0-US)
      CTS_HK = - SV * 0.5 * SV * (1.0/H) / (1.0-US)
      CTS_H = SV * 0.5 * HS * (HK - 1.0) * (1.0/H**2) / (1.0-US)
      CTS_HS = SV * 0.5 * (1.0 + (1.0 - HK) / H) / (1.0 - US)
      CTS_US = CTS / (1.0-US)

      RETURN
      END


      SUBROUTINE CQT(H , HK , HS , US, SV, CQ, CQ_H , CQ_HK , CQ_HS , CQ_US, CQ_SV)
!      Gael de Oliveira
!      CQT is the square root of the total shear stress coefficient CTT, composed of:
!         A term due to suction CTS
!         A term due to the rest CTZ
!      CQ = equilibrium wake layer shear coefficient (Ctau)EQ ** 1/2

!      Start by calling auxilliary functions for each contribution to shear stress
      CALL CTAUZERO(H , HK , HS , US, CTZ, CTZ_H , CTZ_HK , CTZ_HS , CTZ_US)
      CALL CTAUSUC(H , HK , HS , US, SV, CTS, CTS_H , CTS_HK , CTS_HS , CTS_US, CTS_SV)

!      Now compose Maximum shear stress and take square root
      CTT = CTS + CTZ
      CQ = SQRT(CTT)

!      Now proceed proceed to derivatives with handy identity we derived in report
      CQ_H  = (0.5 / CQ) * (CTS_H + CTZ_H)
      CQ_HK = (0.5 / CQ) * (CTS_HK + CTZ_HK)
      CQ_HS = (0.5 / CQ) * (CTS_HS + CTZ_HS)
      CQ_US = (0.5 / CQ) * (CTS_US + CTZ_US)
      CQ_SV = (0.5 / CQ) * (CTS_SV)

      RETURN
      END


      SUBROUTINE VSUCBLSET_BACKWARDS(VSUCTION, IBLS, ISS)
      INCLUDE 'RFOIL.f90'
      INTEGER IBLS, ISS
      REAL VSUCTION
      ! Use of auxiliary function is to avoid potential namespace conflicts in BLVAR
      ! This a tweaked approach... but seemed a good idea...

      IBLTETMP = IBLTE(ISS)
      ! Check that we are before trailing edge (that is not in wake)
      IF (IBLS .LT. IBLTETMP) THEN
      ! Set VSUCBL in BL indexing scheme
        VSUCBL(IBLS,ISS) = VSUCTION
        ! Set VSUC in panel indexing scheme to avoid loss of data in VSUCBLSET
        I = IPAN(IBLS,ISS)
        VSUC(I) = VSUCTION
      ENDIF

      RETURN
      END

      SUBROUTINE CFTCKAYS( CFI , SV, CFTC, CFTC_CFI, CFTC_SV)
      IMPLICIT REAL (A-Z)
!      CFTC is value of correction to CF for accounting for transpiration according to merchant
!      Gael de Oliveira

      ABSSV = ABS(SV)
      ! Suction is large enough for calculation to be stable
      IF (ABSSV .GE. 0.0001) THEN
        ! Before separation
        IF (CFI .GT. 0.00001) THEN
          ! Suction Parameter
          BFI     =  2 * SV /  CFI
!          BFI_CFI = -2 * SV /  CFI**2

          CFTC     = 2 * SV / (exp(BFI) -1)
!          CFTC_BF = - exp(BFI) /  (exp(BFI) -1)**2
!          CFTC_SV  = 2 / (exp(BFI) -1)**2 + SV * CFTC_BF
          CFTC_SV  = 2 / (exp(BFI) -1)**2 * ( exp(BFI) - 1 - SV * (2 / CFI) * exp(BFI))
          CFTC_CFI = BFI**2 * exp(BFI) / (exp(BFI) -1)**2

        ENDIF

        ! Handle post separation
        IF (CFI .LE. 0.00001) THEN
          CFTC     = CFI
          CFTC_CFI = 1
          CFTC_SV  = 0.0
        END IF
      ENDIF

      ! Suction is too small for numerical calculation
      IF (ABSSV .LT. 0.0001) THEN
        CFTC     = CFI  ! Approximation
        CFTC_CFI = 1
!        CFTC_SV  = -1   ! Analytical Limit
        CFTC_SV  = 0.0
      ENDIF

      RETURN
      END

! mlcode -->
!     du25.air re3e6
!     6   rms: 0.7186E-05   max: -.1280E-03   D at   49  1
!     a =  0.000       CL =  0.4415
!     Cm = -0.1306      CD =  0.00720     CDp =  0.00246
      SUBROUTINE CFTGMERCHANT( HK, RT, MSQ, US, SV, CF, CF_HK, CF_RT , CF_MSQ, CF_US, CF_SV)
      IMPLICIT REAL (A-H,M,O-Z)
      INCLUDE 'PLASMA_OPTIONS.f90'
      
      CALL CFTGMERCHANT2( HK, RT, MSQ, US, SV, CF0, CF0_HK, CF0_RT , CF0_MSQ, CF0_US, CF0_SV)
      
!      CF     = CF0
!      CF_HK  = CF0_HK
!      CF_RT  = CF0_RT
!      CF_MSQ = CF0_MSQ
!      CF_US  = CF0_US
!      CF_SV  = CF0_SV
      
      CALL XD5_FROM_H(HK, THMIN, THMAX, XD5, XD5_H)
      CALL BST_D5( XD5, TA1C, TA2C, TA3C, TA4C, TA5C, TA6C, BD5, BD5_XD5)
      CALL SD5_FROM_H(BD5, BD5_XD5, XD5_H, SD5, SD5_H)
      
      CF     = SD5   * (CF0 + TDCF) - TDCF
      CF_HK  = SD5_H * (CF0 + TDCF) + SD5  * CF0_HK
      CF_RT  =                        SD5  * CF0_RT
      CF_MSQ =                        SD5  * CF0_MSQ
      CF_US  =                        SD5  * CF0_US
      CF_SV  =                        SD5  * CF0_SV
      
!      WRITE(*,1000) TA1C, TA2C, TA3C
!      WRITE(*,1010) TA4C, TA5C, TA6C
!      WRITE(*,1020) HK, THMIN, THMAX
!      WRITE(*,1030) XD5, BD5, BD5_XD5
!      WRITE(*,1040) XD5_H, SD5, SD5_H
      
! 1000 FORMAT(' TA1C=', F12.4, '   TA2C=', F12.4, '    TA3C=', F12.4)
! 1010 FORMAT(' TA4C=', F12.4, '   TA5C=', F12.4, '    TA6C=', F12.4)
! 1020 FORMAT('    H=', F12.4, '  THMIN=', F12.4, '   THMAX=', F12.4)
! 1030 FORMAT('  XD5=', F12.4, '     D5=', F12.4, '  D5_XD5=', F12.4)
! 1040 FORMAT('XD5_H=', F12.4, '    SD5=', F12.4, '   SD5_H=', F12.4)
 
      RETURN
      END
! <-- mlcode

      SUBROUTINE CFTGMERCHANT2( HK, RT, MSQ, US, SV, CF, CF_HK, CF_RT , CF_MSQ, CF_US, CF_SV)
      IMPLICIT REAL (A-H,M,O-Z)
!     Gael de Oliveira
!      (...)I represents impermeable contribution to value (...)
!     (...)P represents permeable contribution to value (...)
!      Make suction options available
       INCLUDE 'RBL_OPTIONS.f90'
!      Call (and prepare to) original Cf functions without suction
!      Rfoil specific code
!        First access the flags for choice of Cf function
      common /TURB_CHOICE/ L_SCC,L_CFT
      logical L_SCC,L_CFT
!        Second: call proper Cf function
       if (L_CFT) then
         CALL CFT_RR( HK, RT, MSQ, CFI, CFI_HK, CFI_RT, CFI_MSQ )
        else
         CALL CFT(    HK, RT, MSQ, CFI, CFI_HK, CFI_RT, CFI_MSQ )
        endif
	 
! --> roughnesscode
!     Get (with UE/UINF innacuracy) the roughness height reynolds number
      TREH = THCR
!     Compute skin friction increase factor due to surface roughness
      CALL CECT(HK , RT, TREH, CEC, CEC_HK, CEC_RT, CEC_TREH)
!     Apply skin friction increase (on skin friction not yet affected by suction)
      CFI     = CFI     * (1.0 + CEC)
      CFI_HK  = CFI_HK  * (1.0 + CEC)  +  CFI * (0.0 + CEC_HK )
      CFI_RT  = CFI_RT  * (1.0 + CEC)  +  CFI * (0.0 + CEC_RT )
      CFI_MSQ = CFI_MSQ * (1.0 + CEC)  +  CFI * (0.0 + CEC_MSQ)
!     Temporary diagnostic display
!      WRITE (*,*) 'CFI=', CFI, '  CEC=', CEC, '  CEC_HK=', CEC_HK
! <-- roughnesscode

!      Equivalent Xfoil specific code
!      CALL CFT( HK, RT, MSQ, CFI, CFI_HK, CFI_RT, CFI_MSQ )
      IF(.NOT. L_CF_KAYS) THEN
  !     And now suction contribution to Cf
        CALL CFTCMERCHANT( US , SV, CFTC, CFTC_US, CFTC_SV)
        ! WRITE (*,*) 'Merchant Cf used'
        ! Ratio_Cfsuc is used to activate or deactive correction of cf due to suction
        CFTC    = CFTC
        CFTC_HK = 0.0
        CFTC_RT = 0.0
        CFTC_MSQ = 0.0
        CFTC_US = CFTC_US
        CFTC_SV = CFTC_SV

        !       Now Compose Cf expression
        CF     = CFI      + RATIO_CFSUC * CFTC
!       And its derivatives
        CF_HK  = CFI_HK   + RATIO_CFSUC * CFTC_HK
        CF_RT  = CFI_RT   + RATIO_CFSUC * CFTC_RT
        CF_MSQ = CFI_MSQ  + RATIO_CFSUC * CFTC_MSQ
        CF_US  = 0.0      + RATIO_CFSUC * CFTC_US
        CF_SV  = 0.0      + RATIO_CFSUC * CFTC_SV

      ENDIF

      IF(L_CF_KAYS) THEN
        CALL CFTCKAYS( CFI , SV, CFTC, CFTC_CFI, CFTC_SV)
        ! WRITE (*,*) 'Kays Cf used'
!        Original:
        CF     = CFTC
        CF_HK  = CFTC_CFI * CFI_HK
        CF_RT  = CFTC_CFI * CFI_RT
        CF_MSQ = CFTC_CFI * CFI_MSQ
        CF_US  = 0.0
        CF_SV  = CFTC_SV

!        With RATIO_CFSUC to control how much the effect of suction is accounted for in CF correlation
!        RATIO_CFSUC = 1              Full Kays CF
!        RATIO_CFSUC = 0              Suction has no effect CF closure relation
        CF     = CFI     + RATIO_CFSUC * (CFTC - CFI)
        CF_HK  = CFI_HK  + RATIO_CFSUC * (CFTC_CFI - 1) * CFI_HK
        CF_RT  = CFI_RT  + RATIO_CFSUC * (CFTC_CFI - 1) * CFI_RT
        CF_MSQ = CFI_MSQ + RATIO_CFSUC * (CFTC_CFI - 1) * CFI_MSQ
        CF_US  = 0.0
        CF_SV  = 0.0    + RATIO_CFSUC * CFTC_SV
!        End of option

      ENDIF

!      WRITE (*,*) 'CF_SV = ' , CF_SV
!       IMPORTANT
!       Record that US_HK is not zero, so CF_HK_TOTAL = CF_HK + CF_US * US_HK
!      WRITE (*,*) 'CFI=',  CFI , '  CFTC=' ,CFTC , '  CF=', CF, '  SV=',SV
      RETURN
      END


      SUBROUTINE XIFSET2(IS)
      INCLUDE 'RFOIL.f90'
      INCLUDE 'RBL.f90'
      INCLUDE 'RBL_OPTIONS.f90'
!     Gael de Oliveira:
!     A slightly modified (by Gael) version of XIFSET subroutine (from Drela) for finding forced transition point in blow mode.
!      X2FORC ,XTR2(2)
!
      IF(XTR2(IS).GE.1.0) THEN
       X2FORC = XSSI(IBLTE(IS),IS)
       RETURN
      ENDIF
!
      IF(IS.EQ.1) THEN
!
!----- set approximate arc length value of forced transition point for SINVRT
       SFORCE = SLE + (    -SLE)*XTR2(IS)
!
!----- calculate actual arc length
       CALL SINVRT(SFORCE,XTR2(IS),X,XP,S,N)
!
!----- set BL coordinate value
       X2FORC = AMIN1( (SST - SFORCE) , XSSI(IBLTE(IS),IS) )
      ELSE
!
       SFORCE = SLE + (S(N)-SLE)*XTR2(IS)
       CALL SINVRT(SFORCE,XTR2(IS),X,XP,S,N)
       X2FORC = AMIN1( (SFORCE - SST) , XSSI(IBLTE(IS),IS) )
!
      ENDIF
!
      RETURN
      END

      SUBROUTINE SVBLOWSET_OLD(IBL , IS)
      INCLUDE 'RFOIL.f90'
      INCLUDE 'RBL.f90'
      INCLUDE 'RBL_OPTIONS.f90'

!     Gael de Oliveira : Module to compute equivalent blowing to a trip wire
! Onl
      ! First we need to find what the Cf at this point, the obstacle point, would be if the obstacle wasn'ty there.
      ! Given the BL expressions are parabolic, it is fair to extrapolate from upstream in the domain.
      ! We extrapolate linearly from the two previous points

! User input
!     D_TRIP = trip wire diameter / airfoil chord
!      D_TRIP = 0.00035


! Extract panel positions of earlier points
!     Reference:
!		A =  two panels upstream of obstacle
!		B =  one panel upstream of behind obstacle
!		C =  position of obstacle
!		D =  one panel after obstacle
      IBL_A = IBL - 2
      IBL_B = IBL - 1
      IBL_D = IBL + 1
      IBL_E = IBL + 2

! Arc distances
      XSI_A = XSSI(IBL_A,IS)
      XSI_B = XSSI(IBL_B,IS)
      XSI_C = XSSI(IBL,IS)
      XSI_D = XSSI(IBL_D,IS)
      XSI_E = XSSI(IBL_E,IS)

! Panel positions
      IP_B  = IPAN(IBL_B,IS)
      IP_C  = IPAN(IBL  ,IS)
      IP_D  = IPAN(IBL_D,IS)

      UEI_A = UEDG(IBL_A,IS)
      UEI_B = UEDG(IBL_B,IS)
      UEI_C = UEDG(IBL,IS)
! Extract Cf's of earlier points
!     TAU(IBL,IS)  = 0.5*R2*U2*U2*CF2
!     R2 = 1 in Xfoil when M=0. So we restrict to incompressible flow for simplicity
      TAU_A = TAU(IBL_A,IS)
      TAU_B = TAU(IBL_B,IS)
      TAU_C = TAU(IBL,IS)

      CF_A = 2 * TAU_A / UEI_A**2
      CF_B = 2 * TAU_B / UEI_B**2
      CF_C = 2 * TAU_C / UEI_C**2

! Obtain size of blowing interval
      ! W_BLOW = (XSI_D - XSI_B) / 2
      ! The interval is divided by two to account for triangular shape of suction distribution once discretized
      W_BLOW = (XSI_E - XSI_D) /2 + (XSI_D - XSI_C) + ( XSI_C - XSI_B) / 2

! Now extrapolate
      CF_C = CF_A + (XSI_C - XSI_A) * (CF_B - CF_A) / (XSI_B - XSI_A)
      CF_C = 2 * TAU_C / UEI_C**2
!      WRITE (*,*) CF_A , CF_B , CF_C , REINF
!      Extrapolation works beautifully on du97w300!

      SV_BLOW = 0.25 * CF_C * REINF * D_TRIP**2 / W_BLOW
!      WRITE (*,*) '  SV_BLOW = ' , SV_BLOW
!      SV_BLOW = D_TRIP / W_BLOW
!     Now decide if this is the point were we blow
!       By default we don't blow
      TR_BLOW = .FALSE.

!     Determine Arc Position of Blowing Point
      CALL XIFSET2(IS)
!      WRITE (*,*) '  SV_BLOW = ' , SV_BLOW , '  IBL = ' , IBL , ' X2FORC = ' , X2FORC

      IF(X2FORC .GT. XSI_C) THEN
        IF(X2FORC .LT. XSI_D) THEN
          TR_BLOW = .TRUE.

          ! Store suction speed for second panel
          SV_BLOWB = SV_BLOW
!          WRITE (*,*) '  SV_BLOW = ' , SV_BLOW , '  IBL = ' , IBL
          WRITE (*,*) '  SV_BLOW = ' , SV_BLOW , '  IBL = ' , IBL , ' X2FORC = ' , X2FORC
        ENDIF
      ENDIF

      IF(X2FORC .GT. XSI_B) THEN
        IF(X2FORC .LT. XSI_C) THEN
          TR_BLOW = .TRUE.
          ! Set suction speed of second panel same as on first panel
          SV_BLOW = SV_BLOWB
!          WRITE (*,*) '  SV_BLOW = ' , SV_BLOW , '  IBL = ' , IBL
          WRITE (*,*) '  SV_BLOW = ' , SV_BLOW , '  IBL = ' , IBL , ' X2FORC = ' , X2FORC
        ENDIF
      ENDIF

      ! Limit SVBLOW to a reasonable value
      SV_BLOW = MIN(SV_BLOW, BMAX)

      RETURN
      END


      SUBROUTINE SVBLOWSET(IBL , IS)
      INCLUDE 'RFOIL.f90'
      INCLUDE 'RBL.f90'
      INCLUDE 'RBL_OPTIONS.f90'

!     Gael de Oliveira : Module to compute equivalent blowing to a trip wire
! Onl
      ! First we need to find what the Cf at this point, the obstacle point, would be if the obstacle wasn'ty there.
      ! Given the BL expressions are parabolic, it is fair to extrapolate from upstream in the domain.
      ! We extrapolate linearly from the two previous points

! User input
!     D_TRIP = trip wire diameter / airfoil chord
!      D_TRIP = 0.00035


! Extract panel positions of earlier points
!     Reference:
!		A =  two panels upstream of obstacle
!		B =  one panel upstream of behind obstacle
!		C =  position of obstacle
!		D =  one panel after obstacle
      IBL_A = IBL - 2
      IBL_B = IBL - 1
      IBL_D = IBL + 1
      IBL_E = IBL + 2

! Arc distances
      XSI_A = XSSI(IBL_A,IS)
      XSI_B = XSSI(IBL_B,IS)
      XSI_C = XSSI(IBL,IS)
      XSI_D = XSSI(IBL_D,IS)
      XSI_E = XSSI(IBL_E,IS)

! Panel positions
      IP_B  = IPAN(IBL_B,IS)
      IP_C  = IPAN(IBL  ,IS)
      IP_D  = IPAN(IBL_D,IS)

      UEI_C = UEDG(IBL,IS)
! Extract Cf's of earlier points
!     TAU(IBL,IS)  = 0.5*R2*U2*U2*CF2
!     R2 = 1 in Xfoil when M=0. So we restrict to incompressible flow for simplicity

      TAU_C = TAU(IBL,IS)

! If we are in exploratory mode, change nothing
      IF (TR_BLOW_EXPLORE) THEN
        CF_C = 0
      ENDIF
! If we are setting the blowing, set the Cf to the value obtained in exploratory mode
      IF (.NOT. TR_BLOW_EXPLORE) THEN
        ! WRITE (*,*) '  CF_REF = ' , CF_REF
        CF_C = CF_REF(IS)
      ENDIF


! Obtain size of blowing interval
      ! W_BLOW = (XSI_D - XSI_B) / 2
      ! The interval is divided by two to account for triangular shape of suction distribution once discretized
      W_BLOW = (XSI_E - XSI_D) /2 + (XSI_D - XSI_C) + ( XSI_C - XSI_B) / 2

!      WRITE (*,*) CF_A , CF_B , CF_C , REINF
!      Extrapolation works beautifully on du97w300!

      SV_BLOW = 0.25 * CF_C * REINF * D_TRIP**2 / W_BLOW
!      WRITE (*,*) '  SV_BLOW = ' , SV_BLOW
!      SV_BLOW = D_TRIP / W_BLOW
!     Now decide if this is the point were we blow
!       By default we don't blow
      TR_BLOW = .FALSE.

!     Determine Arc Position of Blowing Point
      CALL XIFSET2(IS)
!      WRITE (*,*) '  SV_BLOW = ' , SV_BLOW , '  IBL = ' , IBL , ' X2FORC = ' , X2FORC

      IF(X2FORC .GT. XSI_C) THEN
        IF(X2FORC .LT. XSI_D) THEN
          TR_BLOW = .TRUE.
          ! Store suction speed for second panel
          SV_BLOWB = SV_BLOW
!          WRITE (*,*) '  SV_BLOW = ' , SV_BLOW , '  IBL = ' , IBL
          ! WRITE (*,*) '  SV_BLOW = ' , SV_BLOW , '  IBL = ' , IBL , ' X2FORC = ' , X2FORC

          ! Store Reference Cf
          IF (TR_BLOW_EXPLORE) THEN
            CF_REF(IS) = 2 * TAU_C / UEI_C**2
          ENDIF

		  ! Store Reference Blowing speed for reporting
          SV_BLOW_REPORTING(IS) = SV_BLOW
          ITR_BLOW_REPORTING(IS) = IBL
        ENDIF
      ENDIF

      IF(X2FORC .GT. XSI_B) THEN
        IF(X2FORC .LT. XSI_C) THEN
          TR_BLOW = .TRUE.
          ! Set suction speed of second panel same as on first panel
          SV_BLOW = SV_BLOWB
!          WRITE (*,*) '  SV_BLOW = ' , SV_BLOW , '  IBL = ' , IBL
!          WRITE (*,*) '  SV_BLOW = ' , SV_BLOW , '  IBL = ' , IBL , ' X2FORC = ' , X2FORC
        ENDIF
      ENDIF

      ! Limit SVBLOW to a reasonable value
      SV_BLOW = MIN(SV_BLOW, BMAX)

!      IF (TR_BLOW) THEN
!        WRITE (*,*) '  SV_BLOW = ' , SV_BLOW , '  IBL = ' , IBL , ' X2FORC = ' , X2FORC , '  TR_BLOW_EXPLORE = ' , TR_BLOW_EXPLORE
!        ! WRITE (*,*) '  CF_REF = ' , CF_REF
!      ENDIF

      IF (TR_BLOW_EXPLORE) THEN
        SV_BLOW = 0
        TR_BLOW = .FALSE.
      ENDIF




      RETURN
      END

      SUBROUTINE XTR_OVERRIDING_CONTROL
      ! A subroutine to define whether the classical forced transition is overriden to give good results with trip-wire-equivalent blowing
      ! Gael de Oliveira
      INCLUDE 'RFOIL.f90'
!      INCLUDE 'RBL.f90'
      INCLUDE 'RBL_OPTIONS.f90'

      ! By default take classical forced transition
      XTR1(1) = XTR1BAK(1)
      XTR1(2) = XTR1BAK(2)

      ! Only if blowing transition is active, and overriding is enabled, force transition to XTR_blowing - DXTR
      IF ((XTR2(1) .LT. 1.0) .OR. (XTR2(2) .LT. 1.0)) THEN
        IF (L_TR_OVERRIDING) THEN
          XTR1(1) = XTR2(1) - DXTR
          XTR1(2) = XTR2(2) - DXTR

          ! Fudge to a minimum transition point of x/c = 0.005 to avoid blowing up
          XTR1(1) = MAX(XTR1(1) , 0.005)
          XTR1(2) = MAX(XTR1(2) , 0.005)
        ENDIF
      ENDIF

      RETURN
      END
	  
	  
! plasmacode -->	  
      SUBROUTINE WXP_FUNCTION(X , LPP, XPP , WXP, WXP_X)
      IMPLICIT REAL (A-Z)
!      WXP is the Plasma Force Field Weighthing Function for the X (body surface) direction 
!      TEMPORARY: Ignore difference between (X)S and X for now! Check line 5052 to see 
!                 how to sort this out
!      Gael de Oliveira

!     Define PI Locally, to avoid namespace fizzling risks with includes
      PIVAL = 3.1415926
!     Set returns at 0 by default and change them if we are in the actuation region	  
      WXP   = 0.0
      WXP_X = 0.0
!     Make numerator and denominator independently to avoid instability 
!     FPexception with 0 length (LPP)
      SINUMER = X - XPP
      SIDENOM = LPP
      IF (SIDENOM .GT. 1E-7) THEN
        IF ((SINUMER .GT. 0.0) .AND. (SINUMER .LT. LPP)) THEN
            WXP   = 0.5 * PIVAL * SIN(PIVAL * SINUMER / SIDENOM)
            WXP_X = 0.5 * (PIVAL**2/SIDENOM) &
                        * COS(PIVAL * SINUMER / SIDENOM)
        ENDIF
      ENDIF
!     Done, we can return

      RETURN
      END
! <-- plasmacode

! plasmacode -->	
      SUBROUTINE XIFSETdemo(IS)
      INCLUDE 'RFOIL.f90'
      INCLUDE 'RBL.f90'
!
      IF(XTR1(IS).GE.1.0) THEN
       XIFORC = XSSI(IBLTE(IS),IS)
       RETURN
      ENDIF
!
      IF(IS.EQ.1) THEN
!
!----- set approximate arc length value of forced transition point for SINVRT
       SFORCE = SLE + (    -SLE)*XTR1(IS)
!
!----- calculate actual arc length
       CALL SINVRT(SFORCE,XTR1(IS),X,XP,S,N)
!
!----- set BL coordinate value
       XIFORC = AMIN1( (SST - SFORCE) , XSSI(IBLTE(IS),IS) )
      ELSE
!
       SFORCE = SLE + (S(N)-SLE)*XTR1(IS)
       CALL SINVRT(SFORCE,XTR1(IS),X,XP,S,N)
       XIFORC = AMIN1( (SFORCE - SST) , XSSI(IBLTE(IS),IS) )
!
      ENDIF
!
      RETURN
      END
	  
! <-- plasmacode

! plasmacode --
      SUBROUTINE XIFSETPLASMArfoilold(IS)
      INCLUDE 'RFOIL.f90'
      INCLUDE 'RBL.f90'
      INCLUDE 'RBL_OPTIONS.f90'
!     Gael de Oliveira:
!     A slightly modified (by Gael) version of XIFSET subroutine (from Drela) for finding leading 
!     edge of plasma actuator (also inspired on XIFSET2)
!     Consider plasma on both side for now [adjust later]
!        INPUT:  XPP1 (actual position in fixed frame)
!       OUTPUT: XSPP1 (current position in fluctuating frame)
!

! Take it out for wake
!      IF(XPP1.GE.1.0) THEN
!!      Wake
!       XSPP = XSSI(IBLTE(IS),IS)
!       RETURN
!      ENDIF
!
      IF(IS.EQ.1) THEN
!      Suction Side

!----- set approximate arc length value of forced transition point for SINVRT
       SPP = SLE + (    -SLE)*XPP
!
!----- calculate actual arc length
       CALL SINVRT(SPP,XPP,X,XP,S,N)
!
!----- set BL coordinate value
       XSPP = AMIN1( (SST - SPP) , XSSI(IBLTE(IS),IS) )
      ELSE
!		  Pressure Side
!
       SPP = SLE + (S(N)-SLE)*XPP
       CALL SINVRT(SPP,XPP,X,XP,S,N)
       XSPP = AMIN1( (SPP - SST) , XSSI(IBLTE(IS),IS) )
!
      ENDIF
!
      RETURN
      END
! <-- plasmacode

! plasmacode -->	  
      SUBROUTINE XIFSETPLASMAxfoil(IS)
!-----------------------------------------------------
!     Sets forced-transition BL coordinate locations.
!-----------------------------------------------------
!     Include needed stuff, Xfoil/Rfoil
!      INCLUDE 'XFOIL.INC'
!      INCLUDE 'XBL.INC'
      INCLUDE 'RFOIL.f90'
      INCLUDE 'RBL.f90'
	  INCLUDE 'PLASMA_OPTIONS.f90'
!
      IF(XPP.GE.1.0) THEN
       XSPP = XSSI(IBLTE(IS),IS)
       RETURN
      ENDIF
!
      CHX = XTE - XLE
      CHY = YTE - YLE
      CHSQ = CHX**2 + CHY**2
!
!---- calculate chord-based x/c, y/c
      DO 10 I=1, N
        W1(I) = ((X(I)-XLE)*CHX + (Y(I)-YLE)*CHY) / CHSQ
        W2(I) = ((Y(I)-YLE)*CHX - (X(I)-XLE)*CHY) / CHSQ
 10   CONTINUE
!
      CALL SPLIND(W1,W3,S,N,-999.0,-999.0)
      CALL SPLIND(W2,W4,S,N,-999.0,-999.0)
!
      IF(IS.EQ.1) THEN
!
!----- set approximate arc length (SPP) of plasma LE point for SINVRT
       SPP = SLE + (S(1)-SLE)*XPP
!
!----- calculate actual arc length (SPP)
       CALL SINVRT(SPP,XPP,W1,W3,S,N)
!
!----- set BL coordinate value of actuator LE position (XSPP)
       XSPP = MIN( (SST - SPP) , XSSI(IBLTE(IS),IS) )
!
!      WRITE(*,1010) SPP
! 1010  FORMAT(/' SPP side 1', F8.4,'  ***')

! 6.b Compute CPHI with care to avoid division by zero (this could be moved out of the loop!)
       IF (LPP .GT. 0.000001) THEN
!        Business as usual if LPP is decently large (ridiculously small!)
!        From report:
!         CPHI = FTP / LPP
!        But something does not make sense here!
!         CPHI = TPP * FTP / LPP
!		 The TPP cancels, and maybe so does the LPP. Maybe this is better, 
!		 try for now, and verify later (lets remember to check this and mention it! 
!		 so that a secure approach can be followed!)
!         CPHI = - FTP
!        Revert to fundamented choice (from report):
         CPHI = FTP / LPP
        ELSE
!       Fudge to zero is LPP is undefined (must be non-zero to avoid division by zero)
         CPHI = 0.0
       ENDIF
	   
       WRITE(*,1011) CPHI
 1011  FORMAT(/'Side 1 , CPHI Value = ', F8.6)


      ELSE
!----- same for bottom side (over)
!
       SPP = SLE + (S(N)-SLE)*XPP
       CALL SINVRT(SPP,XPP,W1,W3,S,N)
       XSPP = MIN( (SPP - SST) , XSSI(IBLTE(IS),IS) )
	   
!       WRITE(*,1011) SPP
! 1011  FORMAT(/' SPP side 2', F8.4,'  ***')
! Put no plasma on side 2
        CPHI = 0.0
        WRITE(*,1012) CPHI
 1012   FORMAT(/'Side 2 , CPHI Value = ', F8.6)
      ENDIF
!

      IF(XSPP .LT. 0.0) THEN
       WRITE(*,1000) IS
 1000  FORMAT(/' ***  Stagnation point is past trip on side',I2,'  ***')
       XIFORC = XSSI(IBLTE(IS),IS)
      ENDIF
!

      RETURN
      END
! <-- plasmacode

! plasmacode --
      SUBROUTINE XIFSETPLASMA(IS)
      INCLUDE 'RFOIL.f90'
      INCLUDE 'RBL.f90'
      INCLUDE 'PLASMA_OPTIONS.f90'
!     Gael de Oliveira:
!     A slightly modified (by Gael) version of XIFSET subroutine (from 
!     Drela) for finding leading edge of plasma actuator (also inspired 
!     on XIFSET2) Consider plasma on both side for now [adjust later]
!        INPUT:  XPP1 (actual position in fixed frame)
!       OUTPUT: XSPP1 (current position in fluctuating frame)
!

! Take it out for wake
!      IF(XPP1.GE.1.0) THEN
!!      Wake
!       XSPP = XSSI(IBLTE(IS),IS)
!       RETURN
!      ENDIF
!
      IF(IS.EQ.1) THEN
!      Suction Side

!----- set approximate arc length value of forced transition point for SINVRT
       SPP = SLE + (    -SLE)*XPP
!
!----- calculate actual arc length
       CALL SINVRT(SPP,XPP,X,XP,S,N)
!
!----- set BL coordinate value
       XSPP = AMIN1( (SST - SPP) , XSSI(IBLTE(IS),IS) )
      ELSE
!		  Pressure Side
!
       SPP = SLE + (S(N)-SLE)*XPP
       CALL SINVRT(SPP,XPP,X,XP,S,N)
       XSPP = AMIN1( (SPP - SST) , XSSI(IBLTE(IS),IS) )
!
      ENDIF
	  
! 6.b Compute CPHI with care to avoid division by zero (this could be moved out of the loop!)
      IF (LPP .GT. 0.000001) THEN
!        Business as usual if LPP is decently large (ridiculously small!)
!        From report:
!         CPHI = FTP / LPP (The TPP cancels)
         CPHI = FTP / LPP
      ELSE
!        Fudge to zero is LPP is undefined (must be non-zero to avoid
!        division by zero):
         CPHI = 0.0
      ENDIF
	   
      IF(IS.NE.ISPP) THEN
!       Plasma only acts on upper side for now
        CPHI = 0.0
      ENDIF
! Print out some diagnostics
!      WRITE(*,1012) IS , XSPP, CPHI
! 1012   FORMAT(/'XIFSETPLASMA: Side ', I2 , ' ,   XSPP=' , F8.6, ' ,   CPHI=', F8.6)
	  
!
      RETURN
      END
! <-- plasmacode



! plasmacode -->
      SUBROUTINE CEITpoly(HK , RT, TPPT, CEI, CEI_HK, CEI_RT, CEI_TPPT)
      IMPLICIT REAL (A-Z)
!      Polynomial fit Version - possibly replaced with trilinear interpolation 
!      closure (slower, but likely better accuracy and hopefully stability!)
!      CEI is the Force Interaction Coefficient, as it was defined in the report
!      HK  is the shape factor, as used for CF (I ignore compressibility for H-HK choice)
!          (anyway, compressibility is usually associated with guns! why would you want to 
!           do such a counterproductive device as a gun?)
!      RT  is the Reynold Theta (momentum thickness)
!      CEI_HK, CEI_RT anf CEI_TPPT derivatives follow usual notation/convention
!      Gael de Oliveira

!      For development purposes only, use a constant representative value of CEI
!	   CEI = 0.1699, obtained by computing the mean of the matlab database
!                    mean(mean(mean(CEI.cei_grid))))
!      CEI = 0.6   , HK=1.6, RT=2000, TPPT=5
      CEI      = 0.6
      CEI_HK   = 0.0
      CEI_RT   = 0.0
      CEI_TPPT = 0.0
!     Get CEI at RT 500 and RT 10000 and interpolate between the two polynomial fits made at these RT	  
      CALL CEIT500(HK , TPPT, CEI5C , CEI5C_HK , CEI5C_TPPT)
      CALL CEIT10000(HK , TPPT, CEI10K, CEI10K_HK, CEI10K_TPPT)
!     Mingle
      LAMBDA = (RT - 500.0) / (10000.0 - 500.0)
      CEI      = CEI5C      * (1.0-LAMBDA) + CEI10K    * LAMBDA
      CEI_HK   = CEI5C_HK   * (1.0-LAMBDA) + CEI10K_HK * LAMBDA
      CEI_TPPT = CEI5C_TPPT * (1.0-LAMBDA) + CEI10K_TPPT*LAMBDA
!	  Get mean slope
!      CEI_RT   = (CEI10K - CEI5C) / (10000 - 500)	  
!     Done, we can return
      IF (RT .GT. 10000) THEN
        WRITE(*,1988) RT
      ENDIF
 1988 FORMAT('CEIT:WARNING: RT =', F12.6)
      RETURN
      END
! <-- plasmacode 
	  
! plasmacode -->
      SUBROUTINE CEIT500(HK , TPPT, CEI, CEI_HK, CEI_TPPT)
      IMPLICIT REAL (A-Z)
!      CEI is the Force Interaction Coefficient, as it was defined in the report
!      HK  is the shape factor, as used for CF (I ignore compressibility for H-HK choice)
!          (anyway, compressibility is usually associated with guns! why would you want to 
!           do such a counterproductive device as a gun?)
!      RT  is the Reynold Theta (momentum thickness)
!      CEI_HK, CEI_RT anf CEI_TPPT derivatives follow usual notation/convention
!      Gael de Oliveira

!      For development purposes only, use a constant representative value of CEI
!	   CEI = 0.1699, obtained by computing the mean of the matlab database
!                    mean(mean(mean(CEI.cei_grid))))
!      CEI = 0.6   , HK=1.6, RT=2000, TPPT=5

! For RT=500
      p00 =      0.9983
      p10 =     -0.8585
      p01 =      0.1558
      p20 =      0.2544
      p11 =    -0.02162
      p02 =    -0.01172
      p30 =     -0.0318
      p21 =   -0.004521
      p12 =    0.004589
      p40 =    0.001437
      p31 =   0.0005442
      p22 =  -0.0003136
	  
      x  = HK
      y  = TPPT

      CEI =    p00   + p10*x        + p01*y         &
        + p20*x**2   + p11*x*y      + p02*y**2      &
        + p30*x**3   + p21*x**2*y   + p12*x*y**2    &
        + p40*x**4   + p31*x**3*y   + p22*x**2*y**2
	 
      CEI_HK =       + p10*x                        &
        + 2*p20*x    + p11  *y                      &
        + 3*p30*x**2 + 2*p21*x*y    + p12  *y**2    &
        + 4*p40*x**3 + 3*p31*x**2*y + 2*p22*x*y**2
	 
      CEI_TPPT =                        p01         &
                     + p11*x        + 2*p02*y       &
                     + p21*x**2     + 2*p12*x*y     &
                     + p31*x**3     + 2*p22*x**2*y
      !TEMPORARY
!     Done, we can return

      RETURN
      END
! <-- plasmacode 
	  
! plasmacode -->
      SUBROUTINE CEIT10000(HK , TPPT, CEI, CEI_HK, CEI_TPPT)
      IMPLICIT REAL (A-Z)
!      CEI is the Force Interaction Coefficient, as it was defined in the report
!      HK  is the shape factor, as used for CF (I ignore compressibility for H-HK choice)
!          (anyway, compressibility is usually associated with guns! why would you want to 
!           do such a counterproductive device as a gun?)
!      RT  is the Reynold Theta (momentum thickness)
!      CEI_HK, CEI_RT anf CEI_TPPT derivatives follow usual notation/convention
!      Gael de Oliveira

!      For development purposes only, use a constant representative value of CEI
!	   CEI = 0.1699, obtained by computing the mean of the matlab database
!                    mean(mean(mean(CEI.cei_grid))))
!      CEI = 0.6   , HK=1.6, RT=2000, TPPT=5

! For RT=10000
      p00 =       1.407
      p10 =     -0.8978
      p01 =   -0.007943
      p20 =      0.1868
      p11 =     0.05688
      p02 =    -0.00152
      p30 =    -0.01613
      p21 =    -0.01374
      p12 =  0.00004121
      p40 =   0.0004755
      p31 =    0.000802
      p22 =  0.00009712
	  
      x  = HK
      y  = TPPT

      CEI =    p00   + p10*x        + p01*y         &
        + p20*x**2   + p11*x*y      + p02*y**2      &
        + p30*x**3   + p21*x**2*y   + p12*x*y**2    &
        + p40*x**4   + p31*x**3*y   + p22*x**2*y**2
	 
      CEI_HK =       + p10*x                        &
        + 2*p20*x    + p11  *y                      &
        + 3*p30*x**2 + 2*p21*x*y    + p12  *y**2    &
        + 4*p40*x**3 + 3*p31*x**2*y + 2*p22*x*y**2
	 
      CEI_TPPT =                        p01         &
                     + p11*x        + 2*p02*y       &
                     + p21*x**2     + 2*p12*x*y     &
                     + p31*x**3     + 2*p22*x**2*y
      !TEMPORARY
!     Done, we can return

      RETURN
      END
! <-- plasmacode 

! --> plasmacode
      SUBROUTINE CEIT(HK , RT, TPPT, CEI, CEI_HK, CEI_RT, CEI_TPPT)
! PROGRAM CEI_EXAMPLE

! TRILINEAR INTERPOLATION OF THE ENERGY INTERACTION COEFFICIENT
! THIS FUNCTION VERSION IS VALIDATED:
!     + for the search algorithm
!     + for the trilinear interpolation
!     + for the hk, rt and tppt derivatives
!
! Gael de Oliveira, Ricardo Pereira, Feb. 2015
! Development Stage / All rights reserved
!
! Acknowledgement: based on the trilinear formula from wikipedia
!           http://en.wikipedia.org/wiki/Trilinear_interpolation
!
! Coding Convention:
!    We follow implicit type rules REAL(A-H,O-Z), INTEGER(I-N)
	
! Declare out-of-bounds variable
     LOGICAL L_OUT_OF_BOUNDS
! datacode -->
! Also add access to verbose mode flags (do not include whole RBL_OPTIONS.f90 file
! to avoid crossnaming errors and facilitate compiler optimizations!)
     REAL EPS_VPAR
     LOGICAL L_PRECISION
     LOGICAL L_VERBOSE
     COMMON /EPSCOMMON/EPS_VPAR , L_PRECISION, L_VERBOSE
! <-- datacode

! Start by including the data (file generated from matlab
! with custom script)
     INCLUDE 'CEI_DATA.f90'
! Typical variable declarations found in the include (together with data)
!    INTEGER N_HK, N_RT, N_TPPT
!    REAL, DIMENSION(1,11) :: HK_RANGE
!    REAL, DIMENSION(1,11) :: RT_RANGE
!    REAL, DIMENSION(1,11) :: TPPT_RANGE
!    REAL, DIMENSION(11,11,11) :: CEI_TABLE

! Complete Include with explicit definition of first tppt point (0)
! To allow for compilation without finit-zero, no runtime overhead involved!
     TPPT_RANGE(1,1) = 0.0
     
! Define CEI evaluation point (for PROGRAM only, in subroutine, it comes as argument!)
!      HK   = 6.5
!      RT   = 1000
!      TPPT = 0.6
      
! Set Search Variables to allow for easy fudging without compromising data integrity
! (Xfoil's Fortran programming style does not specify intent!)
! This does not allow warning but is probaly more efficient than if constructs, and
! much more compact to write! Assumes range vectors are growing monotonically!
      HK_S   = MAX(HK     , HK_RANGE(1,1)       )
      HK_S   = MIN(HK_S   , HK_RANGE(1,N_HK)    )
      
      RT_S   = MAX(RT     , RT_RANGE(1,1)       )
      RT_S   = MIN(RT_S   , RT_RANGE(1,N_RT)    )
      
      TPPT_S = MAX(TPPT   , TPPT_RANGE(1,1)     )
      TPPT_S = MIN(TPPT_S , TPPT_RANGE(1,N_TPPT))

! Diagnostic Output 
!      WRITE(*,*) 'REQUESTED POINT:'
!      WRITE(*,1100) HK  , RT  , TPPT
!      WRITE(*,*) 'SEARCHED POINT:'
!      WRITE(*,1100) HK_S, RT_S, TPPT_S 
!      WRITE(*,*) 'START SEARCH:'

! SEARCH ALGORITHM
! Look for upper nearest neighbour of desired CEI point in (HK,RT,TPPT) space
! This code is a bit a hard to understand (deep nesting) but fairly efficient for a
! simple algorithm (in terms of theorethical computational complexity, might still
! pose efficiency issues on architectures with inneficient execution control, but
! I seem to remember that x64 is fairly good at that!)
      DO 100 I_HK=2,N_HK 
        ! Check if we are in proper HK interval
        IF (HK_S .LE. HK_RANGE(1,I_HK)) THEN
          ! Loop through RT only if that is the case (and not otherwise)
          DO 200 I_RT=2,N_RT
            ! Check if we are in proper RT interval
            IF (RT_S .LE. RT_RANGE(1,I_RT)) THEN
              ! Loop through TPPT only if that is the case (and not otherwise)
              DO 300 I_TPPT=2,N_TPPT
                ! Check if we are in proper TPPT interval
                IF (TPPT_S .LE.  TPPT_RANGE(1,I_TPPT)) THEN
                  ! Write some output
!                  WRITE(*,1000) I_HK, I_RT, I_TPPT
!                  WRITE(*,1100) HK_RANGE(1,I_HK) , RT_RANGE(1,I_RT), TPPT_RANGE(1,I_TPPT)
                  GOTO 500
                ENDIF
 300          CONTINUE
            ENDIF
 200      CONTINUE
        ENDIF
 100  CONTINUE
 
! CYCLE BREAK LABEL (place at first computation when commenting festivities!) 
! 500  WRITE(*,*) 'FOUND UPPER NEIGHBOUR!!!'
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Search complete, go on with trilinear Interpolation
 ! Store lower bound (_0 subscript at I-1 index)
 500  HK_0   = HK_RANGE(1,I_HK-1)
      RT_0   = RT_RANGE(1,I_RT-1)
      TPPT_0 = TPPT_RANGE(1,I_TPPT-1)
 ! Store upper bound (_1 subscript at I   index)
      HK_1   = HK_RANGE(1,I_HK)
      RT_1   = RT_RANGE(1,I_RT)
      TPPT_1 = TPPT_RANGE(1,I_TPPT)
 ! Compute Scaled Differences (in independent variable space)
      HK_D   = (  HK_S - HK_0   ) / (  HK_1 - HK_0  )
      RT_D   = (  RT_S - RT_0   ) / (  RT_1 - RT_0  )
      TPPT_D = (TPPT_S - TPPT_0 ) / (TPPT_1 - TPPT_0)
 ! Store 8 surrounding value points
 ! (Notice beautiful mnemonic with right->left binary counting <3)
      V000   = CEI_TABLE(I_HK-1 , I_RT-1 , I_TPPT-1)
      V100   = CEI_TABLE(I_HK   , I_RT-1 , I_TPPT-1)
      V010   = CEI_TABLE(I_HK-1 , I_RT   , I_TPPT-1)
      V110   = CEI_TABLE(I_HK   , I_RT   , I_TPPT-1)
      V001   = CEI_TABLE(I_HK-1 , I_RT-1 , I_TPPT  )
      V101   = CEI_TABLE(I_HK   , I_RT-1 , I_TPPT  )
      V011   = CEI_TABLE(I_HK-1 , I_RT   , I_TPPT  )
      V111   = CEI_TABLE(I_HK   , I_RT   , I_TPPT  )
 ! Evaluate two bit coefficient set (fix plane in HK dimension)
      C00    = V000 * (1-HK_D)  +  V100 * HK_D
      C10    = V010 * (1-HK_D)  +  V110 * HK_D
      C01    = V001 * (1-HK_D)  +  V101 * HK_D
      C11    = V011 * (1-HK_D)  +  V111 * HK_D
 ! Make single bit coefficient set (fix line in HK-RT dimension)
      C0     = C00  * (1-RT_D)  +   C10 * RT_D
      C1     = C01  * (1-RT_D)  +   C11 * RT_D
! Make single point! (fix point in HK-RT-TPPT dimension)
      C      = C0   * (1-TPPT_D)+    C1 * TPPT_D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Process into CEI and derivatives
      CEI      = C
! TPPT derivative is XXeasyXX (that's why I fucked up on it, right?)   
! Two step chain rule to avoid errors!
!     First derivative of C to TPPT_D (scaled diff, not same as TPPT)
      C_TPPT_D = C1 - C0
!     Derivative of TPPT_D scaled difference to TPPT_S
      TPPT_D_TPPTS = (  1        ) / (TPPT_1 - TPPT_0)
!     Combination via chain rule
      C_TPPTS  = C_TPPT_D * TPPT_D_TPPTS
      CEI_TPPT = C_TPPTS
! RT derivative by chain rule 
!     via C0 and C1 derivatives
      C_C0      =        (1-TPPT_D)
      C_C1      =                            TPPT_D
!     and C0/1_RT_D derivatives
      C0_RT_D   = C10 - C00
      C1_RT_D   = C11 - C01
!     and RT_D_RTS derivatives
      RT_D_RTS  = (  1             ) / (  RT_1 - RT_0  )
!     before combination with chain rule
      C_RTS     = C_C0 * C0_RT_D * RT_D_RTS + &
                  C_C1 * C1_RT_D * RT_D_RTS
      CEI_RT    = C_RTS
! HK derivative is cumbersome, but also goes via chain rule
!     first on derivatives of single bit coefficient to double bit coefs
      C0_C00    =        (1-RT_D)
      C0_C10    =                            RT_D
      
      C1_C01    =        (1-RT_D)
      C1_C11    =                            RT_D
!     then on derivatives of double bit coefs to HK_D (scaled HK differences)
      C00_HK_D  = V100 - V000
      C10_HK_D  = V110 - V010
      C01_HK_D  = V101 - V001
      C11_HK_D  = V111 - V011
!     and finally (almost finally) on HK_D to HK(S)
      HK_D_HKS  = (  1             ) / (  HK_1 - HK_0  )
!     get explicit derivatives of double bit coefs to HK(S) (05 will get rid of
!     outcoming inneficiency, in principle!) 
      C00_HKS   = C00_HK_D * HK_D_HKS
      C10_HKS   = C10_HK_D * HK_D_HKS
      C01_HKS   = C01_HK_D * HK_D_HKS
      C11_HKS   = C11_HK_D * HK_D_HKS
!     and now do the same for single bit coefs (climbing the complexity ladder step
!     by step)
      C0_HKS    = C0_C00*C00_HKS + C0_C10*C10_HKS
      C1_HKS    = C1_C01*C01_HKS + C1_C11*C11_HKS
!     complete with C derivatives
      C_HKS     = C_C0*C0_HKS + C_C1*C1_HKS
      CEI_HK    = C_HKS
!     Done!

! Now check and deal with out of bounds events
      L_OUT_OF_BOUNDS = .FALSE.

! Now Fudge Derivatives to 0 for out of bounds conditions (and warn)
!     For HK first (Shape Factor)
      IF (ABS(HK-HK_S) .GT. 0.001) THEN
          CEI_HK = 0
          L_OUT_OF_BOUNDS = .TRUE.
      ENDIF
!     Then for RT (Reynolds Theta) 
      IF (ABS(RT-RT_S) .GT. 1) THEN
          CEI_RT = 0
          L_OUT_OF_BOUNDS = .TRUE.
      ENDIF
!     And finally for TPPTS
      IF (ABS(TPPT-TPPT_S) .GT. 0.001) THEN
          CEI_TPPT = 0
          L_OUT_OF_BOUNDS = .TRUE.
      ENDIF

! datacode -->
!     Write it out if we are out-of-bounds AND verbose is on:
      IF (L_VERBOSE) THEN
        IF (L_OUT_OF_BOUNDS) THEN
          WRITE(*,*   ) 'CEI_CLOSURE:WARNING:OUT_OF_BOUNDS'
          WRITE(*,*   ) 'Requested/Fudged:'
          WRITE(*,1100) HK, RT, TPPT
          WRITE(*,1100) HK_S, RT_S, TPPT_S
          WRITE(*,*   ) 'Result:'
          WRITE(*,1600) CEI , CEI_HK, CEI_RT, CEI_TPPT
        ENDIF
      ENDIF
! <-- datacode
      
! Display some diagnostics now
!      WRITE(*,*   ) 'TRILINEAR INTERPOLATION STEPS:'
!      WRITE(*,1000) I_HK, I_RT, I_TPPT
!      WRITE(*,1100) HK_0, RT_0, TPPT_0
!      WRITE(*,1100) HK_1, RT_1, TPPT_1
!      WRITE(*,1200) HK_D, RT_D, TPPT_D
!      WRITE(*,1300) V000, V100, V010, V110
!      WRITE(*,1310) V001, V101, V011, V111
!      WRITE(*,1400) C00 , C10 , C01 , C11 
!      WRITE(*,*   ) 'RESULT:'
!      WRITE(*,1600) CEI , CEI_HK, CEI_RT, CEI_TPPT
 
! 1000 FORMAT('I_HK=', I12  , '   I_RT=', I12  , '   I_TPPT=', I12)
 1100 FORMAT('  HK=', F12.4, '     RT=', F12.1, '     TPPT=', F12.4)
! 1200 FORMAT(' dHK=', F12.4, '    dRT=', F12.4, '    dTPPT=', F12.4)
! 1300 FORMAT('V000=', F12.4, '   V100=', F12.4, '     V010=', F12.4 , '     V110=', F12.4)
! 1310 FORMAT('V001=', F12.4, '   V101=', F12.4, '     V011=', F12.4 , '     V111=', F12.4)
! 1400 FORMAT(' C00=', F12.4, '    C10=', F12.4, '      C01=', F12.4 , '      C11=', F12.4)
 1600 FORMAT(' CEI=', F12.4, ' CEI_HK=', F12.4, '   CEI_RT=', F12.8 , ' CEI_TPPT=', F12.4)

 !END PROGRAM
      RETURN
      END
! <-- plasmacode

! mlcode -->
      SUBROUTINE SD5_FROM_H(BD5, BD5_XD5, XD5_H, SD5, SD5_H)
      IMPLICIT REAL (A-Z) 
      ! Make dimensional shape function from bernstein polynomials
      SD5   = BD5
      SD5_H = BD5_XD5 * XD5_H
      
      RETURN
      END
	  
      SUBROUTINE XD5_FROM_H(H, THMIN, THMAX, XD5, XD5_H)
      IMPLICIT REAL (A-Z) 
      ! Map H intervention region into unit interval
      XD5   = (H - THMIN) / (THMAX - THMIN)
      XD5_H =  1.0        / (THMAX - THMIN)
      ! Fudge for values of H above intervention region
      IF(H .GT. THMAX) THEN
             XD5   = 1.0
             XD5_H = 0.0
      ENDIF
      ! Fudge for values of H below intervention region
      IF(H .LT. THMIN) THEN
             XD5   = 0.0
             XD5_H = 0.0
      ENDIF
      RETURN
      END
	  
      SUBROUTINE BST_D5( X, A1, A2, A3, A4, A5, A6, BD5, BD5_X)
      IMPLICIT REAL (A-Z) 
      
      ! Make Bernstein polynomial basis for degree 5 (order 6)
      BD5_R0 =        X**5
      BD5_R1 =  5.0 * X**4 * (1.0-X)
      BD5_R2 = 10.0 * X**3 * (1.0-X)**2
      BD5_R3 = 10.0 * X**2 * (1.0-X)**3
      BD5_R4 =  5.0 * X    * (1.0-X)**4
      BD5_R5 =               (1.0-X)**5
      ! Combine basis linearly
      BD5    = A1*BD5_R0 + A2*BD5_R1 + A3*BD5_R2 + A4*BD5_R3 + A5*BD5_R4 + A6*BD5_R5
      
      ! Now make derivatives of polynomials from Bernstein basis of degree 5 (order 6)
      BD5_R0_X =  1.0*(5.0 *X**4                                     )
      BD5_R1_X =  5.0*(4.0 *X**3 *(1.0-X)    -  1.0 *X**4            )
      BD5_R2_X = 10.0*(3.0 *X**2 *(1.0-X)**2 -  2.0 *X**3 *(1.0-X)   )
      BD5_R3_X = 10.0*(2.0 *X    *(1.0-X)**3 -  3.0 *X**2 *(1.0-X)**2)
      BD5_R4_X =  5.0*(1.0       *(1.0-X)**4 -  4.0 *X    *(1.0-X)**3)
      BD5_R5_X =  1.0*(                      -  5.0       *(1.0-X)**4)
      ! Combine derivatives of basis linearly
      BD5_X    = A1*BD5_R0_X + A2*BD5_R1_X + A3*BD5_R2_X + A4*BD5_R3_X + A5*BD5_R4_X + A6*BD5_R5_X
      
!      BD5_R0_X = 
      
!      BD5   = A1*X +     A2*X**2
!      BD5_X = A1   + 2.0*A2*X
      
      RETURN
      END
! <-- mlcode