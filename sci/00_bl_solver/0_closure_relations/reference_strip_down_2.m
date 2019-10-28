      SUBROUTINE HKIN( H, MSQ, HK, HK_H, HK_MSQ )
      REAL MSQ
!
!---- calculate kinematic shape parameter (assuming air)
!     (from Whitfield )
      HK     = (H - 0.29*MSQ)/(1.0 + 0.113*MSQ)
      HK_H   =  1.0          /(1.0 + 0.113*MSQ)
      HK_MSQ = (-.29 - 0.113*HK) / (1.0 + 0.113*MSQ)
!
      RETURN
      END


      SUBROUTINE CFRT(CFR, CF, A, CFR_A, CFR_CF, CFR_D, ETA)
      DOUBLE PRECISION AD, C1, C2, DFC, DC2, D
!
!---- Cross-flow skin friction (Johnston)
      AD= DBLE(A)
      C1= DBLE(CF)
      C2= DBLE(CFR)
      D = DBLE(ETA)
!
      IF (ABS(CF).GT.0.01) THEN
       CFR    = 0.
       CFR_CF = 0.
       CFR_A  = 0.
       CFR_D  = 0.
      ELSE
       DFC    = 1.D0 - 0.5D0*D*AD*C2/((C1*C1 + C2*C2)**(0.75D0))
       DC2    = - (C2 - D*AD*(C1*C1+C2*C2)**(0.25D0) + &
                   AD*(C1) ) / DFC
       C2     = C2 + DC2
       CFR    = C2
       CFR_A  = (D*(C1*C1+C2*C2)**(0.25D0) - (C1) ) / DFC
       CFR_CF = (0.5D0*D*AD*C1/(C1*C1+C2*C2)**(0.75D0) - A ) &
                / DFC
       CFR_D  = - ( AD*(C1*C1+C2*C2)**(0.25D0) ) / DFC
      ENDIF
!
      RETURN
      END



      SUBROUTINE DITR(HS, US, CFR, CFX, A, ST, DIR, DIR_HS, &
                      DIR_US, DIR_CFR, DIR_CF, DIR_ST, DIR_A)
!
!---- Cross-flow turbulent dissipation coefficient (2 CD/H*)
      crsq   = abs(cfr)**1.5
      DIR    =  (6.1707*crsq  + A*A*ST*ST*(1.0-US) ) * 2.0/HS
      DIR_HS = -(6.1707*crsq  + A*A*ST*ST*(1.0-US) ) * 2.0/HS**2
      DIR_US =  (             - A*A*ST*ST          ) * 2.0/HS
      DIR_CF =  (0.                                ) * 2.0/HS
      DIR_CFR=  (6.1707*1.5*sqrt(abs(cfr))         ) * 2.0/HS
      DIR_ST =  (              2.0*A*A*ST*(1.0-US) ) * 2.0/HS
      DIR_A  =  (             2.0*A*ST*ST*(1.0-US) ) * 2.0/HS
!
      RETURN
      CF= CFX
      IF (ABS(CFX).LT.1.E-5 .AND. ABS(CFR).LT.1.E-5) CF= 1.
      DIR    =  ( 0.5*CFR*CFR*US/CF + A*A*ST*ST*(1.0-US) ) * 2.0/HS
      DIR_HS = -( 0.5*CFR*CFR*US/CF + A*A*ST*ST*(1.0-US) ) * 2.0/HS**2
      DIR_US =  ( 0.5*CFR*CFR   /CF - A*A*ST*ST          ) * 2.0/HS
      DIR_CF =  (-0.5*CFR*CFR*US/(CF*CF)                 ) * 2.0/HS
      DIR_CFR=  (     CFR    *US/CF                      ) * 2.0/HS
      DIR_ST =  (                    2.0*A*A*ST*(1.0-US) ) * 2.0/HS
      DIR_A  =  (                   2.0*A*ST*ST*(1.0-US) ) * 2.0/HS
!
      RETURN
      END



      SUBROUTINE CFTGMERCHANT_old( HK, RT, MSQ, US, SV, CF, CF_HK, CF_RT , CF_MSQ, CF_US, CF_SV)
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
!      Equivalent Xfoil specific code
!      CALL CFT( HK, RT, MSQ, CFI, CFI_HK, CFI_RT, CFI_MSQ )
        CALL CFTCMERCHANT( US , SV, CFTC, CFTC_US, CFTC_SV)
        ! Ratio_Cfsuc is used to activate or deactive correction of cf due to suction
        CFTC    = RATIO_CFSUC * CFTC
        CFTC_US = RATIO_CFSUC * CFTC_US
        CFTC_SV = RATIO_CFSUC * CFTC_SV

!       Now Compose Cf expression
      CF     = CFI      + RATIO_CFSUC * CFTC
!       And its derivatives
      CF_HK  = CFI_HK
      CF_RT  = CFI_RT
      CF_MSQ = CFI_MSQ
      CF_US  = CFTC_US
      CF_SV  = CFTC_SV
!      WRITE (*,*) 'CF_SV = ' , CF_SV
!       IMPORTANT
!       Record that US_HK is not zero, so CF_HK_TOTAL = CF_HK + CF_US * US_HK
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

      SUBROUTINE CFTGMERCHANT( HK, RT, MSQ, US, SV, CF, CF_HK, CF_RT , CF_MSQ, CF_US, CF_SV)
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