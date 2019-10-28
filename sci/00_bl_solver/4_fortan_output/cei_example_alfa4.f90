PROGRAM CEI_EXAMPLE
! TRILINEAR INTERPOLATION OF THE ENERGY INTERACTION COEFFICIENT
! THIS FUNCTION VERSION IS VALIDATED:
!     + for the search algorithm
!     + for the trilinear interpolation
!     + for the hk, rt and tppt derivatives
!
! Gael de Oliveira, Ricardo Pereira, Feb. 2015
! Development Stage / All rights reserved
! We follow implicit type rules REAL(A-H,O-Z), INTEGER(I-N)
	
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
     
! Define CEI evaluation point
      HK   = 2
      RT   = 1000
      TPPT = 0.6
      
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

      WRITE(*,*) 'REQUESTED POINT:'
      WRITE(*,1100) HK  , RT  , TPPT
      WRITE(*,*) 'SEARCHED POINT:'
      WRITE(*,1100) HK_S, RT_S, TPPT_S 
      WRITE(*,*) 'START SEARCH:'
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
                  WRITE(*,1000) I_HK, I_RT, I_TPPT
                  WRITE(*,1100) HK_RANGE(1,I_HK) , RT_RANGE(1,I_RT), TPPT_RANGE(1,I_TPPT)

                  GOTO 500
                ENDIF
 300          CONTINUE
            ENDIF
 200      CONTINUE
        ENDIF
 100  CONTINUE
 
 
 500  WRITE(*,*) 'FOUND UPPER NEIGHBOUR!!!'
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Search complete, go on with trilinear Interpolation
 ! Store lower bound (_0 subscript at I-1 index)
      HK_0   = HK_RANGE(1,I_HK-1)
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

! Display some diagnostics now
      WRITE(*,*   ) 'TRILINEAR INTERPOLATION STEPS:'
!      WRITE(*,1000) I_HK, I_RT, I_TPPT
      WRITE(*,1100) HK_0, RT_0, TPPT_0
      WRITE(*,1100) HK_1, RT_1, TPPT_1
      WRITE(*,1200) HK_D, RT_D, TPPT_D
      WRITE(*,1300) V000, V100, V010, V110
      WRITE(*,1310) V001, V101, V011, V111
      WRITE(*,1400) C00 , C10 , C01 , C11 
      WRITE(*,*   ) 'RESULT:'
      WRITE(*,1600) CEI , CEI_HK, CEI_RT, CEI_TPPT
 
 1000 FORMAT('I_HK=', I12  , '   I_RT=', I12  , '   I_TPPT=', I12)
 1100 FORMAT('  HK=', F12.4, '     RT=', F12.1, '     TPPT=', F12.4)
 1200 FORMAT(' dHK=', F12.4, '    dRT=', F12.4, '    dTPPT=', F12.4)
 1300 FORMAT('V000=', F12.4, '   V100=', F12.4, '     V010=', F12.4 , '     V110=', F12.4)
 1310 FORMAT('V001=', F12.4, '   V101=', F12.4, '     V011=', F12.4 , '     V111=', F12.4)
 1400 FORMAT(' C00=', F12.4, '    C10=', F12.4, '      C01=', F12.4 , '      C11=', F12.4)
 1600 FORMAT(' CEI=', F12.4, ' CEI_HK=', F12.4, '   CEI_RT=', F12.8 , ' CEI_TPPT=', F12.4)

END PROGRAM