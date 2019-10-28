PROGRAM CEI_EXAMPLE
! Gael de Oliveira, Ricardo Pereira, Feb. 2015
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
	
! Define CEI evaluation point
      HK   = 8.0
      RT   = 1000.0
      TPPT = 0.5 
      
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
                  ! Store lower bound (_0 subscript at I-1 index)
                  HK_0   = HK_RANGE(1,I_HK-1)
                  RT_0   = RT_RANGE(1,I_RT-1)
                  TPPT_0 = TPPT_RANGE(1,I_RT-1)
                  ! Store upper bound (_1 subscript at I   index)
                  
                  
                  ! And leave loop
                  GOTO 500
                ENDIF
 300          CONTINUE
            ENDIF
 200      CONTINUE
        ENDIF
 100  CONTINUE
 
 
 500  WRITE(*,*) 'FOUND UPPER NEIGHBOUR!!!'
 
      WRITE(*,1000) I_HK, I_RT, I_TPPT
      WRITE(*,1100) HK_RANGE(1,I_HK) , RT_RANGE(1,I_RT), TPPT_RANGE(1,I_TPPT)
 
 1000 FORMAT('I_HK=', I12  , '   I_RT=', I12  , '   I_TPPT=', I12)
 1100 FORMAT('  HK=', F12.4, '     RT=', F12.1, '     TPPT=', F12.4)

END PROGRAM