      PROGRAM TEST 
!                                                               8.2.93  
!***********************************************************************
!                                                                      *
! Test program for the  SUBROUTINE QUAROM                              *
!                                                                      *
!                       3.1                                            *
! Evaluate the integral  I  3.0*SIN(X) DX  approximately using the     *
!                       0.0                                            *
!                                                                      *
! Romberg method.                                                      *
!                                                                      *
! Test results:                                                        *
!                                                                      *
![ ROMBERG QUADRATURE                                                 ]*
![ ------------------                                                 ]*
![                                                                    ]*
![ A =  0.0000000000D+00                                              ]*
![ B =  3.1000000000D+00                                              ]*
![ EPS =  5.0000000000D-08                                            ]*
![ MAXIMAL NUMBER OF COLUMNS        =    4                            ]*
![ STEP SIZE                      H =  2.5833333333D-02               ]*
![ COMPUTED VALUE OF INTEGRAL       =  5.9974054508D+00               ]*
![ EXACT VALUE                      =  5.9974054508D+00               ]*
![ ERROR COMMITTED                  =  6.2172489379D-15               ]*
![ ERROR ESTIMATE                   =  3.7783109974D-12               ]*
![ ERROR INVOLVED IN ERROR ESTIMATE =  3.7720937485D-12               ]*
!                                                                      *
!***********************************************************************
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION EL (100) 
      EXTERNAL FUNK 
      PI = 4.0D0 * DATAN (1.0D0) 
      ALOES = - 3.0D+00 * (DCOS (3.1D+00) - 1.0D+00) 
!                                                                       
      A = 0.0D+00 
      B = 3.1D+00 
      EPS = 0.5D-07 
      N = 6 
      H = (B - A) / 15.0D+00 
!                                                                       
      CALL QUAROM (A, B, EPS, N, H, FUNK, EL, RES, FEH, IFEHL) 
      IF (IFEHL.EQ.0) THEN 
         DIFF = DABS (RES - ALOES) 
         DF = DABS (DIFF - FEH) 
         WRITE ( *, 100) A, B, EPS, N, H, RES, ALOES, DIFF, FEH, DF 
      ELSE 
         WRITE ( * , * ) 'ERROR: IFEHL=', IFEHL 
      ENDIF 
      STOP 
  100 FORMAT(1X, 'C[ ROMBERG QUADRATURE', T72, ']*', /,                 &
     &       1X, 'C[ ', 18('-'), T72, ']*', /,                          &
     &       1X, 'C[', T72, ']*', /,                                    &
     &       1X, 'C[ A = ', 1PD17.10 , T72, ']*', /,                    &
     &       1X, 'C[ B = ', D17.10, T72, ']*', /,                       &
     &       1X, 'C[ EPS = ', D17.10, T72, ']*', /,                     &
     &       1X, 'C[ MAXIMAL NUMBER OF COLUMNS        = ', I4,          &
     &           T72, ']*', /,                                          &
     &       1X, 'C[ STEP SIZE                      H = ', D17.10,      &
     &           T72, ']*', /,                                          &
     &       1X, 'C[ COMPUTED VALUE OF INTEGRAL       = ', D17.10,      &
     &           T72, ']*', /,                                          &
     &       1X, 'C[ EXACT VALUE                      = ', D17.10,      &
     &           T72, ']*', /,                                          &
     &       1X, 'C[ ERROR COMMITTED                  = ', D17.10,      &
     &           T72, ']*', /,                                          &
     &       1X, 'C[ ERROR ESTIMATE                   = ', D17.10,      &
     &           T72, ']*', /,                                          &
     &       1X, 'C[ ERROR INVOLVED IN ERROR ESTIMATE = ', D17.10,      &
     &           T72, ']*')                                             
      END PROGRAM TEST                              
!                                                                       
      DOUBLEPRECISION FUNCTION FUNK (X) 
      DOUBLEPRECISION X 
      FUNK = 3.D0 * DSIN (X) 
      RETURN 
      END FUNCTION FUNK                             
