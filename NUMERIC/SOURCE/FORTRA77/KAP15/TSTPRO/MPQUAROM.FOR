      PROGRAM TEST
C                                                               8.2.93
C***********************************************************************
C                                                                      *
C Test program for the  SUBROUTINE QUAROM                              *
C                                                                      *
C                       3.1                                            *
C Evaluate the integral  I  3.0*SIN(X) DX  approximately using the     *
C                       0.0                                            *
C                                                                      *
C Romberg method.                                                      *
C                                                                      *
C Test results:                                                        *
C                                                                      *
C[ ROMBERG QUADRATURE                                                 ]*
C[ ------------------                                                 ]*
C[                                                                    ]*
C[ A =  0.0000000000D+00                                              ]*
C[ B =  3.1000000000D+00                                              ]*
C[ EPS =  5.0000000000D-08                                            ]*
C[ MAXIMAL NUMBER OF COLUMNS        =    4                            ]*
C[ STEP SIZE                      H =  2.5833333333D-02               ]*
C[ COMPUTED VALUE OF INTEGRAL       =  5.9974054508D+00               ]*
C[ EXACT VALUE                      =  5.9974054508D+00               ]*
C[ ERROR COMMITTED                  =  6.2172489379D-15               ]*
C[ ERROR ESTIMATE                   =  3.7783109974D-12               ]*
C[ ERROR INVOLVED IN ERROR ESTIMATE =  3.7720937485D-12               ]*
C                                                                      *
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EL(100)
      EXTERNAL FUNK
      PI=4.0D0*DATAN(1.0D0)
      ALOES=-3.0D+00*(DCOS(3.1D+00)-1.0D+00)
C
      A=0.0D+00
      B=3.1D+00
      EPS=0.5D-07
      N=6
      H=(B-A)/15.0D+00
C
      CALL QUAROM(A,B,EPS,N,H,FUNK,EL,RES,FEH,IFEHL)
      IF (IFEHL .EQ. 0) THEN
        DIFF=DABS(RES-ALOES)
        DF=DABS(DIFF-FEH)
        WRITE(*,100) A,B,EPS,N,H,RES,ALOES,DIFF,FEH,DF
         ELSE
           WRITE(*,*) 'ERROR: IFEHL=',IFEHL
         ENDIF
      STOP
  100 FORMAT(1X, 'C[ ROMBERG QUADRATURE', T72, ']*', /,
     +       1X, 'C[ ', 18('-'), T72, ']*', /,
     +       1X, 'C[', T72, ']*', /,
     +       1X, 'C[ A = ', 1PD17.10 , T72, ']*', /,
     +       1X, 'C[ B = ', D17.10, T72, ']*', /,
     +       1X, 'C[ EPS = ', D17.10, T72, ']*', /,
     +       1X, 'C[ MAXIMAL NUMBER OF COLUMNS        = ', I4,
     +           T72, ']*', /,
     +       1X, 'C[ STEP SIZE                      H = ', D17.10,
     +           T72, ']*', /,
     +       1X, 'C[ COMPUTED VALUE OF INTEGRAL       = ', D17.10,
     +           T72, ']*', /,
     +       1X, 'C[ EXACT VALUE                      = ', D17.10,
     +           T72, ']*', /,
     +       1X, 'C[ ERROR COMMITTED                  = ', D17.10,
     +           T72, ']*', /,
     +       1X, 'C[ ERROR ESTIMATE                   = ', D17.10, 
     +           T72, ']*', /,
     +       1X, 'C[ ERROR INVOLVED IN ERROR ESTIMATE = ', D17.10,
     +           T72, ']*')
      END
C
      DOUBLE PRECISION FUNCTION FUNK (X)
      DOUBLE PRECISION X
      FUNK = 3.D0*DSIN(X)
      RETURN
      END
