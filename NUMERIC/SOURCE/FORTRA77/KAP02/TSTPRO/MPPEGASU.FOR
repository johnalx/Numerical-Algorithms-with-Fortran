      PROGRAM TEST
C                                                          ( Thomas Meuser )
C***************************************************************************
C                                                                          *
C   Testprogram for the subroutine  PEGASU.                                *
C   Compute a zero of the function FKT in the interval [A,B].              *
C   Necessary assumption:  FKT(A) * FKT(B) < 0.0 .                         *
C                                                                          *
C   Test function: F(X)=COS(X)-X*EXP(X); test interval: [-1,1];            *
C   test results:                                                          *
C                                                                          *
C[  EXAMPLE SUBROUTINE PEGASU                                             ]*
C[  =========================                                             ]*
C[                                                                        ]*
C[  FUNCTION: F(X)= COS(X)-X*EXP(X)                                       ]*
C[  GIVEN INTERVAL [ -1.0,  1.0]                                          ]*
C[                                                                        ]*
C[  RESULTS:                                                              ]*
C[                                                                        ]*
C[  - ERROR PARAMETERS:  .500000D-14   .000000D+00                        ]*
C[  - FUNCTIONAL ZERO AT  X =  .51775736368246D+00                        ]*
C[  - FINAL INCL. INTERVAL [ .51775736368246D+00, .51775736368246D+00]    ]*
C[  - NUMBER OF ITERATIONS = 10                                           ]*
C[                                                                        ]*
C[  NOTE: ABSERR < 4.0 * MACHINE CONSTANT                                 ]*
C                                                                          *
C***************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION FKT,A,B,ABSERR,RECERR,XSI,X1,X2
      INTEGER NMAX,IFEHL,IANZ
      CHARACTER *20 FVONX
      EXTERNAL FKT
C
C     Initialize; define another example here if desired 
C
      DATA FVONX /'COS(X)-X*EXP(X)'/
      DATA A /-1.0D0/
      DATA B /1.0D0/
      DATA ABSERR /.5D-14/
      DATA RECERR /0.D0/
      DATA NMAX /99/
C
C     Output of the tested example
C
      WRITE(*,100)
      WRITE(*,110) FVONX
      WRITE(*,115) A,B
C
      CALL PEGASU (FKT,A,B,ABSERR,RECERR,NMAX,XSI,X1,X2,IANZ,IFEHL)
      WRITE(*,120) ABSERR,RECERR,XSI,X1,X2,IANZ
      IFEHL=IFEHL+3
C
C     Put out error message according to  IFEHL (-2/-1/0/1/2/3)
C
      GOTO (10,20,30,40,50,60) IFEHL
   10 WRITE(*,130)
      STOP
   20 WRITE(*,140)
      STOP
   30 WRITE(*,150)
      STOP
   40 WRITE(*,160)
      STOP
   50 WRITE(*,170)
      STOP
   60 WRITE(*,180)
      STOP
C
  100 FORMAT(1X,'C[  EXAMPLE SUBROUTINE PEGASU',T76,']*',/,
     +       1X,'C[  ',25('='),T76,']*',/,
     +       1X,'C[',T76,']*')
  110 FORMAT(1X,'C[  FUNCTION: F(X)= ',A,T76,']*')
  115 FORMAT(1X,'C[  GIVEN INTERVAL [',F5.1,',',F5.1,']',T76,']*',/,
     +       1X,'C[',T76,']*')
  120 FORMAT(1X,'C[  RESULTS:',T76,']*',/,
     +       1X,'C[',T76,']*',/,
     +       1X,'C[  - ERROR PARAMETERS: ',D12.6,2X,D12.6,T76,']*',/,
     +       1X,'C[  - FUNCTIONAL ZERO AT  X = ',D20.14,T76,']*',/,
     +       1X,'C[  - FINAL INCL. INTERVAL [',D20.14,',',D20.14,']',
     +          T76,']*',/,
     +       1X,'C[  - NUMBER OF ITERATIONS = ',I2,T76,']*',/,
     +       1X,'C[',T76,']*')
  130 FORMAT(1X,'C[  TERMINATED DUE TO INVALID ERROR BOUNDS',
     +          T76,']*')
  140 FORMAT(1X,'C[  TERMINATED SINCE F(A)*F(A)>=0.0.',T76,']*')
  150 FORMAT(1X,'C[  A OR B ARE ZEROS:',T76,']*')
  160 FORMAT(1X,'C[  NOTE: ABSERR < 4.0 * MACHINE CONSTANT',T76,']*')
  170 FORMAT(1X,'C[  NOTE: BREAK OFF CRITERION MET',T76,']*')
  180 FORMAT(1X,'C[  TERMINATED DUE TO EXCESSIVELY MANY STEPS',T76,']*')
      END
C
C
      DOUBLE PRECISION FUNCTION FKT(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      FKT=DCOS(X)-X*DEXP(X)
      RETURN
      END
