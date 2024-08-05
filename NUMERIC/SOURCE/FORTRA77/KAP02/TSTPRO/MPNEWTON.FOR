      PROGRAM TEST
C                                                12.12.1991  ( Dubois Guido )
C****************************************************************************
C                                                                           *
C   Testprogram for the subroutines NEWPSZ, NEWPOL and NEWMOD.              *
C   We compute approximations for the zeros of functions.                   *
C                                                                           *
C   We test with the following functions and zeros:                         *
C                                                                           *
C   F(X)=A*X**N+B*X+C                                                       *
C                                                                           *
C   with:  A=1.D0, B=1.D0, C=1.D-4, N=5                                     *
C                                                                           *
C[  SOLUTION:                                                              ]*
C[  =========                                                              ]*
C[                                                                         ]*
C[  STARTING VALUE: X0 =    .50    NMAX =  99     RELERR =  .0D+00         ]*
C[                  A =  .10D+01   B =  .10D+01   C =  .10D-03   N =  5    ]*
C[                                                                         ]*
C[  STARTING FUNCTIONAL VALUE F(X0)   =  .53135000000000D+00               ]*
C[  1st DERIVATIVE AT START   F'(X0)  =  .13125000000000D+01               ]*
C[  2nd DERIVATIVE AT START   F''(X0) =  .25000000000000D+01               ]*
C[                                                                         ]*
C[  STARTING FUNCTIONAL VALUE F(X0)   =  .53135000000000D+00               ]*
C[  1st DERIVATIVE AT START   F'(X0)  =  .13125000000000D+01               ]*
C[                                                                         ]*
C[                                                                         ]*
C[  ABSERR          X                    F           ITANZ IFEHL J         ]*
C[                                                                         ]*
C[  .5D-12 -.10000000000000D-03  .35535679115278D-20    3    2      NEWPSZ ]*
C[  .5D-12 -.10000000000000D-03  .35535679115278D-20    5    2   1  NEWMOD ]*
C[  .5D-12 -.10000000000000D-03  .00000000000000D+00    3    2      NEWPOL ]*
C                                                                           *
C   The results are displayed on screen and can be sent to a file.          *
C   The iterates of the subroutine NEWMOD are stored in a file              *
C   [TAPE4 (CYBER 175), FORT4 (IBM), or interactively (MS-DOS)].            *
C                                                                           *
C****************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A,B,C,AA(0:5),WORK(0:5),X0,ABSERR,RELERR,FKT,
     +                 FSTR,F2STR,PN,PNSTR
      INTEGER N,NMAX
      COMMON /AB/ N,A,B,C
      EXTERNAL FKT,FSTR,F2STR
C
C     Initialize test values, set up test examples and adjust parameters
C
      DATA NMAX /99/
      DATA X0 /.5D0/
      DATA ABSERR /.5D-12/
      DATA RELERR /0.0D0/
      DATA AA /6*0.D0/
      N=5
      A=1.D0
      B=1.D0
      C=1.D-4
C
C     Output of test data
C
      WRITE(*,100)
      WRITE(*,101) X0,NMAX,RELERR,A,B,C,N
C
C     Compute functional values and derivatives
C
      WRITE(*,105) FKT(X0),FSTR(X0)
      WRITE(*,106) F2STR(X0)
      AA(N)=A
      AA(1)=B
      AA(0)=C
C
      CALL NEWPOH(AA,N,X0,PN,PNSTR,WORK)
      WRITE(*,105) PN,PNSTR
      WRITE(*,110)
C
      CALL NEWPSZ(FKT,FSTR,X0,NMAX,ABSERR,RELERR,X,F,ITANZ,IFEHL)
      WRITE(*,115) ABSERR,X,F,ITANZ,IFEHL
C
      CALL NEWMOD(FKT,FSTR,F2STR,X0,NMAX,4,ABSERR,RELERR,X,F,J,
     +            ITANZ,IFEHL)
      WRITE(*,116) ABSERR,X,F,ITANZ,IFEHL,J
C
      CALL NEWPOL(AA,N,X0,NMAX,ABSERR,RELERR,X,PN,ITANZ,IFEHL,WORK)
      WRITE(*,117) ABSERR,X,PN,ITANZ,IFEHL
      STOP
C
  100 FORMAT(1X,'C[  SOLUTION:',T77,']*',/,
     +       1X,'C[  =========',T77,']*',/,1X,'C[',T77,']*')
  101 FORMAT(1X,'C[  STARTING VALUE: X0 = ',F6.2,4X,'NMAX = ',I3,
     +       5X,'RELERR = ',D7.1,T77,']*',/,
     +       1X,'C[',18X,'A = ',D8.2,3X,'B = ',D8.2,3X,'C = ',D8.2,
     +       3X,'N = ',I2,T77,']*',/,
     +       1X,'C[',T77,']*')
  105 FORMAT(1X,'C[  STARTING FUNCTIONAL VALUE',T32,'F(X0)   = ',
     +       D20.14,T77,']*',/,
     +       1X,'C[  1st DERIVATIVE AT START',T32,'F''(X0)  = ',
     +       D20.14,T77,']*')
  106 FORMAT(1X,'C[  2nd DERIVATIVE AT START',T32,'F''''(X0) = ',
     +       D20.14,T77,']*',/,1X,'C[',T77,']*')
  110 FORMAT(1X,'C[',T77,']*',/,
     +       1X,'C[',T77,']*',/,
     +       1X,'C[  ABSERR',10X,'X',20X,'F',11X,'ITANZ',
     +       ' IFEHL J',T77,']*',/,
     +       1X,'C[',T77,']*')
  115 FORMAT(1X,'C[',1X,D7.1,2(1X,D20.14),2X,I3,4X,I1,6X,'NEWPSZ',
     +       T77,']*')
  116 FORMAT(1X,'C[',1X,D7.1,2(1X,D20.14),2X,I3,4X,I1,2X,I2,2X,'NEWMOD',
     +       T77,']*')
  117 FORMAT(1X,'C[',1X,D7.1,2(1X,D20.14),2X,I3,4X,I1,6X,'NEWPOL',
     +       T77,']*')
      END
C
C
      DOUBLE PRECISION FUNCTION FKT(X)
C
C     function to be tested
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /AB/ N,A,B,C
      S=1.D0
      DO 10 I=1,N
         S=S*X
   10 CONTINUE
      FKT=A*S+B*X+C
      RETURN
      END
C
C
      DOUBLE PRECISION FUNCTION FSTR(X)
C
C     first derivative of  FKT(X)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /AB/ N,A,B,C
      S=1.D0
      DO 10 I=1,N-1
         S=S*X
   10 CONTINUE
      FSTR=N*A*S+B
      RETURN
      END
C
C
      DOUBLE PRECISION FUNCTION F2STR(X)
C
C     2nd derivative of FKT(X)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /AB/ N,A,B,C
      S=1.D0
      DO 10 I=1,N-2
         S=S*X
   10 CONTINUE
      F2STR=N*(N-1)*A*S
      RETURN
      END
