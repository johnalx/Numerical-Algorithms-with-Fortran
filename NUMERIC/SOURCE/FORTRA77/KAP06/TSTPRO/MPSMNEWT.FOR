      PROGRAM TEST
C                                                                    10/09/87
C                                                ( Thomas Meuser )
C*****************************************************************************
C                                                                            *
C     Testprogramm for subroutine SMNEWT.                                    *
C     Solve a nonlinear system of equations                                  *
C             F(1)(X(1),...,X(N)) = 0                                        *
C             F(2)(X(1),...,X(N)) = 0                                        *
C             -----------------------                                        *
C             F(N)(X(1),...,X(N)) = 0                                        *
C     with damped Newton iterations.                                         *
C                                                                            *
C     The test example produces:                                             *
C                                                                            *
C[  TEST EXAMPLE:                                                           ]*
C[  ==============                                                          ]*
C[  SYSTEM OF EQUATIONS TO BE SOLVED:                                       ]*
C[                                                                          ]*
C[          F(X) = X(1)**2 - X(2) - 1                  = 0                  ]*
C[          F(X) = (X(1)-2)**2 + (X(2)-0.5)**2 - 1     = 0                  ]*
C[                                                                          ]*
C[  STARTING VECTOR:       .000D+00  .000D+00                               ]*
C[                                                                          ]*
C[                                                                          ]*
C[  THE ITERATION STEPS ARE SAVED ON TAPE 5 !                               ]*
C[                                                                          ]*
C[  REQUIRED PARAMETERS:                                                    ]*
C[                                                                          ]*
C[    ATTENUATION =  0                                                      ]*
C[    MAX. 100 ITERATIONS                                                   ]*
C[    ERROR BOUND:  .100D-05                                                ]*
C[                                                                          ]*
C[  SOLUTION:                                                               ]*
C[                                                                          ]*
C[       1.06735D+00     1.39228D-01                                        ]*
C[                                                                          ]*
C[  TERMINATED AFTER   6 ITERATION STEPS                                    ]*
C                                                                            *
C     Other test runs with different data are possible.                      *
C                                                                            *
C*****************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (N=2,IDF=N)
      EXTERNAL FX,DFX
      CHARACTER *80 FVONX(N)
      INTEGER IWORK(N)
      DOUBLE PRECISION X(N),DF(IDF,N),WORK(4*N),F(N)
C
C     Initialize test example, change if desired
C
      DATA ITAPE,MAXIT,KMAX,EPS/5,100,0,1D-06/
      DATA (X(I),I=1,N)/0.D0,0.D0/
C
C     Prepare for output
C
      DATA (FVONX(I),I=1,N)/'X(1)**2 - X(2) - 1',
     +                      '(X(1)-2)**2 + (X(2)-0.5)**2 - 1'/
C
C     Put out test example
C
      WRITE(*,2000)
      WRITE(*,2010)  (FVONX(I),I=1,N)
      WRITE(*,2020)  (X(I),I=1,N)
      WRITE(*,2030)  ITAPE
      WRITE(*,2040)  KMAX,MAXIT,EPS
C
      OPEN(ITAPE,FILE='TAPE5')
C
      CALL SMNEWT(FX,DFX,N,MAXIT,IFEHL,KMAX,ITAPE,
     +            EPS,RNORM2,F,X,DF,IDF,IWORK,WORK)
C
      CLOSE(ITAPE)
C
C     Put out error message in  IFEHL (1,2,3)
C
      IF(IFEHL.EQ.1) THEN
        WRITE(*,1000)
      ELSEIF(IFEHL.EQ.2) THEN
        WRITE(*,1100)
      ELSEIF(IFEHL.EQ.3) THEN
        WRITE(*,1200)
      ELSE
C
C     Put out solution
C
         WRITE(*,2050) (X(I),I=1,N)
         WRITE(*,2060) MAXIT
      ENDIF
      STOP
C
 1000 FORMAT (1X,'AFTER MAX STEPS ERROR BOUND WAS NOT REACHED!')
 1100 FORMAT (1X,'ERROR IN THE SOLVING THE LINEAR SYSTEM OF',
     +           'EQUATIONS (JACOBI MATRIX IS SINGULAR!)')
 1200 FORMAT (1X,'INCORRECT INPUT PARAMETER!')
 2000 FORMAT (1X,'C[',2X,'TEST EXAMPLE:',T78,']*',/,
     +        1X,'C[',2X,14('='),T78,']*',/,
     +        1X,'C[',2X,'SYSTEM OF EQUATIONS TO BE SOLVED:',T78,']*'
     +           ,/,1X,'C[',T78,']*')
 2010 FORMAT (1X,'C[',10X,'F(X) = ',A35,' = 0',T78,']*')
 2020 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[',2X,'STARTING VECTOR:',5X,2(D10.3),T78,']*',/,
     +        1X,'C[',T78,']*')
 2030 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[',2X,'THE ITERATION STEPS ARE SAVED ON TAPE',
     +           I2,1X,'!',T78,']*')
 2040 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[',2X,'REQUIRED PARAMETERS:',T78,']*',/,
     +        1X,'C[',T78,']*',/,
     +        1X,'C[',4X,'ATTENUATION =',I3,T78,']*',/,
     +        1X,'C[',4X,'MAX.',I4,1X,'ITERATIONS',T78,']*',/,
     +        1X,'C[',4X,'ERROR BOUND:',D10.3,T78,']*')
 2050 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[',2X,'SOLUTION:',T78,']*',/,
     +        1X,'C[',T78,']*',/,
     +        1X,'C[',2X,2(1X,1PD15.5),T78,']*')
 2060 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[',2X,'TERMINATED AFTER',I4,1X,
     +           'ITERATION STEPS ',T78,']*')
      END
C
C
      SUBROUTINE FX (N,X,F)
      DOUBLE PRECISION X(N),F(N)
      F(1)=X(1)*X(1)-X(2) -1.D0
      F(2)=(X(1)-2.D0)**2+(X(2)-0.5D0)**2-1.D0
      RETURN
      END
C
C
      SUBROUTINE DFX (N,X,DF,IDF)
      DOUBLE PRECISION DF(IDF,N), X(N)
      DF(1,1) = 2.D0*X(1)
      DF(1,2) = -1.D0
      DF(2,1) = 2.D0*X(1)-4.D0
      DF(2,2) = 2.D0*X(2)-1.D0
      RETURN
      END
