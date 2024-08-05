      PROGRAM TEST
C                                                                    93/02/11
C*****************************************************************************
C                                                                            *
C     Testprogramm for subroutine  NLESYS.                                   *
C     Solve a nonlinear system of equations                                  *
C             F(1)(X(1),...,X(N)) = 0                                        *
C             F(2)(X(1),...,X(N)) = 0                                        *
C             -----------------------                                        *
C             F(M)(X(1),...,X(N)) = 0                                        *
C     with the damped Newton iteration method. Here the Jacobi matrix is     *
C     approximated by central difference quotients.                          *  
C                                                                            *
C     The test example produces:                                             *
C                                                                            *
C[  TEST EXAMPLE:                                                           ]*
C[  ==============                                                          ]*
C[  SYSTEM OF EQUATIONS TO BE SOLVED:                                       ]*
C[                                                                          ]*
C[          F(X) = X(0)**2 - X(1) - 1                  = 0                  ]*
C[          F(X) = (X(0)-2)**2 + (X(1)-0.5)**2 - 1     = 0                  ]*
C[                                                                          ]*
C[  STARTING VECTOR:       .000D+00  .000D+00                               ]*
C[                                                                          ]*
C[                                                                          ]*
C[  REQUIRED PARAMETERS:                                                    ]*
C[                                                                          ]*
C[    ATTENUATION =  4                                                      ]*
C[    MAX. 100 ITERATIONS                                                   ]*
C[    ERROR BOUND:  .100D-05                                                ]*
C[                                                                          ]*
C[  SOLUTION:                                                               ]*
C[                                                                          ]*
C[       1.06735D+00     1.39228D-01                                        ]*
C[                                                                          ]*
C[  TERMINATED AFTER   6 ITERATION STEPS                                    ]*
C                                                                            *
C     Weitere Testl„ufe sind durch Modifizierung der Werte in den            *
C     DATA-Anweisungen und EXTERNAL-SUBROUTINES jederzeit m”glich.           *
C                                                                            *
C*****************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (N=1,M=1,LDDF=M)
      EXTERNAL FX,DFX
      CHARACTER *80 FVONX(0:M)
      DOUBLE PRECISION X(0:N),DF(0:LDDF,0:N+1),WORK(0:M),S(0:N),
     +                 D(0:N),F(0:N)
C
C     Initialize test example. Exchanging example is possible.
C
      DATA MAXIT,KMAX,EPS/100,4,1.D-06/
      DATA (X(I),I=0,N)/0.0D0,0.0D0/
C
C     Prepare for output
C
      DATA (FVONX(I),I=0,M)/'X(0)**2 - X(1) - 1',
     +'(X(0)-2)**2 + (X(1)-0.5)**2 - 1'/
C
C     Put out test example
C
      WRITE(*,2000)
      WRITE(*,2010)  (FVONX(I),I=0,M)
      WRITE(*,2020)  (X(I),I=0,N)
      WRITE(*,2040)  KMAX,MAXIT,EPS
C
      CALL NLESYS(FX,M,N,.FALSE.,DFX,DF,LDDF,MAXIT,EPS,
     +            KMAX,0,X,F,RNORM2,IERR,D,S,WORK)
C
      IF (IERR.NE.0) THEN
C
C     Put out error message
C
        WRITE(*,2030) 'ERROR: IERR=',IERR
      ELSE
C
C     Put out solution
C
        WRITE(*,2050) (X(I),I=0,N)
        WRITE(*,2060) MAXIT
      ENDIF
      STOP
C
C
 2000 FORMAT (1X,'C[',2X,'TEST EXAMPLE:',T78,']*',/,
     +        1X,'C[',2X,14('='),T78,']*',/,
     +        1X,'C[',2X,'SYSTEM OF EQUATIONS TO BE SOLVED:',T78,']*'
     +           ,/,1X,'C[',T78,']*')
 2010 FORMAT (1X,'C[',10X,'F(X) = ',A35,' = 0',T78,']*')
 2020 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[',2X,'STARTING VECTOR:',5X,2(D10.3),T78,']*',/,
     +        1X,'C[',T78,']*')
 2030 FORMAT (1X,'C[',2X,A,I1,T78,']*')
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
C
      SUBROUTINE FX (X,N,F,M)
      DOUBLE PRECISION X(0:N),F(0:M)
      F(0)=X(0)*X(0)-X(1) -1.D0
      F(1)=(X(0)-2.D0)**2.D0+(X(1)-0.5D0)**2.D0-1.D0
      RETURN
      END
C
C
      SUBROUTINE DFX (X,M,N,DF,LDDF)
      DOUBLE PRECISION DF(0:LDDF,0:N), X(0:N)
      DF(0,0) = 2.D0*X(0)
      DF(0,1) = -1.D0
      DF(1,0) = 2.D0*X(0)-4.D0
      DF(1,1) = 2.D0*X(1)-1.D0
      RETURN
      END
