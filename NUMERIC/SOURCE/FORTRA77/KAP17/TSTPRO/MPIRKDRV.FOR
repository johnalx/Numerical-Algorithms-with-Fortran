C                                                (Thomas Meuser)     7/12/87
C*****************************************************************************
C                                                                            *
C     Test program for the subroutines IRKDRV, IRKCOE and IMRUKU.            *
C     Solve an initial value problem for a system of ordinary DE's:          *
C                                                                            *
C                     Y(1)' = F1(X,Y(1),Y(2),...,Y(N))                       *
C                     Y(2)' = F2(X,Y(1),Y(2),...,Y(N))                       *
C                                  ...                                       *
C                     Y(N)' = FN(X,Y(1),Y(2),...,Y(N))                       *
C                                                                            *
C     using the  RUNGE-KUTTA method.                                         *
C                                                                            *
C     Test results:                                                          *
C[                                                                          ]*
C[ TEST EXAMPLE:                                                            ]*
C[ =============                                                            ]*
C[ GIVEN SYSTEM OF DIFFERENTIAL EQUATIONS:                                  ]*
C[ ---------------------------------------                                  ]*
C[   Y( 1)' = -Y2                           WITH WEIGHT G =   .100D+01      ]*
C[   Y( 2)' = Y1                            WITH WEIGHT G =   .100D+01      ]*
C[   Y( 3)' = Y5 - Y2*Y3                    WITH WEIGHT G =   .100D+01      ]*
C[   Y( 4)' = Y1*Y4 - Y6                    WITH WEIGHT G =   .100D+01      ]*
C[   Y( 5)' = -Y3 - Y2*Y5                   WITH WEIGHT G =   .100D+01      ]*
C[   Y( 6)' = Y4 + Y1*Y6                    WITH WEIGHT G =   .100D+01      ]*
C[                                                                          ]*
C[ XK =   .000D+00 THE FOLLOWING INITIAL VALUES ARE PROVIDED:               ]*
C[                                                                          ]*
C[   Y( 1) =   .100D+01                                                     ]*
C[   Y( 2) =   .000D+00                                                     ]*
C[   Y( 3) =   .000D+00                                                     ]*
C[   Y( 4) =   .100D+01                                                     ]*
C[   Y( 5) =   .272D+01                                                     ]*
C[   Y( 6) =   .000D+00                                                     ]*
C[                                                                          ]*
C[ REQUIRED PARAMETERS:                                                     ]*
C[ --------------------                                                     ]*
C[ -ABSOLUTE ERROR =     -.1000000000D+01                                   ]*
C[ -RELATIVE ERROR =      .1000000000D-07                                   ]*
C[ -THE NODES FOR IRKV ARE STORED IN FILE 9                                 ]*
C[ -BEFORE SOLVING THE IVP A NODE FILE IS CREATED BY IRKCOE                 ]*
C[                                                                          ]*
C[ SOLUTION:                                                                ]*
C[ ---------                                                                ]*
C[ AT  XENDE=  .628D+01                                                     ]*
C[ THE FUNCTION VALUES ARE COMPUTED AS FOLLOWS:                             ]*
C[                                                                          ]*
C[   Y( 1) =      .1000000000D+01                                           ]*
C[   Y( 2) =      .5449934193D-10                                           ]*
C[   Y( 3) =     -.5372572278D-08                                           ]*
C[   Y( 4) =      .1000000003D+01                                           ]*
C[   Y( 5) =      .2718281828D+01                                           ]*
C[   Y( 6) =     -.1845522259D-08                                           ]*
C[                                                                          ]*
C[ PROCESS TERMINATED AS PLANNED!                                           ]*
C                                                                            *
C  Remark: For the given initial values  Y(0)=(1,0,0,1,e,0), the exact       *
C          solution is given by:                                             *
C                                                                            *
C          Y1 = COS(X)                                                       *
C          Y2 = SIN(X)                                                       *
C          Y3 = SIN(X) * EXP( COS(X) )                                       *
C          Y4 = COS(X) * EXP( SIN(X) )                                       *
C          Y5 = COS(X) * EXP( COS(X) )                                       *
C          Y6 = SIN(X) * EXP( SIN(X) )                                       *
C                                                                            *
C          The test problem is designed to compute the function values       *
C          at XENDE = 2*PHI.                                                 *
C                                                                            *
C     Further tests can be performed for other data and  EXTERNAL-SUBROUTINES*
C                                                                            *
C*****************************************************************************
C
      PROGRAM TEST
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (N=6,NMAX=12)
      DIMENSION G(N),Y0(N),YQ(N),WORK1(8*NMAX+5*N-2),WORK21(N,N),
     +          WORK22(N+NMAX-1,NMAX-1),WORK23(2*N+NMAX,NMAX)
      EXTERNAL DGL
      CHARACTER *80 YSTRI(N)
C
C     Initialize data. Use different data idf desired.
C     Problem specific we use  XENDE = 2*PHI and Y0(5) = e = 2.71... her.
C
      DATA X0,XENDE,Y0/0.D0,0.D0,1.D0,0.D0,0.D0,1.D0,0.D0,0.D0/
      DATA G/6*1.D0/
      DATA IWAEHL,ISTDAT,IAUS,IPROT/1,9,0,0/
      DATA EPSM,EPS/-1.0D00,1.D-08/
      XENDE = 8.D0*DATAN(1.D0)
      Y0(5) = DEXP(1.D0)
C
C     Prepare for output
C
      DATA (YSTRI(I),I=1,N)/'-Y2','Y1','Y5 - Y2*Y3','Y1*Y4 - Y6',
     +                    '-Y3 - Y2*Y5','Y4 + Y1*Y6'/
C
C     output of test example
C
      WRITE(*,2000)
      WRITE(*,2010)  (I,YSTRI(I),G(I),I=1,N)
      WRITE(*,2020)  X0,(I,Y0(I),I=1,N)
      WRITE(*,2030)  EPSM,EPS,ISTDAT
      IF (IAUS.GT.0) WRITE (*,2040) IAUS
      IF (IPROT.GT.0) WRITE (*,2045) IPROT
      IF (IWAEHL .EQ. 0) THEN
         WRITE(*,2050)
      ELSE
         WRITE(*,2055)
      ENDIF
C
C
      CALL IRKDRV (DGL,N,NMAX,IWAEHL,ISTDAT,IAUS,IPROT,IFEHL,EPSM,EPS,
     +          G,X0,XENDE,Y0,YQ,WORK1,WORK21,WORK22,WORK23)
C
C
      IF (IFEHL .EQ. 1) THEN
C
C     put out error code in  IFEHL (0/1/2)
C
          WRITE(*,2080)'ERROR: NO CONVERGENCE; TRY TO INCREASE NMAX ?'
      ELSEIF (IFEHL .EQ. 2) THEN
          WRITE(*,2080)'ERROR: INCORRECT INPUT PARAMETERS!'
      ELSE
C
C     put out solution
C
          WRITE(*,2060)XENDE,(I,YQ(I),I=1,N)
          WRITE(*,2080)'PROCESS TERMINATED AS PLANNED!'
      ENDIF
      STOP
C
C
 2000 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[ TEST EXAMPLE:',T78,']*',/,
     +        1X,'C[ ',13('='),T78,']*',/,
     +        1X,'C[ GIVEN SYSTEM OF DIFFERENTIAL EQUATIONS:',
     +           T78,']*',/,
     +        1X,'C[ ',39('-'),T78,']*')
 2010 FORMAT (1X,'C[',3X,'Y(',I2,')'' = ',A30,'WITH WEIGHT G = ',
     *           D10.3,T78,']*')
 2020 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[ XK = ',D10.3,' THE FOLLOWING INITIAL VALUES ',
     +           'ARE PROVIDED:',T78,']*',/,
     +         1X,'C[',T78,']*',/,
     +        (1X,'C[',3X,'Y(',I2,') = ',D10.3,T78,']*'))
 2030 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[ REQUIRED PARAMETERS:',T78,']*',/,
     +        1X,'C[ ',20('-'),T78,']*',/,
     +        1X,'C[ -ABSOLUTE ERROR = ',D20.10,T78,']*',/,
     +        1X,'C[ -RELATIVE ERROR = ',D20.10,T78,']*',/,
     +        1X,'C[ -THE NODES FOR IRKV ARE STORED IN FILE ',I1,
     +           T78,']*')
 2040 FORMAT (1X,'C[ -OUTPUT OF INTERIM SOLUTION IS STORED IN FILE',I1,
     +           T78,']*')
 2045 FORMAT (1X,'C[ -A RUNNING PROTOCOL IS KEPT IN FILE ',I1,T78,']*')
 2050 FORMAT (1X,'C[ -ONLY SOLUTION OF THE IVP AT THE NODES MUST ',
     +           'BE PRESENT',T78,']*')
 2055 FORMAT (1X,'C[ -BEFORE SOLVING THE IVP A NODE FILE IS CREATED ',
     +           'BY IRKCOE',T78,']*')
 2060 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[ SOLUTION:',T78,']*',/,
     +        1X,'C[ ',9('-'),T78,']*',/,
     +        1X,'C[ AT  XENDE=',D10.3,T78,']*',/,
     +        1X,'C[ THE FUNCTION VALUES ARE COMPUTED AS FOLLOWS:',
     +           T78,']*',/,
     +        1X,'C[',T78,']*',/,
     +        (1X,'C[',3X,'Y(',I2,') = ',D20.10,T78,']*'))
 2080 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[',1X,A,T78,']*')
      END
C
C
      SUBROUTINE DGL (N,X,Y,F)
      DOUBLE PRECISION X,Y(N),F(N)
      F(1)=-Y(2)
      F(2)=Y(1)
      F(3)=Y(5)-Y(2)*Y(3)
      F(4)=Y(1)*Y(4)-Y(6)
      F(5)=-Y(3)-Y(2)*Y(5)
      F(6)=Y(4)+Y(1)*Y(6)
      RETURN
      END
