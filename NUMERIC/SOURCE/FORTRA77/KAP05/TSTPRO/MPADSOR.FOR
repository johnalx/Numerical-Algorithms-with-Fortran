      PROGRAM TEST
C
C*****************************************************************************
C                                                                            *
C     Testprogram for the subroutine  ADSOR.                                 *
C     Solve an inhomogeneous system of linear equations  A*X = B with        *
C     det A  not zero iteratively via the SOR method.                        *
C                                                                            *
C     The output for the test example is as follows:                         *
C                                                                            *
C[  EXAMPLE:                                                                ]*
C[  ========                                                                ]*
C[  COEFFICIENT MATRIX A:                                                   ]*
C[                                                                          ]*
C[   .40000D+01  -.10000D+01  -.10000D+01   .00000D+00                      ]*
C[  -.10000D+01   .40000D+01   .00000D+00  -.10000D+01                      ]*
C[  -.10000D+01   .00000D+00   .40000D+01  -.10000D+01                      ]*
C[   .00000D+00  -.10000D+01  -.10000D+01   .40000D+01                      ]*
C[                                                                          ]*
C[  RIGHT HAND SIDE:                                                        ]*
C[                                                                          ]*
C[   .00000D+00   .00000D+00   .10000D+04   .10000D+04                      ]*
C[                                                                          ]*
C[  STARTING VECTOR X:                                                      ]*
C[                                                                          ]*
C[   .00000D+00   .00000D+00   .00000D+00   .00000D+00                      ]*
C[                                                                          ]*
C[                                                                          ]*
C[                                                                          ]*
C[                                                                          ]*
C[  INPUT PARAMETERS:                                                       ]*
C[    KADAPT =    4                                                         ]*
C[    EPS =  .10000000000000D-08                                            ]*
C[    KMAX =  100                                                           ]*
C[    IMETH  =    0                                                         ]*
C[    ISWITC =    0                                                         ]*
C[    OMEGA =  .10000000000000D+01                                          ]*
C[                                                                          ]*
C[                                                                          ]*
C[  OUTPUT PARAMETERS:                                                      ]*
C[    OMEGA =  .10739144783362D+01                                          ]*
C[    ITNUMB =   13                                                         ]*
C[    IERR =    1                                                           ]*
C[                                                                          ]*
C[                                                                          ]*
C[  SOLUTION VECTOR X:                                                      ]*
C[   .12500D+03   .12500D+03   .37500D+03   .37500D+03                      ]*
C[                                                                          ]*
C[  RESIDUAL VECTOR RES:                                                    ]*
C[   .79888D-08  -.30867D-09  -.30867D-09  -.53348D-09                      ]*
C[                                                                          ]*
C[                                                                          ]*
C[  INPUT PARAMETERS:                                                       ]*
C[    KADAPT =    4                                                         ]*
C[    EPS =  .10000000000000D-08                                            ]*
C[    KMAX =  100                                                           ]*
C[    IMETH  =    1                                                         ]*
C[    ISWITC =    0                                                         ]*
C[    OMEGA =  .10739144783362D+01                                          ]*
C[                                                                          ]*
C[                                                                          ]*
C[  OUTPUT PARAMETERS:                                                      ]*
C[    OMEGA =  .10739144783362D+01                                          ]*
C[    ITNUMB =   10                                                         ]*
C[    IERR =    1                                                           ]*
C[                                                                          ]*
C[                                                                          ]*
C[  SOLUTION VECTOR X:                                                      ]*
C[   .12500D+03   .12500D+03   .37500D+03   .37500D+03                      ]*
C[                                                                          ]*
C[  RESIDUAL VECTOR RES:                                                    ]*
C[   .15031D-07  -.56438D-08  -.44269D-08   .10584D-09                      ]*
C[                                                                          ]*
C[                                                                          ]*
C[  INPUT PARAMETERS:                                                       ]*
C[    KADAPT =    4                                                         ]*
C[    EPS =  .10000000000000D-08                                            ]*
C[    KMAX =  100                                                           ]*
C[    IMETH  =    2                                                         ]*
C[    ISWITC =    0                                                         ]*
C[    OMEGA =  .10739144783362D+01                                          ]*
C[                                                                          ]*
C[                                                                          ]*
C[  OUTPUT PARAMETERS:                                                      ]*
C[    OMEGA =  .10000000000000D+01                                          ]*
C[    ITNUMB =   17                                                         ]*
C[    IERR =    1                                                           ]*
C[                                                                          ]*
C[                                                                          ]*
C[  SOLUTION VECTOR X:                                                      ]*
C[   .12500D+03   .12500D+03   .37500D+03   .37500D+03                      ]*
C[                                                                          ]*
C[  RESIDUAL VECTOR RES:                                                    ]*
C[   .43656D-07   .10914D-07   .10914D-07   .00000D+00                      ]*
C[                                                                          ]*
C[                                                                          ]*
C[  INPUT PARAMETERS:                                                       ]*
C[    KADAPT =    4                                                         ]*
C[    EPS =  .10000000000000D-08                                            ]*
C[    KMAX =  100                                                           ]*
C[    IMETH  =    0                                                         ]*
C[    ISWITC =    1                                                         ]*
C[    OMEGA =  .10000000000000D+01                                          ]*
C[                                                                          ]*
C[                                                                          ]*
C[  OUTPUT PARAMETERS:                                                      ]*
C[    OMEGA =  .10739144783362D+01                                          ]*
C[    ITNUMB =   13                                                         ]*
C[    IERR =    1                                                           ]*
C[                                                                          ]*
C[                                                                          ]*
C[  SOLUTION VECTOR X:                                                      ]*
C[   .12500D+03   .12500D+03   .37500D+03   .37500D+03                      ]*
C[                                                                          ]*
C[  RESIDUAL VECTOR RES:                                                    ]*
C[   .79888D-08  -.30867D-09  -.30867D-09  -.53348D-09                      ]*
C[                                                                          ]*
C[                                                                          ]*
C[  INPUT PARAMETERS:                                                       ]*
C[    KADAPT =    4                                                         ]*
C[    EPS =  .10000000000000D-08                                            ]*
C[    KMAX =  100                                                           ]*
C[    IMETH  =    1                                                         ]*
C[    ISWITC =    1                                                         ]*
C[    OMEGA =  .10739144783362D+01                                          ]*
C[                                                                          ]*
C[                                                                          ]*
C[  OUTPUT PARAMETERS:                                                      ]*
C[    OMEGA =  .10739144783362D+01                                          ]*
C[    ITNUMB =   10                                                         ]*
C[    IERR =    1                                                           ]*
C[                                                                          ]*
C[                                                                          ]*
C[  SOLUTION VECTOR X:                                                      ]*
C[   .12500D+03   .12500D+03   .37500D+03   .37500D+03                      ]*
C[                                                                          ]*
C[  RESIDUAL VECTOR RES:                                                    ]*
C[   .15031D-07  -.56438D-08  -.44269D-08   .10584D-09                      ]*
C[                                                                          ]*
C[                                                                          ]*
C[  INPUT PARAMETERS:                                                       ]*
C[    KADAPT =    4                                                         ]*
C[    EPS =  .10000000000000D-08                                            ]*
C[    KMAX =  100                                                           ]*
C[    IMETH  =    2                                                         ]*
C[    ISWITC =    1                                                         ]*
C[    OMEGA =  .10739144783362D+01                                          ]*
C[                                                                          ]*
C[                                                                          ]*
C[  OUTPUT PARAMETERS:                                                      ]*
C[    OMEGA =  .10000000000000D+01                                          ]*
C[    ITNUMB =   17                                                         ]*
C[    IERR =    1                                                           ]*
C[                                                                          ]*
C[                                                                          ]*
C[  SOLUTION VECTOR X:                                                      ]*
C[   .12500D+03   .12500D+03   .37500D+03   .37500D+03                      ]*
C[                                                                          ]*
C[  RESIDUAL VECTOR RES:                                                    ]*
C[   .43656D-07   .10914D-07   .10914D-07   .00000D+00                      ]*
C[                                                                          ]*
C[                                                                          ]*
C[  INPUT PARAMETERS:                                                       ]*
C[    KADAPT =    4                                                         ]*
C[    EPS =  .10000000000000D-08                                            ]*
C[    KMAX =  100                                                           ]*
C[    IMETH  =    0                                                         ]*
C[    ISWITC =    2                                                         ]*
C[    OMEGA =  .10000000000000D+01                                          ]*
C[                                                                          ]*
C[                                                                          ]*
C[  OUTPUT PARAMETERS:                                                      ]*
C[    OMEGA =  .10739144783362D+01                                          ]*
C[    ITNUMB =   13                                                         ]*
C[    IERR =    1                                                           ]*
C[                                                                          ]*
C[                                                                          ]*
C[  SOLUTION VECTOR X:                                                      ]*
C[   .12500D+03   .12500D+03   .37500D+03   .37500D+03                      ]*
C[                                                                          ]*
C[  RESIDUAL VECTOR RES:                                                    ]*
C[   .79888D-08  -.30867D-09  -.30867D-09  -.53348D-09                      ]*
C[                                                                          ]*
C[                                                                          ]*
C[  INPUT PARAMETERS:                                                       ]*
C[    KADAPT =    4                                                         ]*
C[    EPS =  .10000000000000D-08                                            ]*
C[    KMAX =  100                                                           ]*
C[    IMETH  =    1                                                         ]*
C[    ISWITC =    2                                                         ]*
C[    OMEGA =  .10739144783362D+01                                          ]*
C[                                                                          ]*
C[                                                                          ]*
C[  OUTPUT PARAMETERS:                                                      ]*
C[    OMEGA =  .10739144783362D+01                                          ]*
C[    ITNUMB =   10                                                         ]*
C[    IERR =    1                                                           ]*
C[                                                                          ]*
C[                                                                          ]*
C[  SOLUTION VECTOR X:                                                      ]*
C[   .12500D+03   .12500D+03   .37500D+03   .37500D+03                      ]*
C[                                                                          ]*
C[  RESIDUAL VECTOR RES:                                                    ]*
C[   .15031D-07  -.56438D-08  -.44269D-08   .10584D-09                      ]*
C[                                                                          ]*
C[                                                                          ]*
C[  INPUT PARAMETERS:                                                       ]*
C[    KADAPT =    4                                                         ]*
C[    EPS =  .10000000000000D-08                                            ]*
C[    KMAX =  100                                                           ]*
C[    IMETH  =    2                                                         ]*
C[    ISWITC =    2                                                         ]*
C[    OMEGA =  .10739144783362D+01                                          ]*
C[                                                                          ]*
C[                                                                          ]*
C[  OUTPUT PARAMETERS:                                                      ]*
C[    OMEGA =  .10000000000000D+01                                          ]*
C[    ITNUMB =   17                                                         ]*
C[    IERR =    1                                                           ]*
C[                                                                          ]*
C[                                                                          ]*
C[  SOLUTION VECTOR X:                                                      ]*
C[   .12500D+03   .12500D+03   .37500D+03   .37500D+03                      ]*
C[                                                                          ]*
C[  RESIDUAL VECTOR RES:                                                    ]*
C[   .43656D-07   .10914D-07   .10914D-07   .00000D+00                      ]*
C[                                                                          ]*
C[                                                                          ]*
C[  INPUT PARAMETERS:                                                       ]*
C[    KADAPT =    4                                                         ]*
C[    EPS =  .10000000000000D-08                                            ]*
C[    KMAX =  100                                                           ]*
C[    IMETH  =    0                                                         ]*
C[    ISWITC =    3                                                         ]*
C[    OMEGA =  .10000000000000D+01                                          ]*
C[                                                                          ]*
C[                                                                          ]*
C[  OUTPUT PARAMETERS:                                                      ]*
C[    OMEGA =  .10739144783362D+01                                          ]*
C[    ITNUMB =   13                                                         ]*
C[    IERR =    1                                                           ]*
C[                                                                          ]*
C[                                                                          ]*
C[  SOLUTION VECTOR X:                                                      ]*
C[   .12500D+03   .12500D+03   .37500D+03   .37500D+03                      ]*
C[                                                                          ]*
C[  RESIDUAL VECTOR RES:                                                    ]*
C[   .79888D-08  -.30867D-09  -.30867D-09  -.53348D-09                      ]*
C[                                                                          ]*
C[                                                                          ]*
C[  INPUT PARAMETERS:                                                       ]*
C[    KADAPT =    4                                                         ]*
C[    EPS =  .10000000000000D-08                                            ]*
C[    KMAX =  100                                                           ]*
C[    IMETH  =    1                                                         ]*
C[    ISWITC =    3                                                         ]*
C[    OMEGA =  .10739144783362D+01                                          ]*
C[                                                                          ]*
C[                                                                          ]*
C[  OUTPUT PARAMETERS:                                                      ]*
C[    OMEGA =  .10739144783362D+01                                          ]*
C[    ITNUMB =   10                                                         ]*
C[    IERR =    1                                                           ]*
C[                                                                          ]*
C[                                                                          ]*
C[  SOLUTION VECTOR X:                                                      ]*
C[   .12500D+03   .12500D+03   .37500D+03   .37500D+03                      ]*
C[                                                                          ]*
C[  RESIDUAL VECTOR RES:                                                    ]*
C[   .15031D-07  -.56438D-08  -.44269D-08   .10584D-09                      ]*
C[                                                                          ]*
C[                                                                          ]*
C[  INPUT PARAMETERS:                                                       ]*
C[    KADAPT =    4                                                         ]*
C[    EPS =  .10000000000000D-08                                            ]*
C[    KMAX =  100                                                           ]*
C[    IMETH  =    2                                                         ]*
C[    ISWITC =    3                                                         ]*
C[    OMEGA =  .10739144783362D+01                                          ]*
C[                                                                          ]*
C[                                                                          ]*
C[  OUTPUT PARAMETERS:                                                      ]*
C[    OMEGA =  .10000000000000D+01                                          ]*
C[    ITNUMB =   17                                                         ]*
C[    IERR =    1                                                           ]*
C[                                                                          ]*
C[                                                                          ]*
C[  SOLUTION VECTOR X:                                                      ]*
C[   .12500D+03   .12500D+03   .37500D+03   .37500D+03                      ]*
C[                                                                          ]*
C[  RESIDUAL VECTOR RES:                                                    ]*
C[   .43656D-07   .10914D-07   .10914D-07   .00000D+00                      ]*
C                                                                            *
C*****************************************************************************
C
C
      DOUBLE PRECISION  A(1:4,1:4), B(1:4), X(1:4),
     *                  RES(1:4), WORK(1:4,1:3),EPS,OMEGA
      OMEGA=1.0D+00
C
C output of the test example
C
      CALL INIT (A,B,X,KADAPT,EPS,KMAX)
      WRITE(*,2000)
      DO 10 I=1,4
        WRITE(*,2010)(A(I,J),J=1,4)
10    CONTINUE
      WRITE(*,2020)(B(I),I=1,4)
      WRITE(*,2030)(X(I),I=1,4)
C
      DO 100 ISWITC=0,3
        DO 200 IMETH=0,2
          CALL INIT (A,B,X,KADAPT,EPS,KMAX)
C
          WRITE(*,3000) KADAPT,EPS,KMAX,IMETH,ISWITC,OMEGA
C
          CALL ADSOR(A,4,4,B,X,KADAPT,EPS,KMAX,IMETH,ISWITC,
     *               OMEGA,WORK,RES,ITNUMB,IERR)
C
          WRITE(*,4000) OMEGA,ITNUMB,IERR
C
          WRITE(*,2040) (X(I),I=1,4)
          WRITE(*,2050) (RES(I),I=1,4)
200     CONTINUE
100   CONTINUE
      STOP
 2000 FORMAT (1X,'C[',2X,'EXAMPLE:',T78,']*',/,
     +        1X,'C[',2X,8('='),T78,']*',/,
     +        1X,'C[',2X,'COEFFICIENT MATRIX A:',T78,']*',/,
     +        1X,'C[',T78,']*')
 2010 FORMAT (1X,'C[',4(1X,D12.5),T78,']*')
 2020 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[',2X,'RIGHT HAND SIDE:',T78,']*',/,
     +        1X,'C[',T78,']*',/,
     +        1X,'C[',4(1X,D12.5),T78,']*',/,
     +        1X,'C[',T78,']*')
 2030 FORMAT (1X,'C[',2X,'STARTING VECTOR X:',T78,']*',/,
     +        1X,'C[',T78,']*',/,
     +        1X,'C[',4(1X,D12.5),T78,']*',/,
     +        1X,'C[',T78,']*',/,1X,'C[',T78,']*')
 2040 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[',2X,'SOLUTION VECTOR X:',T78,']*',/,
     +        1X,'C[',4(1X,D12.5),T78,']*')
 2050 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[',2X,'RESIDUAL VECTOR RES:',T78,']*',/,
     +        1X,'C[',4(1X,D12.5),T78,']*')
 3000 FORMAT (1X,'C[',T78,']*',/,
     *        1X,'C[',T78,']*',/,
     *        1X,'C[',2X,'INPUT PARAMETERS:',T78,']*',/,
     *        1X,'C[',4X,'KADAPT = ',I4,T78,']*',/,
     *        1X,'C[',4X,'EPS = ',D20.14,T78,']*',/,
     *        1X,'C[',4X,'KMAX = ',I4,T78,']*',/,
     *        1X,'C[',4X,'IMETH  = ',I4,T78,']*',/,
     *        1X,'C[',4X,'ISWITC = ',I4,T78,']*',/,
     *        1X,'C[',4X,'OMEGA = ',D20.14,T78,']*',/,
     *        1X,'C[',T78,']*')
 4000 FORMAT (1X,'C[',T78,']*',/,
     *        1X,'C[',2X,'OUTPUT PARAMETERS:',T78,']*',/,
     *        1X,'C[',4X,'OMEGA = ',D20.14,T78,']*',/,
     *        1X,'C[',4X,'ITNUMB = ',I4,T78,']*',/,
     *        1X,'C[',4X,'IERR = ',I4,T78,']*',/,
     *        1X,'C[',T78,']*' )
      END
C
C
      SUBROUTINE INIT (A,B,X,KADAPT,EPS,KMAX)
C
C  Initialize system matrix A and right hand side B; 
C  supply starting vector  X for the iteration
C
      DOUBLE PRECISION A(1:4,1:4), B(1:4), X(1:4),EPS
      KADAPT=4
      EPS=1.0D-09
      KMAX=100
      DO 10 I=1,4
        DO 20 J=1,4
           IF (I.EQ.J) THEN
             A(I,J)=4.0D+00
           ELSEIF ((I+J).EQ.5) THEN
             A(I,J)=0.0D+00
           ELSE
             A(I,J)=-1.0D+00
           ENDIF
20      CONTINUE
          IF (I.GT.2) THEN
            B(I)=1.0D+03
          ELSE
            B(I)=0.0D+00
          ENDIF
          X(I)=0.0D+00
10    CONTINUE
      RETURN
      END
