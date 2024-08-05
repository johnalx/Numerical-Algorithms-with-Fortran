      PROGRAM TEST
C
C*********************************************************************
C                                                                    *
C  Testprogram for the  SUBROUTINE FRKFSY                            *
C                                                                    *
C  Test example:                                                     *
C  We solve one initial value problem using the method of            *
C  RUNGE-KUTTA-FEHLBERG.                                             *
C                                                                    *
C  Test problem:                                                     *
C               Y1'=Y2                        Y1(0)=1   I=[0,20]     *
C               Y2'=-Y1/(SQRT(Y1^2+Y3^2))^3   Y2(0)=0                *
C               Y3'=Y4                        Y3(0)=0                *
C               Y4'=-Y3/(SQRT(Y1^2+Y3^2))^3   Y4(0)=1                *
C                                                                    *
C                                                                    *
C  Results:                                                          *
C                                                                    *
C[                                                                  ]*
C[                                                                  ]*
C[ ABSERR= .00000000E+00    RELERR= .10000000E-05                   ]*
C[                                                                  ]*
C[                                                                  ]*
C[ IFEHL=  1     X=  20.000                                         ]*
C[                                                                  ]*
C[          SOLUTION:I APPROXIMATION        I EXACT                 ]*
C[          ---------+----------------------+---------------------- ]*
C[          Y1(X)    I  .40807628764144E+00 I  .40808206181339E+00  ]*
C[          Y2(X)    I -.91294791088434E+00 I -.91294525072763E+00  ]*
C[          Y3(X)    I  .91294748136748E+00 I  .91294525072763E+00  ]*
C[          Y4(X)    I  .40807650152524E+00 I  .40808206181339E+00  ]*
C                                                                    *
C*********************************************************************
C
C declarations
C
      EXTERNAL DGL
      DOUBLE PRECISION Y(4),YEX(4)
      DOUBLE PRECISION X,BETA,H,HMX,ABSERR,RELERR
C
C initialize
C
      ABSERR=0.0D+00
      RELERR=1.0D-06
      WRITE(*,2000) ABSERR,RELERR
      N=4
      X=0.0D+00
      H=0.01D+00
      HMX=0.1D+00
      BETA=20.D+00
      Y(1)=1.0D+00
      Y(2)=0.0D+00
      Y(3)=0.0D+00
      Y(4)=1.0D+00
C
C begin computations
C
      CALL FRKFSY(X,BETA,N,Y,DGL,H,HMX,ABSERR,RELERR,IFEHL)
      CALL EXAKT(X,N,YEX)
      WRITE(*,2100) IFEHL,X
C
C put out results
C
      DO 40 I=1,N
         WRITE(*,2200) I,Y(I),YEX(I)
40    CONTINUE
      STOP
C
C Format declarations
C
2000  FORMAT(1X,'C[',T70,']*',/,1X,'C[',T70,']*',/,
     +       1X,'C[',1X,'ABSERR=',E14.8,'    RELERR=',E14.8,T70,']*',/,
     +       1X,'C[',T70,']*')
2100  FORMAT(1X,'C[',T70,']*',/,
     +       1X,'C[',1X,'IFEHL=',I3,5X,'X=',F8.3,T70,']*',/,
     +       1X,'C[',T70,']*',/,
     +       1X,'C[',10X,'SOLUTION:I APPROXIMATION',8X,'I EXACT',
     +          T70,']*',/,
     +       1X,'C[',10X,9('-'),'+',22('-'),'+',22('-'),T70,']*')
2200  FORMAT(1X,'C[',10X,'Y',I1,'(X)',4X,'I ',E20.14,' I ',E20.14,
     +          T70,']*')
      END
C
C Test problem
C
      SUBROUTINE DGL(X,Y,F)
      DOUBLE PRECISION Y(4),F(4),X
      F(1)=Y(2)
      F(2)=-Y(1)/(SQRT(Y(1)*Y(1)+Y(3)*Y(3)))**3.D+00
      F(3)=Y(4)
      F(4)=-Y(3)/(SQRT(Y(1)*Y(1)+Y(3)*Y(3)))**3.D+00
      RETURN
      END
C
      SUBROUTINE EXAKT(X,N,YEXAKT)
      DOUBLE PRECISION X,YEXAKT(N)
      YEXAKT(1)=COS(X)
      YEXAKT(2)=-SIN(X)
      YEXAKT(3)=SIN(X)
      YEXAKT(4)=COS(X)
      RETURN
      END
