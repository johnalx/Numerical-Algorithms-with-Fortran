      PROGRAM TEST
C      
C*****************************************************************
C                                                                *
C  This program tests the subroutine GEAR4 on the stiff IVP      *
C        Y1' = Y2                           , Y1(0) =  1         *
C        Y2' = DK*Y(1) + (DK - 1)*Y(2)      , Y2(0) = -1         *
C  with DK = -10000 . The exact solution of the problem is       *
C        Y1(X) = EXP(X)  and  Y2(X) = -EXP(-X) .                 *
C  We test this problem with varying values for  EPS... .        *
C                                                                *
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(2), EPS(2,4)
      EXTERNAL DGL
      YEX(X) = EXP(-X)
      DATA EPS / 1.D-6,1.D-6,1.D-10,1.D-10,1.D-11,0.D0,0.D0,1.D-11 /
      DO 50 K = 1, 4
        X    =  0.D0
        Y(1) =  1.D0
        Y(2) = -1.D0
        H    =  1.D-2
        WRITE(*,*) 'EPSABS = ',EPS(1,K), '  EPSREL = ',EPS(2,K)
        DO 10 I = 5, 50, 5
          XX = 1.D-1*DBLE(I)
          CALL GEAR4(X, H, Y, 2, DGL, XX, EPS(1,K), EPS(2,K),
     #               50000, NN, IFEHL)
          WRITE(*,*) 'ERRORPARAMETER:',IFEHL,'  NUSED:',NN
          WRITE(*,*) '    X:',X,'  Y:',Y(1)
          WRITE(*,*) 'ERROR:',ABS(Y(1)-YEX(X))
  10    CONTINUE
        WRITE(*,*)
  50  CONTINUE
      STOP
      END
C      
      SUBROUTINE DGL(X, Y, N, F)
      DOUBLE PRECISION Y(1:N), F(1:N), X, DK
      PARAMETER( DK = -1.D4 )
      F(1) = Y(2)
      F(2) = DK*Y(1) + (DK - 1.D0)*Y(2)
      RETURN
      END
