      PROGRAM TESTG
C                                                                     95/05/01
C*****************************************************************************
C                                                                            *
C   This program tests the initial value problem solver routine IVP with     *
C   varying accuracy demands.                                                *
C                                                                            *
C   Test IVP for various values of k:                                        *
C       Y1' = Y2                 ,   Y1(0) =  1                              *
C       Y2' = k*Y1 + (k-1)*Y2    ,   Y2(0) = -1 ;                            *
C   exact solution  Y1(X) = EXP(-X) .                                        *
C                                                                            *
C   The output consists of the maximal error of Y1 at 10 equidistant nodes   *
C   in the interval  [0 , 5] as well in case of the Prince-Dormand formula   *
C   the result of the stiffness tests as follows:                            *
C                                                                            *
C[                                                                          ]*
C[  K           EPS        FEHL_RK2/3  FEHL_ENGL   FEHL_PR_DO               ]*
C[ ----------  ----------  ----------  ----------  ----------               ]*
C[ -5.000D+01   5.000D-05   1.273D-06   5.294D-06   2.887D-06               ]*
C[                          318         102          91  FCT-EVALUAT.       ]*
C[ -5.000D+01   5.000D-12   FEHL =  5                                       ]*
C[ -5.000D+01   5.000D-12   0.000D+00   2.575D-14   3.772D-14               ]*
C[                        50004        3420        2247  FCT-EVALUAT.       ]*
C[ -5.000D+03   5.000D-05                           FEHL = -2               ]*
C[ -5.000D+03   5.000D-05                           FEHL = -1               ]*
C[ -5.000D+03   5.000D-05                           FEHL = -1               ]*
C[ -5.000D+03   5.000D-05                           FEHL = -1               ]*
C[ -5.000D+03   5.000D-05                           FEHL = -1               ]*
C[ -5.000D+03   5.000D-05                           FEHL = -1               ]*
C[ -5.000D+03   5.000D-05                           FEHL = -1               ]*
C[ -5.000D+03   5.000D-05                           FEHL = -1               ]*
C[ -5.000D+03   5.000D-05                           FEHL = -1               ]*
C[ -5.000D+03   5.000D-05                           FEHL = -1               ]*
C[ -5.000D+03   5.000D-05   9.168D-08   2.311D-09   1.521D-08               ]*
C[                        41766       73386       58373  FCT-EVALUAT.       ]*
C[ -5.000D+03   5.000D-12   FEHL =  5                                       ]*
C[ -5.000D+03   5.000D-12                           FEHL = -1               ]*
C[ -5.000D+03   5.000D-12                           FEHL = -1               ]*
C[ -5.000D+03   5.000D-12                           FEHL = -1               ]*
C[ -5.000D+03   5.000D-12                           FEHL = -1               ]*
C[ -5.000D+03   5.000D-12                           FEHL = -1               ]*
C[ -5.000D+03   5.000D-12                           FEHL = -1               ]*
C[ -5.000D+03   5.000D-12                           FEHL = -1               ]*
C[ -5.000D+03   5.000D-12                           FEHL = -1               ]*
C[ -5.000D+03   5.000D-12                           FEHL = -1               ]*
C[ -5.000D+03   5.000D-12                           FEHL = -1               ]*
C[ -5.000D+03   5.000D-12   0.000D+00   6.770D-16   1.123D-15               ]*
C[                        50004      112686      105763  FCT-EVALUAT.       ]*
C                                                                            *
C   The  Prince-Dormand formula diagnoses stiffness (negative entry in IERR) *
C   according to the maximal value of k. Different compilers will generate   *
C   slightly different data; The magnitudes, however must be similar.        *
C                                                                            *
C*****************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (N=2)
      DIMENSION YA(N), FEHL(3), NUS(3)
      COMMON /LAMBDA/ELL
      EXTERNAL DGL
      YEX(X) = EXP(-X)
C
C *** Initialize ***  
C
      A    =  0.0D0
      B    =  5.0D0
      YANF =  1.0D0
      YANF1= -1.0D0
C
      WRITE(*,900)
      DX   = 0.5D0
C
C*****************************************************************
C     Loops with different values of k and Epsilon               *
C*****************************************************************
C
      DO 200 KAH = 0, 2, 2
      ELL = -50.0D0*10.0D0**KAH
      EA  = 5.0D2
      ER  = 5.0D2
      DO 100 IGEN = 4, 12, 8
        EA      = 1.0D-7*EA
        ER      = 1.0D-7*ER
C        
        XK      = A
        HK      = 0.1D0
        YA(1)   = YANF
        YA(2)   = YANF1
        XEND    = A
        FEHL(1) = 0.0D0
        NUS(1)  = 0
  10      XEND  = XEND + DX
          CALL IVP(XK,HK,YA,N,DGL,XEND,EA,ER,0,50000,NUSS,IFEHL)
          NUS(1) = NUS(1) + NUSS
          IF (IFEHL .NE. 0) THEN
            WRITE(*,991) ELL,ER,IFEHL
          ELSE
            FEHL(1) = MAX(FEHL(1),ABS(YA(1)-YEX(XK)))
            IF (XEND .LT. B) GOTO 10
          ENDIF
C
        XK      = A
        HK      = 0.1D0
        YA(1)   = YANF
        YA(2)   = YANF1
        XEND    = A
        FEHL(2) = 0.0D0
        NUS(2)  = 0
  20    CONTINUE
          XEND  = XEND + DX
          CALL IVP(XK,HK,YA,N,DGL,XEND,EA,ER,1,50000,NUSS,IFEHL)
          NUS(2) = NUS(2) + NUSS
          IF (IFEHL .NE. 0) THEN
            WRITE(*,992) ELL,ER,IFEHL
          ELSE
            FEHL(2) = MAX(FEHL(2),ABS(YA(1)-YEX(XK)))
            IF (XEND .LT. B) GOTO 20
          ENDIF
C
        XK      = A
        HK      = 0.1D0
        YA(1)   = YANF
        YA(2)   = YANF1
        XEND    = A
        FEHL(3) = 0.0D0
        NUS(3)  = 0
  30    CONTINUE
          XEND = XEND + DX
          CALL IVP(XK,HK,YA,N,DGL,XEND,EA,ER,-1,50000,NUSS,IFEHL)
          NUS(3) = NUS(3) + NUSS
          IF (IFEHL .NE. 0)
     #       WRITE(*,993) ELL,ER,IFEHL
          IF(IFEHL.LE.0 .AND. IFEHL.GE.-2) THEN
            FEHL(3) = MAX(FEHL(3),ABS(YA(1)-YEX(XK)))
            IF (XEND .LT. B) GOTO 30
          ENDIF
C
C
        WRITE(*,910) ELL,ER,FEHL
        WRITE(*,920) NUS
 100  CONTINUE
 200  CONTINUE
      S T O P
 900  FORMAT(1X,'C[',T78,']*',/,
     #       1X,'C[  K',11X,'EPS',8X,'FEHL_RK2/3',2X,'FEHL_ENGL',3X,
     #          'FEHL_PR_DO',T78,']*',/,
     #       1X,'C[ ',5('----------',2X),T78,']*')
 910  FORMAT(1X,'C[ ',5(1PD10.3,2X),T78,']*')
 920  FORMAT(1X,'C[ ',18X,3(I10,2X),'FCT-EVALUAT.',T78,']*')
 991  FORMAT(1X,'C[ ',2(1PD10.3,2X),' FEHL = ',I2,T78,']*')
 992  FORMAT(1X,'C[ ',2(1PD10.3,2X),12X,' FEHL = ',I2,T78,']*')
 993  FORMAT(1X,'C[ ',2(1PD10.3,2X),24X,' FEHL = ',I2,T78,']*')
      END
C
C
      SUBROUTINE DGL(X,Y,N,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(N),F(N)
      COMMON /LAMBDA/ ELL
      F(1) = Y(2)
      F(2) = ELL*Y(1) + (ELL-1.0D0)*Y(2)
      R E T U R N
      END
