      PROGRAM TEST
C
C***********************************************************************
C                Test program for the subroutines                      *
C                       TSPANA und SPLFVD                              *
C----------------------------------------------------------------------*
C     required subroutines    :  ISPLTR, ISPLPE, CYTSY,  CYTSYP,       *
C                                CYTSYS, TSPANA, SPLFVD, MACHPD        *
C----------------------------------------------------------------------*
C                                                                      *
C     This program computes a transformed interpolating cubic spline.  *
C                                                                      *
C     The results are stored on TAPE3.                                 *
C                                                                      *
C[  TEST PROGRAM MPTSPANA:                                            ]*
C[  ======================                                            ]*
C[                                                                    ]*
C[  A: RESULTS AT THE NODES PHIN(I), I=0(1)12:                        ]*
C[  ------------------------------------------                        ]*
C[                                                                    ]*
C[  1) FUNCTION VALUE, 1ST AND 2ND DERIVATIVES OF THE SPLINE S(PHI):  ]*
C[                                                                    ]*
C[   I           S                 S1                 S2              ]*
C[  -----------------------------------------------------------       ]*
C[   0    .2933428029D+02    .3779969600D+01    .8510569578D+01       ]*
C[   1    .3290136775D+02    .3267260373D+01   -.9858293560D+01       ]*
C[   2    .3296209945D+02   -.2906287557D+01   -.8249423363D+01       ]*
C[   3    .3010813843D+02   -.7116818986D+01   -.6745614741D+01       ]*
C[   4    .2815137652D+02   -.6102141060D+01    .1409753282D+02       ]*
C[   5    .2754087871D+02    .4879969117D+01    .2132722453D+02       ]*
C[   6    .2977414986D+02   -.4420916643D+01   -.5334956947D+02       ]*
C[   7    .2850438563D+02    .6518158127D+01    .1183004158D+03       ]*
C[   8    .3159905062D+02    .4830853186D+00   -.1587974669D+03       ]*
C[   9    .2902585055D+02   -.8149656355D+01    .9481283603D+02       ]*
C[  10    .2933428029D+02    .3661080070D+01   -.2244158579D+02       ]*
C[  11    .2888771365D+02   -.2171521543D+01    .4337113619D+01       ]*
C[  12    .2933428029D+02    .3779969600D+01    .8510569578D+01       ]*
C[                                                                    ]*
C[                                                                    ]*
C[  2) COORDINATES OF THE CURVE K:                                    ]*
C[                                                                    ]*
C[   I          XK                 YK                                 ]*
C[  ----------------------------------------                          ]*
C[   0    .3200000000D+02    .1500000000D+02                          ]*
C[   1    .1600000000D+02    .3200000000D+02                          ]*
C[   2   -.6000000000D+01    .3100000000D+02                          ]*
C[   3   -.1800000000D+02    .1800000000D+02                          ]*
C[   4   -.2000000000D+02    .1000000000D+02                          ]*
C[   5   -.2000000000D+02   -.7000000000D+01                          ]*
C[   6   -.1300000000D+02   -.2200000000D+02                          ]*
C[   7   -.4000000000D+01   -.2600000000D+02                          ]*
C[   8    .4000000000D+01   -.3100000000D+02                          ]*
C[   9    .1200000000D+02   -.2800000000D+02                          ]*
C[  10    .2100000000D+02   -.2500000000D+02                          ]*
C[  11    .3300000000D+02   -.1100000000D+02                          ]*
C[  12    .3200000000D+02    .1500000000D+02                          ]*
C[                                                                    ]*
C[                                                                    ]*
C[  3) 1ST DERIVATIVE AND CURVATURE OF THE CURVE K:                   ]*
C[                                                                    ]*
C[   I          C1         IFEHL1          CKR        IFEHL2          ]*
C[  --------------------------------------------------------          ]*
C[   0   -.2440536782D+01    0       .2471355355D-01    0             ]*
C[   1   -.4132690074D+00    0       .3951435524D-01    0             ]*
C[   2    .5166770245D+00    0       .3795872935D-01    0             ]*
C[   3    .2445729015D+01    0       .4089191871D-01    0             ]*
C[   4    .7604006784D+01    0       .1966921345D-01    0             ]*
C[   5   -.9923080886D+01    0       .9997504304D-02    0             ]*
C[   6   -.6363028662D+00    0       .9218190030D-01    0             ]*
C[   7   -.6871592320D+00    0      -.9898459184D-01    0             ]*
C[   8   -.9476803374D-01    0       .1906304986D+00    0             ]*
C[   9    .5008953944D+00    0      -.6483734863D-01    0             ]*
C[  10    .4144120056D+00    0       .5982821411D-01    0             ]*
C[  11    .2878057104D+01    0       .2955984647D-01    0             ]*
C[  12   -.2440536782D+01    0       .2471355355D-01    0             ]*
C[                                                                    ]*
C[                                                                    ]*
C[  B: 1ST DERIVATIVE OF K IN A REGION WITH LARGE SLOPE:              ]*
C[  ----------------------------------------------------              ]*
C[                                                                    ]*
C[   I          PHI                C1            IFEHL1               ]*
C[  ---------------------------------------------------               ]*
C[   1    .2564493500D+01    .3602794001D+06       0                  ]*
C[   2    .2564494500D+01    .4142842253D+06       0                  ]*
C[   3    .2564495500D+01    .4873340500D+06       0                  ]*
C[   4    .2564496500D+01    .5916599791D+06       0                  ]*
C[   5    .2564497500D+01    .7528193508D+06       0                  ]*
C[   6    .2564498500D+01    .1034639208D+07       0                  ]*
C[   7    .2564499500D+01    .1653708045D+07       0                  ]*
C[   8    .2564500500D+01    .4117201599D+07       0                  ]*
C[   9    .2564501500D+01   -.8408008470D+07       0                  ]*
C[  10    .2564502500D+01   -.2080077659D+07       0                  ]*
C[  11    .2564503500D+01   -.1186848193D+07       0                  ]*
C[  12    .2564504500D+01   -.8303002430D+06       0                  ]*
C[  13    .2564505500D+01   -.6384885186D+06       0                  ]*
C[  14    .2564506500D+01   -.5186685533D+06       0                  ]*
C[  15    .2564507500D+01   -.4367140936D+06       0                  ]*
C[  16    .2564508500D+01   -.3771249402D+06       0                  ]*
C[  17    .2564509500D+01   -.3318451167D+06       0                  ]*
C[  18    .2564510500D+01   -.2962729178D+06       0                  ]*
C[  19    .2564511500D+01   -.2675887401D+06       0                  ]*
C[  20    .2564512500D+01   -.2439685405D+06       0                  ]*
C[                                                                    ]*
C[                                                                    ]*
C[  MEANING OF PARAMETER IFEHL1:                                      ]*
C[  ----------------------------                                      ]*
C[                                                                    ]*
C[  IFEHL1 = 0 : ALL OK                                               ]*
C[  IFEHL1 = 1 : DUE TO HN (= DENOMINATOR IN THE DEFINING EQUATION    ]*
C[               FOR C1) = 0 WE COULD NOT COMPUTE C1                  ]*
C[  IFEHL1 = 2 : DUE TO HN <= 4*MACHINE CONSTANT THE VALUE FOR C1 IS  ]*
C[               IMPRECISE                                            ]*
C                                                                      *
C***********************************************************************
C                                                                      *
C     Author      :  Gnter Palm                                       *
C     Date        :  5.20.1988                                         *
C     Source code :  FORTRAN 77                                        *
C                                                                      *
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (N=12)
C
      DIMENSION X(0:N),Y(0:N),PHIN(0:N),R(0:N)
      DIMENSION B(0:N),C(0:N),D(0:N)
      DIMENSION HILF(5*N+1)
      DIMENSION XK(0:N),YK(0:N),C1(0:N),CKR(0:N)
      INTEGER IFEHL1(0:N),IFEHL2(0:N)
C
      DATA X / 32.0D0,  16.0D0, -6.0D0,- 18.0D0, -20.0D0,
     +        -20.0D0, -13.0D0, -4.0D0,   4.0D0,  12.0D0,
     +         21.0D0,  33.0D0, 32.0D0/
      DATA Y / 15.0D0,  32.0D0,  31.0D0,  18.0D0,  10.0D0,
     +         -7.0D0, -22.0D0, -26.0D0, -31.0D0, -28.0D0,
     +        -25.0D0, -11.0D0,  15.0D0/
      DATA JV /-1/
C
      CALL ISPLTR (N,X,Y,JV,PX,PY,R,B,C,D,PHIN,PHID,HILF,IFEHL)
      IF (IFEHL .NE. 0) THEN
         WRITE (*,1000) IFEHL
         STOP
      ENDIF
      OPEN (3,FILE='TAPE3')
      REWIND (3)
      WRITE (3,2000)
      DO 100 K = 0,N,1
         CALL TSPANA (PHIN(K),N,PHIN,R,B,C,D,PHID,PX,PY,S,S1,S2,S3,
     +                XK(K),YK(K),C1(K),CKR(K),IFEHL1(K),IFEHL2(K))
         WRITE (3,2010) K, S, S1, S2
         IF (IFEHL1(K) .EQ. 1) C1(K) = 9.99D+09
  100 CONTINUE
      WRITE (3,2020)
      DO 110 K = 0,N,1
         WRITE (3,2030) K, XK(K), YK(K)
  110 CONTINUE
      WRITE (3,2040)
      DO 120 K = 0,N,1
         WRITE (3,2050) K, C1(K), IFEHL1(K), CKR(K), IFEHL2(K)
  120 CONTINUE
      WRITE (3,2060)
      PHI = 2.5644925D0
      H   = 0.000001D0
      DO 130 K = 1,20,1
         PHI = PHI + H
         CALL TSPANA (PHI,N,PHIN,R,B,C,D,PHID,PX,PY,S,S1,S2,S3,
     +                XK(0),YK(0),C1(0),CKR(0),MARKE1,MARKE2)
         IF (MARKE1 .EQ. 1) C1(0) = 9.99D+09
         WRITE (3,2070) K, PHI, C1(0), MARKE1
  130 CONTINUE
      WRITE (3,2080)
C
      WRITE (*,*)' PROGRAM RAN SUCCESSFULLY.'
      WRITE (*,*)' RESULTS ON TAPE3.'
      STOP
C
 1000 FORMAT(' ERROR WHEN COMPUTING COEFFICIENTS.',/,
     +       ' ERROR PARAMETER =',I4)
 2000 FORMAT('C[  TEST PROGRAM MPTSPANA:',T71,']*',/,
     +       'C[  ',22('='),T71,']*',/,
     +       'C[',T71,']*',/,
     +       'C[  A: RESULTS AT THE NODES ',
     +       'PHIN(I), I=0(1)12:',T71,']*',/,
     +       'C[  ',42('-'),T71,']*',/,
     +       'C[',T71,']*',/,
     +       'C[  1) FUNCTION VALUE, 1ST AND 2ND DERIVATIVES OF THE ',
     +          'SPLINE S(PHI):',T71,']*',/,
     +       'C[',T71,']*',/,
     +       'C[   I',11X,'S',17X,'S1',17X,'S2',T71,']*',/,
     +       'C[  ',59('-'),T71,']*')
 2010 FORMAT('C[  ',I2,3(3X,D16.10),T71,']*')
 2020 FORMAT('C[',T71,']*',/,
     +       'C[',T71,']*',/,
     +       'C[  2) COORDINATES OF THE CURVE K: ',T71,']*',/,
     +       'C[',T71,']*',/,
     +       'C[   I',10X,'XK',17X,'YK',T71,']*',/,
     +       'C[  ',40('-'),T71,']*')
 2030 FORMAT('C[  ',I2,2(3X,D16.10),T71,']*')
 2040 FORMAT('C[',T71,']*',/,
     +       'C[',T71,']*',/,
     +       'C[  3) 1ST DERIVATIVE AND CURVATURE OF THE CURVE K:',
     +       T71,']*',/,
     +       'C[',T71,']*',/,
     +       'C[   I',10X,'C1',9X,'IFEHL1',10X,'CKR',8X,'IFEHL2',
     +       T71,']*',/,
     +       'C[  ',56('-'),T71,']*')
 2050 FORMAT('C[  ',I2,2(3X,D16.10,4X,I1,3X),T71,']*')
 2060 FORMAT('C[',T71,']*',/,
     +       'C[',T71,']*',/,
     +       'C[  B: 1ST DERIVATIVE OF K IN A REGION WITH LARGE SLOPE:',
     +       T71,']*',/,
     +       'C[  ',52('-'),T71,']*',/,
     +       'C[',T71,']*',/,
     +       'C[   I',10X,'PHI',16X,'C1',12X,'IFEHL1',T71,']*',/,
     +       'C[  ',51('-'),T71,']*')
 2070 FORMAT('C[  ',I2,2(3X,D16.10),7X,I1,T71,']*')
 2080 FORMAT('C[',T71,']*',/,
     +       'C[',T71,']*',/,
     +       'C[  MEANING OF PARAMETER IFEHL1:',T71,']*',/,
     +       'C[  ',28('-'),T71,']*',/,
     +       'C[',T71,']*',/,
     +       'C[  IFEHL1 = 0 : ALL OK',T71,']*',/,
     +       'C[  IFEHL1 = 1 : DUE TO HN (= DENOMINATOR IN THE ',
     +       'DEFINING EQUATION',T71,']*',/,
     +       'C[',15X,'FOR C1) = 0 WE COULD NOT COMPUTE C1 ',T71,']*',/,
     +       'C[  IFEHL1 = 2 : DUE TO HN <= ',
     +       '4*MACHINE CONSTANT THE VALUE FOR C1 IS',T71,']*',/,
     +       'C[',15X,'IMPRECISE',T71,']*')
      END
