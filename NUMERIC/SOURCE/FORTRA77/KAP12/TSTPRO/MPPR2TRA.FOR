      PROGRAM TEST
C
C*****************************************************************
C                                                                *
C     Test program for the subroutines PR2TRA and APPRT2         *
C----------------------------------------------------------------*
C  required subroutiens     : FCT,PR2TRA,ALPHA2,GAMMA2,NEXT2,    *
C                             APPRT2,CEPSPM,ZSPMMK,PCOSOL,       *
C                             PCOLTG,SESSPM,SCAPRO,VECMWC,       *
C                             ABSSUM,INDMAX,VECADD,VECXCH        *
C----------------------------------------------------------------*
C                                                                *
C  In the unit square defined by (0,0),(1,1) we generate nodes   *
C  for a square submesh using the function                       *  
C            FCT = SIN(PI*X) * COS(PI*Y).                        *
C  We compute the coefficients of the spline function via PR2TRA *
C  for various derivative orders and use APPRT2 to interpolate   *
C  points on the surface above the square (0.42,0.42),(0.48,0.48)*
C  We then compare with teh exact function values.               *
C                                                                *
C  Results on a PC with MICROSOFT FORTRAN 5.0 compiler:          
C                                                                *
C[                                                              ]*
C[ INTERPOLATION WITH VARIOUS DERIVATIVE ORDERS AND COMPARISON  ]*
C[ WITH EXACT FUNCTION VALUES.                                  ]*
C[                                                              ]*
C[   X    Y    F(X,Y)     2        3        4        5          ]*
C[ ======================================================       ]*
C[  .42  .42  .240877  .240635  .240520  .240668  .240842       ]*
C[  .42  .44  .181494  .181345  .181115  .181243  .181442       ]*
C[  .42  .46  .121396  .121306  .121085  .121178  .121347       ]*
C[  .42  .48  .060818  .060775  .060645  .060693  .060789       ]*
C[  .44  .42  .244285  .243691  .243782  .244020  .244254       ]*
C[  .44  .44  .184062  .183659  .183569  .183768  .184012       ]*
C[  .44  .46  .123113  .122860  .122725  .122867  .123066       ]*
C[  .44  .48  .061678  .061555  .061466  .061539  .061650       ]*
C[  .46  .42  .246729  .245840  .246122  .246422  .246701       ]*
C[  .46  .44  .185904  .185280  .185331  .185577  .185856       ]*
C[  .46  .46  .124345  .123948  .123903  .124077  .124299       ]*
C[  .46  .48  .062295  .062102  .062055  .062145  .062268       ]*
C[  .48  .42  .248199  .247117  .247531  .247866  .248173       ]*
C[  .48  .44  .187012  .186242  .186393  .186665  .186965       ]*
C[  .48  .46  .125086  .124593  .124612  .124804  .125041       ]*
C[  .48  .48  .062667  .062426  .062411  .062510  .062639       ]*
C[ ======================================================       ]*
C[ MAX. ABS. ERROR:    .001082  .000668  .000347  .000052       ]*
C                                                                *
C  By changing the parameter MMAX and the DO 60 loop, one can    *
C  compute values for other derivative orders. Changing FCT is   *
C  advised for further tests.                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  Author      : Hartmut Turowski                                *
C  Date        : 6.8.1988                                        *
C  Source code : FORTRAN 77                                      *
C                                                                *
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C..
C..   set maximal derivative order (MMAX)
C..   set up number of nodes (NNX)
C..
      PARAMETER (MMAX=5,NN=5,NNX=(NN+1)*(NN+1))
C..
C..   declarations
C..
      DIMENSION TA((NNX+MMAX*(MMAX+1)/2)*
     F             (3+NNX+MMAX*(MMAX+1)/2)/2),
     F          TX(NNX), TY(NNX), TZ(NNX),
     F          TC(NNX+MMAX*(MMAX+1)/2), FEHLER(MMAX-1),
     F          TERG(NNX,MMAX-1),XX(NN*NN),YY(NN*NN),ZZ(NN*NN),
     F          WK((NNX+MMAX*(MMAX+1)/2)*
     F            ((NNX+MMAX*(MMAX+1)/2)+1)/2)
      INTEGER IWORK (NNX+MMAX*(MMAX+1)/2)

C..
C..   loop to compute values for various derivative orders
C..
      DO 60 M = 2,MMAX
C..
C..   create nodes
C..
      H = 0.2D0
      K = 0
      DO 20 I= 0,NN
         DO 10 J= 0,NN
            K = K + 1
            TX(K)= DBLE(I) * H
            TY(K)= DBLE(J) * H
            TZ(K)= FCT (TX(K),TY(K))
10       CONTINUE
20    CONTINUE
C..
C..   number of nodes
C..
      NX = K
C..
C..   initialize intermediate values, where we shall compare 
C..   the result
C..
      H = 0.2D-1
      K = 0
      DO 40 I= 1, 4
         DO 30 J= 1, 4
            K = K + 1
            XX(K)= H * DBLE(I) + 0.4D0
            YY(K)= H * DBLE(J) + 0.4D0
            ZZ(K)= FCT (XX(K),YY(K))
30       CONTINUE
40    CONTINUE
C..
C..   compute coefficients
C..
          CALL PR2TRA (NX,TX,TY,TZ,M,MARKE,TC,TA,IWORK,WK
     F                                     ,XQUER,YQUER,R)
          IF (MARKE .NE. 1) THEN
             WRITE (*,'(///,1X,''MATRIX IS SINGULAR'')')
             STOP
          END IF
C..
C..   compare results
C..
          FEHLER (M-1) = 0.0D0
          DO 50 I=1,K
              XHILF = XX(I)
              YHILF = YY(I)

              CALL APPRT2(XHILF,YHILF,NX,M,TX,TY,TC,Z,
     F                                  XQUER,YQUER,R)
C..
C..   compute error
C..
              FEHLER (M-1) = DMAX1 (FEHLER(M-1),DABS(Z-ZZ(I)))
              TERG(I,M-1) = Z
50        CONTINUE
60    CONTINUE
      WRITE (*,1000)(M,M=2,MMAX)
      DO 70 I = 1,K
         WRITE (*,1100) XX(I),YY(I),ZZ(I),
     F                  (TERG (I,M),M = 1, MMAX-1)
70    CONTINUE
      WRITE (*,1200) (FEHLER (I),I = 1, MMAX-1)
      STOP
1000  FORMAT (1X,'C[',T66,']*',/,1X,
     +       'C[ INTERPOLATION WITH VARIOUS DERIVATIVE ORDERS AND ',
     +       'COMPARISON',T66,']*',/,
     +       1X,'C[ WITH EXACT FUNCTION VALUES.',
     +       T66,']*',/,1X,'C[',T66,']*',/,1X,'C[',
     +       3X,'X',4X,'Y',4X,'F(X,Y)',1X,4(4X,I1,4X),T66,']*',/,
     +       1X,'C[',1X,54('='),T66,']*')
1100  FORMAT (1X,'C[',1X,2(F4.2,1X),5(F8.6,1X),T66,']*')
1200  FORMAT (1X,'C[',1X,54('='),T66,']*',/,
     +        1X,'C[',1X,'MAX. ABS. ERROR:',3X,4(F8.6,1X),T66,']*')
      END
C
C
      DOUBLE PRECISION FUNCTION FCT (X,Y)
C
C*****************************************************************
C                                                                *
C  COMPUTE FUNCTION VALUES
C                                                                *
C*****************************************************************
C
      DOUBLE PRECISION PI,X,Y
      PI = 3.1415926D0
      FCT = DSIN(PI*X)*DCOS(PI*Y)
      RETURN
      END
