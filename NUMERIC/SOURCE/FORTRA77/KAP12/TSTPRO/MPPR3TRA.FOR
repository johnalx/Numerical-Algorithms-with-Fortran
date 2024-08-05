      PROGRAM TEST
C
C*****************************************************************
C                                                                *
C  Test program for the subroutines PR3TRA and APPTR3            *
C----------------------------------------------------------------*
C  required subroutines     : FCT,PROB3,ALPHA3,GAMMA3,NEXT3,     *
C                             APPTR3,CEPSPM,ZSPMMK,PCOSOL,       *
C                             PCOLTG,SESSPM,SCAPRO,VECMWC,       *
C                             ABSSUM,INDMAX,VECADD,VECXCH        *
C----------------------------------------------------------------*
C                                                                *
C  For the unit cube defined by (0,0,0),(1,1,1) we create a three*
C  dimensional mesh of nodes using the function                  *
C        FUNCTION FCT = SIN(PI*X) * SIN(PI*Y) * COS(PI*Z).       *
C  The nodes are transformed onto the unit sphere which may add  *
C  to the accuracy of the computations.                          *
C  We find the spline coefficients for various derivative orders *
C  using PROB3 in the subcube defined by (0.42,0.42,0.42) and    *
C  (0.44,0.44,0.44). Then we interpolate points on the surface   *
C  using APPTR3 and compare with the exact function values.      *
C                                                                *
C  Results on a PC with  MICROSOFT FORTRAN compiler V5.0 :       *
C                                                                *
C[                                                              ]*
C[ INTERPOLATION WITH VARIOUS DERIVATIVE ORDERS AND COMPARISON  ]*
C[ WITH EXACT FUNCTION VALUES.                                  ]*
C[                                                              ]*
C[   X    Y    Z    F(X,Y)     2        3        4              ]*
C[ ==================================================           ]*
C[  .42  .42  .42  .233309  .223174  .235782  .236707           ]*
C[  .42  .42  .44  .175792  .190602  .202108  .202854           ]*
C[  .42  .44  .42  .236610  .262484  .277021  .277949           ]*
C[  .42  .44  .44  .178280  .198618  .210185  .211021           ]*
C[  .44  .42  .42  .236610  .258160  .272428  .273374           ]*
C[  .44  .42  .44  .178280  .197642  .209140  .209978           ]*
C[  .44  .44  .42  .239958  .263132  .276995  .278099           ]*
C[  .44  .44  .44  .180802  .197101  .208042  .208999           ]*
C[ ==================================================           ]*
C[ MAX. ABS. ERROR:         .025873  .040411  .041338           ]*
C                                                                *
C  By changing the parameter MMAX and the DO 60 loop, one can    *
C  compute values for other derivative orders. Changing FCT is   *
C  advised for further tests.                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  Author      : Hartmut Turowski                                *
C  Date        : 12.8.1989                                       *
C  Source code : FORTRAN 77                                      *
C                                                                *
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C..
C..   SET MAXIMAL DERIVATIVE ORDER (MMAX)
C..   and number of nodes (NNX)
C..
      PARAMETER (MMAX=4,NN=4,NNX=(NN+1)*(NN+1)*(NN+1))
C..
C..   deklarations
C..
      DIMENSION TA((NNX+MMAX*(MMAX+1)*(MMAX+2)/6)*
     F           (3+NNX+MMAX*(MMAX+1)*(MMAX+2)/6)/2),
     F           TX(NNX), TY(NNX), TZ(NNX), TF(NNX),
     F           HX(NNX), HY(NNX), HZ(NNX),
     F           TC(NNX+MMAX*(MMAX+1)*(MMAX+2)/6), FEHLER(MMAX-1),
     F           TERG(NNX,MMAX-1),XX(NNX),YY(NNX),ZZ(NNX),
     F           FF(NNX), WK((NNX+MMAX*(MMAX+1)*(MMAX+2)/6)*
     F                       (NNX+MMAX*(MMAX+1)*(MMAX+2)/6+1)/2)
      DIMENSION IWORK (NNX+MMAX*(MMAX+1)*(MMAX+2)/6)
C..
C..   create nodes
C..
      H = 0.25D0
      L = 0
      DO 10 I= 0,NN
         DO 10 J= 0,NN
            DO 10 K= 0,NN
               L = L + 1
               TX(L)= DBLE(I) * H
               TY(L)= DBLE(J) * H
               TZ(L)= DBLE(K) * H
               TF(L)= FCT (TX(L),TY(L),TZ(L))
10    CONTINUE
C..
C..   number of nodes
C..
      NX = L
C..
C..   set up intermediate values for comparisons
C..
      H = 0.02D0
      L = 0
      DO 20 I= 1, 2
         DO 20 J= 1, 2
            DO 20 K= 1, 2
               L = L + 1
               XX(L)= H * DBLE(I) + 0.4D0
               YY(L)= H * DBLE(J) + 0.4D0
               ZZ(L)= H * DBLE(K) + 0.4D0
               FF(L)= FCT (XX(L),YY(L),ZZ(L))
20    CONTINUE
C..
C..   loop to find intermediate function values approximately
C..   for various derivative orders
C..
      DO 60 M = 2,MMAX
C..
C..   compute coefficients
C..
          DO 30 I=1,NX
             HX (I) = TX (I)
             HY (I) = TY (I)
             HZ (I) = TZ (I)
30        CONTINUE
          CALL PR3TRA (NX,HX,HY,HZ,TF,M,MARKE,TC,TA,IWORK,WK,
     F                                   XQUER,YQUER,ZQUER,R)
          IF (MARKE .NE. 1) THEN
             WRITE (*,'(///,1X,''MATRIX IS SINGULAR'')')
             STOP
          END IF
C..
C..   compare values
C..
          FEHLER (M-1) = 0.0
          DO 50 I=1,L
              XXX = XX(I)
              YYY = YY(I)
              ZZZ = ZZ(I)
              F = APPTR3(XXX,YYY,ZZZ,NX,M,HX,HY,HZ,TC,
     F                           XQUER,YQUER,ZQUER,R )
C..
C..   compute error
C..
              FEHLER (M-1) = DMAX1 (FEHLER(M-1),DABS(F-FF(I)))
              TERG(I,M-1) = F
50        CONTINUE
60    CONTINUE
      WRITE (*,1000)(M,M=2,MMAX)
      DO 70 I = 1,L
         WRITE (*,1100) XX(I),YY(I),ZZ(I),FF(I),
     F                  (TERG (I,M),M = 1, MMAX-1)
70    CONTINUE
      WRITE (*,1200) (FEHLER (I),I = 1, MMAX-1)
      STOP
1000  FORMAT (1X,'C[',T66,']*',/,1X,
     +       'C[ INTERPOLATION WITH VARIOUS DERIVATIVE ORDERS AND ',
     +       'COMPARISON',T66,']*',/,
     +       1X,'C[ WITH EXACT FUNCTION VALUES.',
     +       T66,']*',/,1X,'C[',T66,']*',/,1X,'C[',
     +       3X,'X',4X,'Y',4X,'Z',4X,'F(X,Y)',1X,3(4X,I1,4X),T66,']*',/,
     +       1X,'C[',1X,50('='),T66,']*')
1100  FORMAT (1X,'C[',1X,3(F4.2,1X),4(F8.6,1X),T66,']*')
1200  FORMAT (1X,'C[',1X,50('='),T66,']*',/,
     +        1X,'C[',1X,'MAX. ABS. ERROR:',8X,3(F8.6,1X),T66,']*')
      END
C
C
      DOUBLE PRECISION FUNCTION FCT (X,Y,Z)
C
C*****************************************************************
C                                                                *
C  compute functional values
C                                                                *
C*****************************************************************
C
      DOUBLE PRECISION PI,X,Y,Z
      PI = 3.1415926D0
      FCT = DSIN(PI*X)*DSIN(PI*Y)*DCOS(PI*Z)
      RETURN
      END
