      PROGRAM TEST 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Test program for the subroutines PR3TRA and APPTR3            *      
!----------------------------------------------------------------*      
!  required subroutines     : FCT,PROB3,ALPHA3,GAMMA3,NEXT3,     *      
!                             APPTR3,CEPSPM,ZSPMMK,PCOSOL,       *      
!                             PCOLTG,SESSPM,SCAPRO,VECMWC,       *      
!                             ABSSUM,INDMAX,VECADD,VECXCH        *      
!----------------------------------------------------------------*      
!                                                                *      
!  For the unit cube defined by (0,0,0),(1,1,1) we create a three*      
!  dimensional mesh of nodes using the function                  *      
!        FUNCTION FCT = SIN(PI*X) * SIN(PI*Y) * COS(PI*Z).       *      
!  The nodes are transformed onto the unit sphere which may add  *      
!  to the accuracy of the computations.                          *      
!  We find the spline coefficients for various derivative orders *      
!  using PROB3 in the subcube defined by (0.42,0.42,0.42) and    *      
!  (0.44,0.44,0.44). Then we interpolate points on the surface   *      
!  using APPTR3 and compare with the exact function values.      *      
!                                                                *      
!  Results on a PC with  MICROSOFT FORTRAN compiler V5.0 :       *      
!                                                                *      
![                                                              ]*      
![ INTERPOLATION WITH VARIOUS DERIVATIVE ORDERS AND COMPARISON  ]*      
![ WITH EXACT FUNCTION VALUES.                                  ]*      
![                                                              ]*      
![   X    Y    Z    F(X,Y)     2        3        4              ]*      
![ ==================================================           ]*      
![  .42  .42  .42  .233309  .223174  .235782  .236707           ]*      
![  .42  .42  .44  .175792  .190602  .202108  .202854           ]*      
![  .42  .44  .42  .236610  .262484  .277021  .277949           ]*      
![  .42  .44  .44  .178280  .198618  .210185  .211021           ]*      
![  .44  .42  .42  .236610  .258160  .272428  .273374           ]*      
![  .44  .42  .44  .178280  .197642  .209140  .209978           ]*      
![  .44  .44  .42  .239958  .263132  .276995  .278099           ]*      
![  .44  .44  .44  .180802  .197101  .208042  .208999           ]*      
![ ==================================================           ]*      
![ MAX. ABS. ERROR:         .025873  .040411  .041338           ]*      
!                                                                *      
!  By changing the parameter MMAX and the DO 60 loop, one can    *      
!  compute values for other derivative orders. Changing FCT is   *      
!  advised for further tests.                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author      : Hartmut Turowski                                *      
!  Date        : 12.8.1989                                       *      
!  Source code : FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!..                                                                     
!..   SET MAXIMAL DERIVATIVE ORDER (MMAX)                               
!..   and number of nodes (NNX)                                         
!..                                                                     
      PARAMETER (MMAX = 4, NN = 4, NNX = (NN + 1) * (NN + 1) * (NN + 1) &
      )                                                                 
!..                                                                     
!..   deklarations                                                      
!..                                                                     
      DIMENSION TA ( (NNX + MMAX * (MMAX + 1) * (MMAX + 2) / 6) *       &
      (3 + NNX + MMAX * (MMAX + 1) * (MMAX + 2) / 6) / 2), TX (NNX),    &
      TY (NNX), TZ (NNX), TF (NNX), HX (NNX), HY (NNX), HZ (NNX),       &
      TC (NNX + MMAX * (MMAX + 1) * (MMAX + 2) / 6), FEHLER (MMAX - 1), &
      TERG (NNX, MMAX - 1), XX (NNX), YY (NNX), ZZ (NNX), FF (NNX),     &
      WK ( (NNX + MMAX * (MMAX + 1) * (MMAX + 2) / 6) * (NNX + MMAX *   &
      (MMAX + 1) * (MMAX + 2) / 6 + 1) / 2)                             
      DIMENSION IWORK (NNX + MMAX * (MMAX + 1) * (MMAX + 2) / 6) 
!..                                                                     
!..   create nodes                                                      
!..                                                                     
      H = 0.25D0 
      L = 0 
      DO 10 I = 0, NN 
         DO 10 J = 0, NN 
            DO 10 K = 0, NN 
               L = L + 1 
               TX (L) = DBLE (I) * H 
               TY (L) = DBLE (J) * H 
               TZ (L) = DBLE (K) * H 
               TF (L) = FCT (TX (L), TY (L), TZ (L) ) 
   10 CONTINUE 
!..                                                                     
!..   number of nodes                                                   
!..                                                                     
      NX = L 
!..                                                                     
!..   set up intermediate values for comparisons                        
!..                                                                     
      H = 0.02D0 
      L = 0 
      DO 20 I = 1, 2 
         DO 20 J = 1, 2 
            DO 20 K = 1, 2 
               L = L + 1 
               XX (L) = H * DBLE (I) + 0.4D0 
               YY (L) = H * DBLE (J) + 0.4D0 
               ZZ (L) = H * DBLE (K) + 0.4D0 
               FF (L) = FCT (XX (L), YY (L), ZZ (L) ) 
   20 CONTINUE 
!..                                                                     
!..   loop to find intermediate function values approximately           
!..   for various derivative orders                                     
!..                                                                     
      DO 60 M = 2, MMAX 
!..                                                                     
!..   compute coefficients                                              
!..                                                                     
         DO 30 I = 1, NX 
            HX (I) = TX (I) 
            HY (I) = TY (I) 
            HZ (I) = TZ (I) 
   30    END DO 
         CALL PR3TRA (NX, HX, HY, HZ, TF, M, MARKE, TC, TA, IWORK, WK,  &
         XQUER, YQUER, ZQUER, R)                                        
         IF (MARKE.NE.1) THEN 
            WRITE ( * , '(///,1X,''MATRIX IS SINGULAR'')') 
            STOP 
         ENDIF 
!..                                                                     
!..   compare values                                                    
!..                                                                     
         FEHLER (M - 1) = 0.0 
         DO 50 I = 1, L 
            XXX = XX (I) 
            YYY = YY (I) 
            ZZZ = ZZ (I) 
            F = APPTR3 (XXX, YYY, ZZZ, NX, M, HX, HY, HZ, TC, XQUER,    &
            YQUER, ZQUER, R)                                            
!..                                                                     
!..   compute error                                                     
!..                                                                     
            FEHLER (M - 1) = DMAX1 (FEHLER (M - 1), DABS (F - FF (I) ) ) 
            TERG (I, M - 1) = F 
   50    END DO 
   60 END DO 
      WRITE ( *, 1000) (M, M = 2, MMAX) 
      DO 70 I = 1, L 
         WRITE ( *, 1100) XX (I), YY (I), ZZ (I), FF (I), (TERG (I, M), &
         M = 1, MMAX - 1)                                               
   70 END DO 
      WRITE ( *, 1200) (FEHLER (I), I = 1, MMAX - 1) 
      STOP 
 1000 FORMAT (1X,'C[',T66,']*',/,1X,                                    &
     &       'C[ INTERPOLATION WITH VARIOUS DERIVATIVE ORDERS AND ',    &
     &       'COMPARISON',T66,']*',/,                                   &
     &       1X,'C[ WITH EXACT FUNCTION VALUES.',                       &
     &       T66,']*',/,1X,'C[',T66,']*',/,1X,'C[',                     &
     &       3X,'X',4X,'Y',4X,'Z',4X,'F(X,Y)',1X,3(4X,I1,4X),T66,']*',/,&
     &       1X,'C[',1X,50('='),T66,']*')                               
 1100 FORMAT (1X,'C[',1X,3(F4.2,1X),4(F8.6,1X),T66,']*') 
 1200 FORMAT (1X,'C[',1X,50('='),T66,']*',/,                            &
     &        1X,'C[',1X,'MAX. ABS. ERROR:',8X,3(F8.6,1X),T66,']*')     
      END PROGRAM TEST                              
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION FCT (X, Y, Z) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  compute functional values                                            
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION PI, X, Y, Z 
      PI = 3.1415926D0 
      FCT = DSIN (PI * X) * DSIN (PI * Y) * DCOS (PI * Z) 
      RETURN 
      END FUNCTION FCT                              
