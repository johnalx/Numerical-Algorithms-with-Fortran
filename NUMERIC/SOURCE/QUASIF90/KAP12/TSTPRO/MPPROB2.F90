      PROGRAM TEST 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Test program for the subroutines PROB2 and APPRX2          *      
!----------------------------------------------------------------*      
!  required subroutines     : FCT, PROB2, ALPHA2, GAMMA2, NEXT2, *      
!                             APPRX2, CEPSPM, ZSPMMK, PCOSOL,    *      
!                             PCOLTG, SESSPM, SCAPRO, VECMWC,    *      
!                             ABSSUM, INDMAX, VECADD, VECXCH     *      
!----------------------------------------------------------------*      
!                                                                *      
!  In the unit square defined by (0,0),(1,1) we generate nodes   *      
!  for a square submesh using the function                       *      
!            FCT = SIN(PI*X) * COS(PI*Y).                        *      
!  We compute the coefficients of the spline function via PROB2  *      
!  for various derivative orders and use APPRX2 to interpolate   *      
!  points on the surface above the square (0.42,0.42),(0.48,0.48)*      
!  We then compare with teh exact function values.               *      
!                                                                *      
!  Results on a PC with MICROSOFT FORTRAN 5.0 compiler:                 
!                                                                *      
![                                                              ]*      
![ INTERPOLATION WITH VARIOUS DERIVATIVE ORDERS AND COMPARISON  ]*      
![ WITH EXACT FUNCTION VALUES.                                  ]*      
![                                                              ]*      
![   X    Y    F(X,Y)     2        3        4        5          ]*      
![ ======================================================       ]*      
![  .42  .42  .240877  .240635  .240520  .240668  .240842       ]*      
![  .42  .44  .181494  .181345  .181115  .181243  .181442       ]*      
![  .42  .46  .121396  .121306  .121085  .121178  .121347       ]*      
![  .42  .48  .060818  .060775  .060645  .060693  .060789       ]*      
![  .44  .42  .244285  .243691  .243782  .244020  .244254       ]*      
![  .44  .44  .184062  .183659  .183569  .183768  .184012       ]*      
![  .44  .46  .123113  .122860  .122725  .122867  .123066       ]*      
![  .44  .48  .061678  .061555  .061466  .061539  .061650       ]*      
![  .46  .42  .246729  .245840  .246122  .246422  .246701       ]*      
![  .46  .44  .185904  .185280  .185331  .185577  .185856       ]*      
![  .46  .46  .124345  .123948  .123903  .124077  .124299       ]*      
![  .46  .48  .062295  .062102  .062055  .062145  .062268       ]*      
![  .48  .42  .248199  .247117  .247531  .247866  .248173       ]*      
![  .48  .44  .187012  .186242  .186393  .186665  .186965       ]*      
![  .48  .46  .125086  .124593  .124612  .124804  .125041       ]*      
![  .48  .48  .062667  .062426  .062411  .062510  .062639       ]*      
![ ======================================================       ]*      
![ MAX. ABS. ERROR:    .001082  .000668  .000347  .000052       ]*      
!                                                                *      
!  By changing the parameter MMAX and the DO 60 loop, one can    *      
!  compute values for other derivative orders. Changing FCT is   *      
!  advised for further tests.                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author      : Hartmut Turowski                                *      
!  Date        : 6.8.1988                                        *      
!  Source code : FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!..                                                                     
!..   set maximal derivative order (MMAX)                               
!..   set up number of nodes (NNX)                                      
!..                                                                     
      PARAMETER (MMAX = 5, NN = 5, NNX = (NN + 1) * (NN + 1) ) 
!..                                                                     
!..   deklarationS                                                      
!..                                                                     
      DIMENSION TA ( (NNX + MMAX * (MMAX + 1) / 2) * (3 + NNX + MMAX *  &
      (MMAX + 1) / 2) / 2), TX (NNX), TY (NNX), TZ (NNX), TC (NNX +     &
      MMAX * (MMAX + 1) / 2), FEHLER (MMAX - 1), TERG (NNX, MMAX - 1),  &
      XX (NN * NN), YY (NN * NN), ZZ (NN * NN), WK ( (NNX + MMAX *      &
      (MMAX + 1) / 2) * ( (NNX + MMAX * (MMAX + 1) / 2) + 1) / 2)       
      INTEGER IWORK (NNX + MMAX * (MMAX + 1) / 2) 
!..                                                                     
!..   GENERATE NODES                                                    
!..                                                                     
      H = 0.2D0 
      K = 0 
      DO 20 I = 0, NN 
         DO 10 J = 0, NN 
            K = K + 1 
            TX (K) = DBLE (I) * H 
            TY (K) = DBLE (J) * H 
            TZ (K) = FCT (TX (K), TY (K) ) 
   10    END DO 
   20 END DO 
!..                                                                     
!..   number of nodes                                                   
!..                                                                     
      NX = K 
!..                                                                     
!..   initialize intermediate points for comparisons                    
!..                                                                     
      H = 0.2D-1 
      K = 0 
      DO 40 I = 1, 4 
         DO 30 J = 1, 4 
            K = K + 1 
            XX (K) = H * DBLE (I) + 0.4D0 
            YY (K) = H * DBLE (J) + 0.4D0 
            ZZ (K) = FCT (XX (K), YY (K) ) 
   30    END DO 
   40 END DO 
                                                                        
!..                                                                     
!..   loop to find apppoximate values for various                       
!..   derivative orders                                                 
!..                                                                     
      DO 60 M = 2, MMAX 
!..                                                                     
!..   compute coefficients                                              
!..                                                                     
         CALL PROB2 (NX, TX, TY, TZ, M, MARKE, TC, TA, IWORK, WK) 
         IF (MARKE.NE.1) THEN 
            WRITE ( * , '(///,1X,''MATRIX IS SINGULAR'')') 
            STOP 
         ENDIF 
!..                                                                     
!..   compare values                                                    
!..                                                                     
         FEHLER (M - 1) = 0.0D0 
         DO 50 I = 1, K 
                                                                        
            CALL APPRX2 (XX (I), YY (I), NX, M, TX, TY, TC, Z) 
!..                                                                     
!..   compute error                                                     
!..                                                                     
            FEHLER (M - 1) = DMAX1 (FEHLER (M - 1), DABS (Z - ZZ (I) ) ) 
            TERG (I, M - 1) = Z 
   50    END DO 
   60 END DO 
      WRITE ( *, 1000) (M, M = 2, MMAX) 
      DO 70 I = 1, K 
         WRITE ( *, 1100) XX (I), YY (I), ZZ (I), (TERG (I, M), M = 1,  &
         MMAX - 1)                                                      
   70 END DO 
      WRITE ( *, 1200) (FEHLER (I), I = 1, MMAX - 1) 
      STOP 
 1000 FORMAT (1X,'C[',T66,']*',/,1X,                                    &
     &        'C[ INTERPOLATION WITH VARIOUS DERIVATIVE ORDERS AND ',   &
     &        'COMPARISON',T66,']*',/,                                  &
     &        1X,'C[ WITH EXACT FUNCTION VALUES.',                      &
     &        T66,']*',/,1X,'C[',T66,']*',/,1X,'C[',                    &
     &        3X,'X',4X,'Y',4X,'F(X,Y)',1X,4(4X,I1,4X),T66,']*',/,      &
     &        1X,'C[',1X,54('='),T66,']*')                              
 1100 FORMAT (1X,'C[',1X,2(F4.2,1X),5(F8.6,1X),T66,']*') 
 1200 FORMAT (1X,'C[',1X,54('='),T66,']*',/,                            &
     &        1X,'C[',1X,'MAX. ABS. ERROR:',3X,4(F8.6,1X),T66,']*')     
      END PROGRAM TEST                              
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION FCT (X, Y) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  compute function values                                              
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION PI, X, Y 
      PI = 3.1415926D0 
      FCT = DSIN (PI * X) * DCOS (PI * Y) 
      RETURN 
      END FUNCTION FCT                              
