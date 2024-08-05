!                                                         4.26.87       
      PROGRAM TEST 
!                                                                       
!***********************************************************************
!                                                                       
!  Testprogram for the  SUBROUTINE GADESM.                              
!                                                                       
!  Test example: Numerische Mathematik fÅr Ingenieure;                  
!                G. Engeln-MÅllges/F. Reutter, 4th ed. 1984;            
!                Example 6.3 (III), p. 213 - 215.                       
!                                                                       
![   I *  X(I)  *  F(I)  *  W(I)  *                                     
![  -------------------------------                                     
![   0 *    .02 *  50.00 *   1.00 *                                     
![   1 *    .10 *  10.00 *   1.00 *                                     
![   2 *    .50 *   1.00 *   1.00 *                                     
![   3 *   1.00 *    .00 *   1.00 *                                     
![                                                                      
![  MAX. DEGREE OF THE LEAST SQUARES POLYNOMIAL, N =  3                 
![                                                                      
![  IFEHL =  0                                                          
![                                                                      
![   J *         C(J)         *                                         
![  ---------------------------                                         
![   0 *      62.981434240377 *                                         
![   1 *    -680.869756236137 *                                         
![   2 *    1609.739229025772 *                                         
![   3 *    -991.850907030012 *                                         
![                                                                      
![  TABLE OF VALUES                                                     
![  ---------------                                                     
![                                                                      
![    .00        62.981434240377                                        
![    .10         9.999999420817                                        
![    .20       -16.737755566964                                        
![    .30       -23.182936431181                                        
![    .40       -15.286648880048                                        
![    .50         1.000001378218                                        
![    .60        19.725908635403                                        
![    .70        34.939967183292                                        
![    .80        40.691071313669                                        
![    .90        31.028115318319                                        
![   1.00         -.000006510973                                        
![   1.10       -58.344399882423                                        
![   1.20      -149.956170504246                                        
![   1.30      -280.786424084656                                        
![   1.40      -456.786266331869                                        
![   1.50      -683.906802954101                                        
![   1.60      -968.099139659566                                        
![   1.70     -1315.314382156479                                        
![   1.80     -1731.503636153057                                        
![   1.90     -2222.618007357514                                        
![   2.00     -2794.608601478065                                        
!                                                                       
!***********************************************************************
!                                                                       
      PARAMETER (IA = 4, N = 3, M = 3) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION X (0:M), F (0:M), W (0:M), C (0:N), A (IA, N + 1),      &
      B (N + 1), Y (N + 1), Z (N + 1)                                   
      DATA X, F, W / .02D0, .1D0, .5D0, 1.D0, 50.D0, 10.D0, 1.D0, 0.D0, &
      4 * 1.D0 /                                                        
      WRITE ( *, 900) 
      DO 10 I = 0, M 
         WRITE ( *, 910) I, X (I), F (I), W (I) 
   10 END DO 
      WRITE ( * , 980) 'MAX. DEGREE OF THE LEAST SQUARES POLYNOMIAL, N =&
     & ', N                                                             
      CALL GADESM (M, X, F, W, IA, N, C, A, B, Y, Z, IFEHL) 
      WRITE ( * , 980) 'IFEHL = ', IFEHL 
      IF (IFEHL.EQ.0) THEN 
         WRITE ( *, 920) 
         DO 20 I = 0, N 
            WRITE ( *, 930) I, C (I) 
   20    END DO 
         WRITE ( *, 940) 
         DO 30 XX = 0., 2.05, .1 
            WRITE ( *, 950) XX, APPROX (N, C, XX) 
   30    END DO 
      ELSE 
         WRITE ( * , 970) '***** ERROR ***** SEE PROGRAM DESCRIPTION !' 
      ENDIF 
      STOP 
  900 FORMAT(1X,'C[',2X,' I *  X(I)  *  F(I)  *  W(I)  *',T78,']*',/,   &
     &       1X,'C[',2X,31('-'),T78,']*')                               
  910 FORMAT(1X,'C[',3X,I1,' *',3(1X,F6.2,' *'),T78,']*') 
  920 FORMAT(1X,'C[',T78,']*',/,                                        &
     &       1X,'C[',2X,' J *         C(J)         *',T78,']*',/,       &
     &       1X,'C[',2X,27('-'),T78,']*')                               
  930 FORMAT(1X,'C[',3X,I1,' * ',F20.12,' *',T78,']*') 
  940 FORMAT(1X,'C[',T78,']*',/,                                        &
     &       1X,'C[',2X,'TABLE OF VALUES',T78,']*',/,                   &
     &       1X,'C[',2X,'---------------',T78,']*',/,                   &
     &       1X,'C[',T78,']*')                                          
  950 FORMAT(1X,'C[',2X,F5.2,3X,F20.12,T78,']*') 
  970 FORMAT(1X,'C[',T78,']*',/,                                        &
     &       1X,'C[',2X,A,T78,']*')                                     
  980 FORMAT(1X,'C[',T78,']*',/,                                        &
     &       1X,'C[',2X,A,I2,T78,']*')                                  
      END PROGRAM TEST                              
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION APPROX (N, C, X) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION C (0:N) 
      APPROX = 0.D0 
      DO 10 I = N, 0, - 1 
         APPROX = APPROX * X + C (I) 
   10 END DO 
      RETURN 
      END FUNCTION APPROX                           
