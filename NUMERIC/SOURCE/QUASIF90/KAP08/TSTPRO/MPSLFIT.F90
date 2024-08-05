!                                                         89.01.14      
      PROGRAM TEST 
!                                                                       
!***********************************************************************
!                                                                       
!  Testprogram for the SUBROUTINES SLFIT, SLPRE.                        
!                                                                       
!  Test example: Numerische Mathematik fÅr Ingenieure;                  
!                G. Engeln-MÅllges/F. Reutter, 5th ed. 1987;            
!                EXAMPLE 6.3 (II), p. 213 - 215.                        
!                                                                       
!  Results from a PC:                                                   
!                                                                       
![   I *  X(I)  *  F(I)  *  W(I)  *                                     
![  -------------------------------                                     
![   0 *    .02 *  50.00 *   1.00 *                                     
![   1 *    .10 *  10.00 *   1.00 *                                     
![   2 *    .50 *   1.00 *   1.00 *                                     
![   3 *   1.00 *    .00 *   1.00 *                                     
![                                                                      
![  NUMBER OF FUNCTIONS, N+1 =  3                                       
![                                                                      
![  IFEHL =  0                                                          
![                                                                      
![   J *         C(J)         *                                         
![  ---------------------------                                         
![   0 *      39.678891125667 *                                         
![   1 *    -136.551407186668 *                                         
![   2 *      97.982953935049 *                                         
![                                                                      
![  TABLE OF VALUES:                                                    
![  ----------------                                                    
![                                                                      
![    .00        39.678891125667                                        
![    .10        27.003579946350                                        
![    .20        16.287927845735                                        
![    .30         7.531934823821                                        
![    .40          .735600880607                                        
![    .50        -4.101073983905                                        
![    .60        -6.978089769717                                        
![    .70        -7.895446476827                                        
![    .80        -6.853144105237                                        
![    .90        -3.851182654945                                        
![   1.00         1.110437874047                                        
![   1.10         8.031717481741                                        
![   1.20        16.912656168135                                        
![   1.30        27.753253933231                                        
![   1.40        40.553510777027                                        
![   1.50        55.313426699524                                        
![   1.60        72.033001700722                                        
![   1.70        90.712235780622                                        
![   1.80       111.351128939222                                        
![   1.90       133.949681176523                                        
![   2.00       158.507892492525                                        
![                                                                      
![  AVERAGE LEAST SQUARES ERROR =  .22038885324468E+02                  
!                                                                       
!***********************************************************************
!                                                                       
      PARAMETER (IA = 3, N = 2, M = 3, IANZ = 21) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION X (0:M), Y (0:M), W (0:M), A (0:IA, 0:N + 1), D (0:N),  &
      C (0:N), FW (1:IANZ), PSI (1:IANZ)                                
      EXTERNAL FKT 
      DATA X, Y, W / .02D0, .1D0, .5D0, 1.D0, 50.D0, 10.D0, 1.D0, 0.D0, &
      4 * 1.D0 /                                                        
      WRITE ( *, 900) 
      DO 10 I = 0, M 
         WRITE ( *, 910) I, X (I), Y (I), W (I) 
   10 END DO 
      WRITE ( * , 980) 'NUMBER OF FUNCTIONS, N+1 = ', N + 1 
      XX = 0. 
      DO 20 I = 1, IANZ 
         PSI (I) = XX 
         XX = XX + .1D0 
   20 END DO 
      CALL SLFIT (X, Y, W, 0, FKT, PSI, IA, M, N, IANZ, A, D, C, FW,    &
      QUADFE, IFEHL)                                                    
      WRITE ( * , 980) 'IFEHL = ', IFEHL 
      IF (IFEHL.EQ.0) THEN 
         WRITE ( *, 920) 
         DO 30 I = 0, N 
            WRITE ( *, 930) I, C (I) 
   30    END DO 
         WRITE ( *, 940) 
         DO 40 I = 1, IANZ 
            WRITE ( *, 950) PSI (I), FW (I) 
   40    END DO 
         WRITE ( * , 960) 'AVERAGE LEAST SQUARES ERROR = ', QUADFE 
      ELSE 
         WRITE ( * , 970) '***** ERROR ***** SEE PROGRAM DESCRIPTION!' 
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
     &       1X,'C[',2X,'TABLE OF VALUES:',T78,']*',/,                  &
     &       1X,'C[',2X,'----------------',T78,']*',/,                  &
     &       1X,'C[',T78,']*')                                          
  950 FORMAT(1X,'C[',2X,F5.2,3X,F20.12,T78,']*') 
  960 FORMAT(1X,'C[',T78,']*',/,                                        &
     &       1X,'C[',2X,A,E20.14,T78,']*')                              
  970 FORMAT(1X,'C[',T78,']*',/,                                        &
     &       1X,'C[',2X,A,T78,']*')                                     
  980 FORMAT(1X,'C[',T78,']*',/,                                        &
     &       1X,'C[',2X,A,I2,T78,']*')                                  
      END PROGRAM TEST                              
!                                                                       
!                                                                       
      SUBROUTINE FKT (X, N, F) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION F (0:N) 
      F (0) = 1.D0 
      F (1) = X 
      F (2) = X * X 
      RETURN 
      END SUBROUTINE FKT                            
