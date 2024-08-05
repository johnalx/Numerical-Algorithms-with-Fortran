      PROGRAM TEST 
!                                                ( Thomas Meuser )      
!***********************************************************************
!                                                                       
!     Testprogram for the subroutines  TRDIG, TRDIGP, TRDIGS.           
!     Solve a linear system   A*X = RS  with strongly nonsingular tridia
!     system matrix A.                                                  
!                                                                       
!     The test example gives the following output:                      
!                                                                       
![  EXAMPLE:                                                            
![  ========                                                            
![  COEFFICIENT MATRIX A:                                               
![                                                                      
![   .50000D+01  -.10000D+01   .00000D+00   .00000D+00                  
![  -.20000D+01   .50000D+01  -.10000D+01   .00000D+00                  
![   .00000D+00  -.20000D+01   .30000D+01  -.10000D+01                  
![   .00000D+00   .00000D+00  -.10000D+01   .50000D+01                  
![                                                                      
![  RIGHT HAND SIDE                                                     
![                                                                      
![   .40000D+01   .20000D+01   .00000D+00   .40000D+01                  
![                                                                      
![                                                                      
![  SOLUTION                                                            
![                                                                      
![   .10000D+01   .10000D+01   .10000D+01   .10000D+01                  
![  STOP. NO ERROR!                                                     
!                                                                       
!     Other tests possible for different data.                          
!                                                                       
!***********************************************************************
!                                                                       
      PARAMETER (N = 4, M = N * N) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION UD (N), HD (N), OD (N), RS (N), X (N), A (N, N) 
!                                                                       
!   Initialize data; change if desired                                  
!                                                                       
      DATA (HD (I), I = 1, N) / 5.D0, 5.D0, 3.D0, 5.D0 / 
      DATA (UD (I), I = 2, N) / - 2.D0, - 2.D0, - 1.D0 / 
      DATA (OD (I), I = 1, N - 1) / - 1.D0, - 1.D0, - 1.D0 / 
      DATA (RS (I), I = 1, N) / 4.D0, 2.D0, 0.D0, 4.D0 / 
!                                                                       
!   Ausgabe des Testbeispiels in Matrixform.                            
!                                                                       
      DATA ( (A (J, I), I = 1, N), J = 1, N) / M * 0.0D0 / 
      A (1, 1) = HD (1) 
      A (1, 2) = OD (1) 
      DO 100 I = 2, N - 1 
         A (I, I) = HD (I) 
         A (I, I + 1) = OD (I) 
         A (I, I - 1) = UD (I) 
  100 END DO 
      A (N, N) = HD (N) 
      A (N, N - 1) = UD (N) 
      WRITE ( *, 2010) 
      DO 110 I = 1, N 
  110 WRITE ( *, 2020) (A (I, J), J = 1, N) 
      WRITE ( *, 2030) (RS (I), I = 1, N) 
!                                                                       
      CALL TRDIG (N, UD, HD, OD, RS, X, MARKE) 
!                                                                       
      WRITE ( *, 2040) (X (I), I = 1, N) 
!                                                                       
!   Output of error message in  MARKE (-1/0/1)                          
!                                                                       
      MARKE = MARKE+2 
      GOTO (10, 20, 30), MARKE 
   10 WRITE ( * , 2050) 'ERROR:  N <= 2 !' 
      STOP 
   20 WRITE ( * , 2050) 'ERROR: MATRIX A NOT STRONGLY NONSINGULAR!' 
      STOP 
   30 WRITE ( * , 2050) 'STOP. NO ERROR!' 
      STOP 
!                                                                       
 2010 FORMAT (1X,'C[',2X,'EXAMPLE:',T78,']*',/,                         &
     &        1X,'C[',2X,8('='),T78,']*',/,                             &
     &        1X,'C[',2X,'COEFFICIENT MATRIX A:',T78,']*',/,            &
     &        1X,'C[',T78,']*')                                         
 2020 FORMAT (1X,'C[',4(1X,D12.5),T78,']*') 
 2030 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'RIGHT HAND SIDE',T78,']*',/,                  &
     &        1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',4(1X,D12.5),T78,']*',/,                           &
     &        1X,'C[',T78,']*')                                         
 2040 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'SOLUTION',T78,']*',/,                         &
     &        1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',4(1X,D12.5),T78,']*')                             
 2050 FORMAT (1X,'C[',2X,A,T78,']*') 
      END PROGRAM TEST                              
