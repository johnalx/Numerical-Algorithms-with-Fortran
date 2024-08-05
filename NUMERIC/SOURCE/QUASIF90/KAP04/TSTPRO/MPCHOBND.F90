      PROGRAM TEST 
!                                                                       
!***********************************************************************
!                                                                       
!     Testprogramm for the subroutine CHOBND to solve a linear system   
!     A*X = RS for a symmetric positive definite  N*N matrix A.         
!                                                                       
!     For the data used, the following output is generated:             
!                                                                       
!                                                                       
![  EXAMPLE:                                                            
![  ========                                                            
![  COEFFICIENT MATRIX A:                                               
![                                                                      
![   .50000D+01   .10000D+01   .10000D+01   .00000D+00                  
![   .10000D+01   .50000D+01   .10000D+01   .10000D+01                  
![   .10000D+01   .10000D+01   .50000D+01   .10000D+01                  
![   .00000D+00   .10000D+01   .10000D+01   .50000D+01                  
![                                                                      
![  RIGHT HAND SIDE                                                     
![                                                                      
![   .10000D+02   .18000D+02   .22000D+02   .25000D+02                  
![                                                                      
![                                                                      
![  SOLUTION                                                            
![                                                                      
![   .10000D+01   .20000D+01   .30000D+01   .40000D+01                  
![  STOP. NO ERROR !                                                    
!                                                                       
!                                                                       
!     Further tests with different data are possible.                   
!                                                                       
!***********************************************************************
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      PARAMETER (N = 4, K = 2, M = N * N) 
      PARAMETER (M2 = N * (N - 1) / 2) 
      DOUBLEPRECISION A (N, N), AP (N, K + 1), Y (N), Z (N), X (N) 
!                                                                       
!   initialize; change data for other test examples                     
!                                                                       
      DATA (Y (I), I = 1, N) / 10.0D+00, 18.0D+00, 22.0D+00, 25.0D+00 / 
!     unteres und oberes Dreieck mit 1, Hauptdiagonale mit 5 vorbesetzen
      DATA ( (A (I, J), I = 1, N), J = 1, N) / M * 1.0D0 / 
                                                                        
      DO 1 I = 1, N 
         A (I, I) = 5.0D0 
    1 END DO 
      A (1, N) = 0.0D0 
      A (N, 1) = 0.0D0 
                                                                        
!                                                                       
!   Output of the test matrix                                           
!                                                                       
      WRITE ( *, 2000) 
      DO 110 I = 1, N 
         WRITE ( *, 2010) (A (I, J), J = 1, N) 
  110 END DO 
      WRITE ( *, 2020) (Y (I), I = 1, N) 
!                                                                       
!   Put matrix A into packed form                                       
!                                                                       
      DO 10 I = 1, N 
         DO 20 J = I, MIN (N, I + K) 
            AP (I, J - I + 1) = A (I, J) 
   20    END DO 
   10 END DO 
!                                                                       
      CALL CHOBND (N, K, AP, Y, X, IFLAG, Z) 
!                                                                       
      WRITE ( *, 2030) (X (I), I = 1, N) 
!                                                                       
!   Output of error code IFLAG (-1/0/1)                                 
!                                                                       
      IF (IFLAG.EQ. - 1) THEN 
      WRITE ( * , 2040) 'STOP. ERROR --- MATRIX NOT POSITIVE DEFINITE!' 
         STOP 
      ELSEIF (IFLAG.EQ.0) THEN 
         WRITE ( * , 2040) 'STOP. ERROR --- MATRIX IS SINGULAR !' 
         STOP 
      ELSE 
         WRITE ( * , 2040) 'STOP. NO ERROR !' 
         STOP 
      ENDIF 
!                                                                       
 2000 FORMAT (1X,'C[',2X,'EXAMPLE:',T78,']*',/,                         &
     &        1X,'C[',2X,8('='),T78,']*',/,                             &
     &        1X,'C[',2X,'COEFFICIENT MATRIX A:',T78,']*',/,            &
     &        1X,'C[',T78,']*')                                         
 2010 FORMAT (1X,'C[',4(1X,D12.5),T78,']*') 
 2020 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'RIGHT HAND SIDE',T78,']*',/,                  &
     &        1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',4(1X,D12.5),T78,']*',/,                           &
     &        1X,'C[',T78,']*')                                         
 2030 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'SOLUTION',T78,']*',/,                         &
     &        1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',4(1X,D12.5),T78,']*')                             
 2040 FORMAT (1X,'C[',2X,A,T78,']*') 
      END PROGRAM TEST                              
