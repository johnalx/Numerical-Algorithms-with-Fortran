      PROGRAM TEST 
!                                                ( Thomas Meuser )      
!***********************************************************************
!                                                                       
!     Testprogram for the subroutines FDIAG, FDIAGP, FDIAGS.            
!     Solve a linear system  A*X = RS with strongly nonsingular five dia
!     matrix A.                                                         
!                                                                       
!     The test example gives the following results:                     
!                                                                       
![  EXAMPLE:                                                            
![  ========                                                            
![  COEFFICIENT MATRIX A:                                               
![                                                                      
![   .30000D+01   .20000D+01   .10000D+01   .00000D+00   .00000D+00     
![   .10000D+01   .60000D+01   .30000D+01   .10000D+01   .00000D+00     
![   .10000D+01   .20000D+01   .30000D+01   .10000D+01   .20000D+01     
![   .00000D+00   .10000D+01   .10000D+01   .40000D+01   .10000D+01     
![   .00000D+00   .00000D+00   .10000D+01   .20000D+01   .50000D+01     
![                                                                      
![  RIGHT HAND SIDE                                                     
![                                                                      
![   .60000D+01   .11000D+02   .90000D+01   .70000D+01   .80000D+01     
![                                                                      
![                                                                      
![  SOLUTION                                                            
![                                                                      
![   .10000D+01   .10000D+01   .10000D+01   .10000D+01   .10000D+01     
![  STOP. NO ERROR!                                                     
!                                                                       
!     Other examples can be created by modifying the test data.         
!                                                                       
!***********************************************************************
!                                                                       
      PARAMETER (N = 5, M = N * N) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION UD2 (N), UD1 (N), HD (N), OD1 (N), OD2 (N),       &
      RS (N), X (N), A (N, N)                                           
!                                                                       
!   Initialize test example; modify data if desired                     
!                                                                       
      DATA (UD2 (I), I = 3, N) / 1.D0, 1.D0, 1.D0 / 
      DATA (UD1 (I), I = 2, N) / 1.D0, 2.D0, 1.D0, 2.D0 / 
      DATA (HD (I), I = 1, N) / 3.D0, 6.D0, 3.D0, 4.D0, 5.D0 / 
      DATA (OD1 (I), I = 1, N - 1) / 2.D0, 3.D0, 1.D0, 1.D0 / 
      DATA (OD2 (I), I = 1, N - 2) / 1.D0, 1.D0, 2.D0 / 
      DATA (RS (I), I = 1, N) / 6.D0, 11.D0, 9.D0, 7.D0, 8.D0 / 
!                                                                       
!   Output of test example in matrix form                               
!                                                                       
      DATA ( (A (J, I), I = 1, N), J = 1, N) / M * 0.0D0 / 
      A (1, 1) = HD (1) 
      A (1, 2) = OD1 (1) 
      A (1, 3) = OD2 (1) 
      A (2, 1) = UD1 (2) 
      A (2, 2) = HD (2) 
      A (2, 3) = OD1 (2) 
      A (2, 4) = OD2 (2) 
      DO 110 I = 3, N - 2 
         A (I, I - 2) = UD2 (I) 
         A (I, I - 1) = UD1 (I) 
         A (I, I) = HD (I) 
         A (I, I + 1) = OD1 (I) 
         A (I, I + 2) = OD2 (I) 
  110 END DO 
      A (N - 1, N - 3) = UD2 (N - 1) 
      A (N - 1, N - 2) = UD1 (N - 1) 
      A (N - 1, N - 1) = HD (N - 1) 
      A (N - 1, N) = OD1 (N - 1) 
      A (N, N - 2) = UD2 (N) 
      A (N, N - 1) = UD1 (N) 
      A (N, N) = HD (N) 
      WRITE ( *, 2000) 
      DO 120 I = 1, N 
  120 WRITE ( *, 2010) (A (I, J), J = 1, N) 
      WRITE ( *, 2020) (RS (I), I = 1, N) 
!                                                                       
      CALL FDIAG (N, UD2, UD1, HD, OD1, OD2, RS, X, MARKE) 
!                                                                       
      WRITE ( *, 2030) (X (I), I = 1, N) 
!                                                                       
!   Output of error code  MARKE (-1,0,1)                                
!                                                                       
      MARKE = MARKE+2 
      GOTO (10, 20, 30), MARKE 
   10 WRITE ( * , 2040) 'ERROR:  N <= 3 !' 
      STOP 
   20 WRITE ( * , 2040) 'ERROR: MATRIX A NOT STRONGLY NONSINGULAR!' 
      STOP 
   30 WRITE ( * , 2040) 'STOP. NO ERROR!' 
      STOP 
!                                                                       
 2000 FORMAT (1X,'C[',2X,'EXAMPLE:',T78,']*',/,                         &
     &        1X,'C[',2X,8('='),T78,']*',/,                             &
     &        1X,'C[',2X,'COEFFICIENT MATRIX A:',T78,']*',/,            &
     &        1X,'C[',T78,']*')                                         
 2010 FORMAT (1X,'C[',5(1X,D12.5),T78,']*') 
 2020 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'RIGHT HAND SIDE',T78,']*',/,                  &
     &        1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',5(1X,D12.5),T78,']*',/,                           &
     &        1X,'C[',T78,']*')                                         
 2030 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'SOLUTION',T78,']*',/,                         &
     &        1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',5(1X,D12.5),T78,']*')                             
 2040 FORMAT (1X,'C[',2X,A,T78,']*') 
      END PROGRAM TEST                              
