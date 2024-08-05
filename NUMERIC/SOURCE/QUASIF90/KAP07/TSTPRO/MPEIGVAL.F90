      PROGRAM TEST 
!                                                                    87/
!                                                ( Thomas Meuser )      
!***********************************************************************
!                                                                       
!     Testprogram for the subroutine EIGVAL.                            
!     Computes the maximal modulus eigenvalue and corresponding eigenvec
!     of a matrix A via vector iteration.                               
!                                                                       
!     The test example produces this result:                            
!                                                                       
![  TEST EXAMPLE:                                                       
![  =============                                                       
![  MATRIX A:                                                           
![                                                                      
![    -.20000D+01   .20000D+01   .20000D+01   .20000D+01                
![    -.30000D+01   .30000D+01   .20000D+01   .20000D+01                
![    -.20000D+01   .00000D+00   .40000D+01   .20000D+01                
![    -.10000D+01   .00000D+00   .00000D+00   .50000D+01                
![                                                                      
![  REQUIRED PARAMETERS:                                                
![                                                                      
![    MAX. 100 ITERATIONS                                               
![    ERROR BOUND:  .100E-13                                            
![                                                                      
![  SOLUTION                                                            
![                                                                      
![    EIGENVALUE:     .40000D+01                                        
![                                                                      
![    EIGENVECTOR:                                                      
![                                                                      
![        .50000D+00      .50000D+00      .50000D+00      .50000D+00    
![  NO ERROR                                                            
!                                                                       
!***********************************************************************
!                                                                       
      PARAMETER (N = 4, ND = N) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A (ND, ND), Z (ND), X (ND), Y (ND) 
!                                                                       
!     Initialize data; modify if desired.                               
!                                                                       
      DATA ( (A (J, I), I = 1, N), J = 1, N) / - 2.0D+00, 2.0D+00,      &
      2.0D+00, 2.0D+00, - 3.0D+00, 3.0D+00, 2.0D+00, 2.0D+00, - 2.0D+00,&
      0.0D+00, 4.0D+00, 2.0D+00, - 1.0D+00, 0.0D+00, 0.0D+00, 5.0D+00 / 
      DATA MO, EPSI / 100, 1.0D-14 / 
!                                                                       
!     Output of test example                                            
!                                                                       
      WRITE ( *, 2000) 
      DO 100 I = 1, N 
  100 WRITE ( *, 2010) (A (I, J), J = 1, N) 
      WRITE ( *, 2020) MO, EPSI 
!                                                                       
!                                                                       
      CALL EIGVAL (A, N, ND, MO, EPSI, X, Y, Z, EW, IER) 
!                                                                       
!                                                                       
      WRITE ( *, 2030) EW 
      WRITE ( *, 2040) (Y (I), I = 1, N) 
!                                                                       
!     Error code in  IER (0/1)                                          
!                                                                       
      IER = IER + 1 
      GOTO (10, 20), IER 
   10 WRITE ( * , 2050) 'NO ERROR' 
      STOP 
   20 WRITE ( * , 2050) 'PROGRAM CRASHED, CHECK INPUT!' 
      STOP 
!                                                                       
!                                                                       
 2000 FORMAT(1X,'C[',2X,'TEST EXAMPLE:',T78,']*',/,                     &
     &       1X,'C[',2X,13('='),T78,']*',/,                             &
     &       1X,'C[',2X,'MATRIX A:',T78,']*',/,                         &
     &       1X,'C[',T78,']*')                                          
 2010 FORMAT(1X,'C[',2X,4(1X,D12.5),T78,']*') 
 2020 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'REQUIRED PARAMETERS:',T78,']*',/,             &
     &        1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',4X,'MAX.',I4,1X,'ITERATIONS',T78,']*',/,          &
     &        1X,'C[',4X,'ERROR BOUND:',E10.3,T78,']*'/,                &
     &        1X,'C[',T78,']*' )                                        
 2030 FORMAT (1X,'C[',2X,'SOLUTION',T78,']*',/,                         &
     &        1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',4X,'EIGENVALUE:',D15.5,T78,']*')                  
 2040 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',4X,'EIGENVECTOR:',T78,']*',/,                     &
     &        1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,4(1X,D15.5),T78,']*' )                         
 2050 FORMAT (1X,'C[',2X,A,T78,']*') 
      END PROGRAM TEST                              
