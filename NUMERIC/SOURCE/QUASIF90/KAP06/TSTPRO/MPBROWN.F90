      PROGRAM TEST 
!                                                                     87
!                                                ( Thomas Meuser )      
!***********************************************************************
!                                                                       
!     Testprogram for subroutine BROWN.                                 
!     Solve a nonlinear system of equations                             
!             F(1)(X(1),...,X(N)) = 0                                   
!             F(2)(X(1),...,X(N)) = 0                                   
!             -----------------------                                   
!             F(N)(X(1),...,X(N)) = 0                                   
!     using the method of Brown.                                        
!                                                                       
!     The test example produces this output;                            
!                                                                       
![  TEST EXAMPLE:                                                       
![  ==============                                                      
![  SYSTEM OF EQUATIONS TO BE SOLVED:                                   
![                                                                      
![          F(X) = X(1)**2 - X(2) - 1                  = 0              
![          F(X) = (X(1)-2)**2 + (X(2)-0.5)**2 - 1     = 0              
![                                                                      
![  STARTING VECTOR:       .000D+00  .000D+00                           
![                                                                      
![                                                                      
![  REQUIRED PARAMETERS:                                                
![                                                                      
![    ERROR BOUND:  .100D-13                                            
![    MACHINE CONSTANT:  .100D-13                                       
![    MAX. 100 ITERATIONS                                               
![                                                                      
![  NO OUTPUT OF THE APPROXIMATE SOLUTION AFTER EACH ITERATION          
![                                                                      
![  SOLUTION:                                                           
![                                                                      
![       1.06735D+00     1.39228D-01                                    
![                                                                      
![  TERMINATED AFTER   8 ITERATION STEPS                                
!                                                                       
!     Further tests are possible by changing the data, N and the        
!     EXTERNAL SUBROUTINE.                                              
!                                                                       
!***********************************************************************
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      PARAMETER (N = 2, NIHF = N, NHF = N) 
      EXTERNAL FKT 
      DOUBLEPRECISION X0 (N), X1 (N), IHF (NIHF, N + 1), HF (NHF, N + 3) 
      CHARACTER(80) FVONX (N) 
!                                                                       
!     Initialize data; other test examples with different data, N and   
!     EXTERNAL SUBROUTINE                                               
!                                                                       
      DATA EPS, EPSM, IAUS, MAXIT / 1D-14, 1D-14, 0, 100 / 
      DATA (X0 (I), I = 1, N) / 0.D0, 0.D0 / 
!                                                                       
!     Initialize output                                                 
!                                                                       
      DATA (FVONX (I) , I = 1, N)  / 'X(1)**2 - X(2) - 1', '(X(1)-2)**2 &
     &+ (X(2)-0.5)**2 - 1' /                                            
!                                                                       
!     Print test example                                                
!                                                                       
      WRITE ( *, 2000) 
      WRITE ( *, 2010) (FVONX (I), I = 1, N) 
      WRITE ( *, 2020) (X0 (I), I = 1, N) 
      WRITE ( *, 2030) EPS, EPSM, MAXIT 
      IF (IAUS.EQ.0) THEN 
         WRITE ( *, 2040) 
      ELSE 
         WRITE ( *, 2045) 
      ENDIF 
!                                                                       
!                                                                       
      CALL BROWN (N, FKT, X0, EPS, EPSM, IAUS, MAXIT, IHF, NIHF, HF,    &
      NHF, X1, ITANZ, IFEHL)                                            
!                                                                       
!                                                                       
      IF (IFEHL.GT.0) THEN 
!                                                                       
!     Put out error message in  IFEHL (1/2)                             
!                                                                       
         GOTO (10, 20), IFEHL 
   10    STOP 'AFTER MAXIT STEPS ERROR BOUND WAS NOT REACHED!' 
   20    STOP 'NO COMPUTATIONS, JACOBI MATRIX IS SINGULAR!' 
      ELSE 
!                                                                       
!     OUTPUT of solution                                                
!                                                                       
         WRITE ( *, 2050) (X1 (I), I = 1, N) 
         WRITE ( *, 2060) ITANZ 
      ENDIF 
!                                                                       
!                                                                       
 2000 FORMAT (1X,'C[',2X,'TEST EXAMPLE:',T78,']*',/,                    &
     &        1X,'C[',2X,14('='),T78,']*',/,                            &
     &        1X,'C[',2X,'SYSTEM OF EQUATIONS TO BE SOLVED:',T78,']*'   &
     &        ,/,1X,'C[',T78,']*')                                      
 2010 FORMAT (1X,'C[',10X,'F(X) = ',A35,' = 0',T78,']*') 
 2020 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'STARTING VECTOR:',5X,2(D10.3),T78,']*',/,     &
     &        1X,'C[',T78,']*')                                         
 2030 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'REQUIRED PARAMETERS:',T78,']*',/,             &
     &        1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',4X,'ERROR BOUND:',D10.3,T78,']*',/,               &
     &        1X,'C[',4X,'MACHINE CONSTANT:',D10.3,T78,']*',/,          &
     &        1X,'C[',4X,'MAX.',I4,1X,'ITERATIONS',T78,']*',/,          &
     &        1X,'C[',T78,']*')                                         
 2040 FORMAT (1X,'C[',2X,'NO OUTPUT OF THE APPROXIMATE SOLUTION ',      &
     &                   'AFTER EACH ITERATION',T78,']*')               
 2045 FORMAT (1X,'C[',2X,'OUTPUT OF THE APPROXIMATE SOLUTION ',         &
     &                   'FOLLOWING EACH ITERATION',T78,']*')           
 2050 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'SOLUTION:',T78,']*',/,                        &
     &        1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,2(1X,1PD15.5),T78,']*')                        
 2060 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'TERMINATED AFTER',I4,1X,                      &
     &                   'ITERATION STEPS ',T78,']*')                   
      END PROGRAM TEST                              
!                                                                       
!                                                                       
      SUBROUTINE FKT (K, X, F) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION X (2) 
      GOTO (1, 2) K 
    1 F = X (1) * X (1) - X (2) - 1.D0 
      RETURN 
    2 F = (X (1) - 2.D0) **2 + (X (2) - 0.5D0) **2 - 1.D0 
      RETURN 
      END SUBROUTINE FKT                            
