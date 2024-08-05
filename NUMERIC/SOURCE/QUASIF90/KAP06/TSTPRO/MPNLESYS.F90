      PROGRAM TEST 
!                                                                    93/
!***********************************************************************
!                                                                       
!     Testprogramm for subroutine  NLESYS.                              
!     Solve a nonlinear system of equations                             
!             F(1)(X(1),...,X(N)) = 0                                   
!             F(2)(X(1),...,X(N)) = 0                                   
!             -----------------------                                   
!             F(M)(X(1),...,X(N)) = 0                                   
!     with the damped Newton iteration method. Here the Jacobi matrix is
!     approximated by central difference quotients.                     
!                                                                       
!     The test example produces:                                        
!                                                                       
![  TEST EXAMPLE:                                                       
![  ==============                                                      
![  SYSTEM OF EQUATIONS TO BE SOLVED:                                   
![                                                                      
![          F(X) = X(0)**2 - X(1) - 1                  = 0              
![          F(X) = (X(0)-2)**2 + (X(1)-0.5)**2 - 1     = 0              
![                                                                      
![  STARTING VECTOR:       .000D+00  .000D+00                           
![                                                                      
![                                                                      
![  REQUIRED PARAMETERS:                                                
![                                                                      
![    ATTENUATION =  4                                                  
![    MAX. 100 ITERATIONS                                               
![    ERROR BOUND:  .100D-05                                            
![                                                                      
![  SOLUTION:                                                           
![                                                                      
![       1.06735D+00     1.39228D-01                                    
![                                                                      
![  TERMINATED AFTER   6 ITERATION STEPS                                
!                                                                       
!     Weitere Testl„ufe sind durch Modifizierung der Werte in den       
!     DATA-Anweisungen und EXTERNAL-SUBROUTINES jederzeit m”glich.      
!                                                                       
!***********************************************************************
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      PARAMETER (N = 1, M = 1, LDDF = M) 
      EXTERNAL FX, DFX 
      CHARACTER(80) FVONX (0:M) 
      DOUBLEPRECISION X (0:N), DF (0:LDDF, 0:N + 1), WORK (0:M),        &
      S (0:N), D (0:N), F (0:N)                                         
!                                                                       
!     Initialize test example. Exchanging example is possible.          
!                                                                       
      DATA MAXIT, KMAX, EPS / 100, 4, 1.D-06 / 
      DATA (X (I), I = 0, N) / 0.0D0, 0.0D0 / 
!                                                                       
!     Prepare for output                                                
!                                                                       
      DATA (FVONX (I) , I = 0, M)  / 'X(0)**2 - X(1) - 1', '(X(0)-2)**2 &
     &+ (X(1)-0.5)**2 - 1' /                                            
!                                                                       
!     Put out test example                                              
!                                                                       
      WRITE ( *, 2000) 
      WRITE ( *, 2010) (FVONX (I), I = 0, M) 
      WRITE ( *, 2020) (X (I), I = 0, N) 
      WRITE ( *, 2040) KMAX, MAXIT, EPS 
!                                                                       
      CALL NLESYS (FX, M, N, .FALSE., DFX, DF, LDDF, MAXIT, EPS, KMAX,  &
      0, X, F, RNORM2, IERR, D, S, WORK)                                
!                                                                       
      IF (IERR.NE.0) THEN 
!                                                                       
!     Put out error message                                             
!                                                                       
         WRITE ( * , 2030) 'ERROR: IERR=', IERR 
      ELSE 
!                                                                       
!     Put out solution                                                  
!                                                                       
         WRITE ( *, 2050) (X (I), I = 0, N) 
         WRITE ( *, 2060) MAXIT 
      ENDIF 
      STOP 
!                                                                       
!                                                                       
 2000 FORMAT (1X,'C[',2X,'TEST EXAMPLE:',T78,']*',/,                    &
     &        1X,'C[',2X,14('='),T78,']*',/,                            &
     &        1X,'C[',2X,'SYSTEM OF EQUATIONS TO BE SOLVED:',T78,']*'   &
     &           ,/,1X,'C[',T78,']*')                                   
 2010 FORMAT (1X,'C[',10X,'F(X) = ',A35,' = 0',T78,']*') 
 2020 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'STARTING VECTOR:',5X,2(D10.3),T78,']*',/,     &
     &        1X,'C[',T78,']*')                                         
 2030 FORMAT (1X,'C[',2X,A,I1,T78,']*') 
 2040 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'REQUIRED PARAMETERS:',T78,']*',/,             &
     &        1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',4X,'ATTENUATION =',I3,T78,']*',/,                 &
     &        1X,'C[',4X,'MAX.',I4,1X,'ITERATIONS',T78,']*',/,          &
     &        1X,'C[',4X,'ERROR BOUND:',D10.3,T78,']*')                 
 2050 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'SOLUTION:',T78,']*',/,                        &
     &        1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,2(1X,1PD15.5),T78,']*')                        
 2060 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'TERMINATED AFTER',I4,1X,                      &
     &           'ITERATION STEPS ',T78,']*')                           
      END PROGRAM TEST                              
!                                                                       
!                                                                       
!                                                                       
      SUBROUTINE FX (X, N, F, M) 
      DOUBLEPRECISION X (0:N), F (0:M) 
      F (0) = X (0) * X (0) - X (1) - 1.D0 
      F (1) = (X (0) - 2.D0) **2.D0 + (X (1) - 0.5D0) **2.D0 - 1.D0 
      RETURN 
      END SUBROUTINE FX                             
!                                                                       
!                                                                       
      SUBROUTINE DFX (X, M, N, DF, LDDF) 
      DOUBLEPRECISION DF (0:LDDF, 0:N), X (0:N) 
      DF (0, 0) = 2.D0 * X (0) 
      DF (0, 1) = - 1.D0 
      DF (1, 0) = 2.D0 * X (0) - 4.D0 
      DF (1, 1) = 2.D0 * X (1) - 1.D0 
      RETURN 
      END SUBROUTINE DFX                            
