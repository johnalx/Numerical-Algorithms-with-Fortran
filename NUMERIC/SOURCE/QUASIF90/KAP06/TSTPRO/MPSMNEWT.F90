      PROGRAM TEST 
!                                                                    10/
!                                                ( Thomas Meuser )      
!***********************************************************************
!                                                                       
!     Testprogramm for subroutine SMNEWT.                               
!     Solve a nonlinear system of equations                             
!             F(1)(X(1),...,X(N)) = 0                                   
!             F(2)(X(1),...,X(N)) = 0                                   
!             -----------------------                                   
!             F(N)(X(1),...,X(N)) = 0                                   
!     with damped Newton iterations.                                    
!                                                                       
!     The test example produces:                                        
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
![  THE ITERATION STEPS ARE SAVED ON TAPE 5 !                           
![                                                                      
![  REQUIRED PARAMETERS:                                                
![                                                                      
![    ATTENUATION =  0                                                  
![    MAX. 100 ITERATIONS                                               
![    ERROR BOUND:  .100D-05                                            
![                                                                      
![  SOLUTION:                                                           
![                                                                      
![       1.06735D+00     1.39228D-01                                    
![                                                                      
![  TERMINATED AFTER   6 ITERATION STEPS                                
!                                                                       
!     Other test runs with different data are possible.                 
!                                                                       
!***********************************************************************
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      PARAMETER (N = 2, IDF = N) 
      EXTERNAL FX, DFX 
      CHARACTER(80) FVONX (N) 
      INTEGER IWORK (N) 
      DOUBLEPRECISION X (N), DF (IDF, N), WORK (4 * N), F (N) 
!                                                                       
!     Initialize test example, change if desired                        
!                                                                       
      DATA ITAPE, MAXIT, KMAX, EPS / 5, 100, 0, 1D-06 / 
      DATA (X (I), I = 1, N) / 0.D0, 0.D0 / 
!                                                                       
!     Prepare for output                                                
!                                                                       
      DATA (FVONX (I) , I = 1, N)  / 'X(1)**2 - X(2) - 1', '(X(1)-2)**2 &
     &+ (X(2)-0.5)**2 - 1' /                                            
!                                                                       
!     Put out test example                                              
!                                                                       
      WRITE ( *, 2000) 
      WRITE ( *, 2010) (FVONX (I), I = 1, N) 
      WRITE ( *, 2020) (X (I), I = 1, N) 
      WRITE ( *, 2030) ITAPE 
      WRITE ( *, 2040) KMAX, MAXIT, EPS 
!                                                                       
      OPEN (ITAPE, FILE = 'TAPE5') 
!                                                                       
      CALL SMNEWT (FX, DFX, N, MAXIT, IFEHL, KMAX, ITAPE, EPS, RNORM2,  &
      F, X, DF, IDF, IWORK, WORK)                                       
!                                                                       
      CLOSE (ITAPE) 
!                                                                       
!     Put out error message in  IFEHL (1,2,3)                           
!                                                                       
      IF (IFEHL.EQ.1) THEN 
         WRITE ( *, 1000) 
      ELSEIF (IFEHL.EQ.2) THEN 
         WRITE ( *, 1100) 
      ELSEIF (IFEHL.EQ.3) THEN 
         WRITE ( *, 1200) 
      ELSE 
!                                                                       
!     Put out solution                                                  
!                                                                       
         WRITE ( *, 2050) (X (I), I = 1, N) 
         WRITE ( *, 2060) MAXIT 
      ENDIF 
      STOP 
!                                                                       
 1000 FORMAT (1X,'AFTER MAX STEPS ERROR BOUND WAS NOT REACHED!') 
 1100 FORMAT (1X,'ERROR IN THE SOLVING THE LINEAR SYSTEM OF',           &
     &           'EQUATIONS (JACOBI MATRIX IS SINGULAR!)')              
 1200 FORMAT (1X,'INCORRECT INPUT PARAMETER!') 
 2000 FORMAT (1X,'C[',2X,'TEST EXAMPLE:',T78,']*',/,                    &
     &        1X,'C[',2X,14('='),T78,']*',/,                            &
     &        1X,'C[',2X,'SYSTEM OF EQUATIONS TO BE SOLVED:',T78,']*'   &
     &           ,/,1X,'C[',T78,']*')                                   
 2010 FORMAT (1X,'C[',10X,'F(X) = ',A35,' = 0',T78,']*') 
 2020 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'STARTING VECTOR:',5X,2(D10.3),T78,']*',/,     &
     &        1X,'C[',T78,']*')                                         
 2030 FORMAT (1X,'C[',T78,']*',/,                                       &
     &        1X,'C[',2X,'THE ITERATION STEPS ARE SAVED ON TAPE',       &
     &           I2,1X,'!',T78,']*')                                    
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
      SUBROUTINE FX (N, X, F) 
      DOUBLEPRECISION X (N), F (N) 
      F (1) = X (1) * X (1) - X (2) - 1.D0 
      F (2) = (X (1) - 2.D0) **2 + (X (2) - 0.5D0) **2 - 1.D0 
      RETURN 
      END SUBROUTINE FX                             
!                                                                       
!                                                                       
      SUBROUTINE DFX (N, X, DF, IDF) 
      DOUBLEPRECISION DF (IDF, N), X (N) 
      DF (1, 1) = 2.D0 * X (1) 
      DF (1, 2) = - 1.D0 
      DF (2, 1) = 2.D0 * X (1) - 4.D0 
      DF (2, 2) = 2.D0 * X (2) - 1.D0 
      RETURN 
      END SUBROUTINE DFX                            
