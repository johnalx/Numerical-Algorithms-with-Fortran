      PROGRAM TEST 
!                                             (Thomas Meuser)     2/11/8
!***********************************************************************
!                                                                       
!     Test program for the subroutines  ORTOGP, ORTPOL                  
!     We compute the nodes and weights of generalized Gaussian quadratur
!                                                                       
!     The test example produces these results:                          
!                                                                       
![                                                                      
![ TEST EXAMPLE:                                                        
![ =============                                                        
![                                                                      
![ GIVEN INTEGRAL VALUES:            I        AINT(I)                   
![                               -------------------------              
![                                   0  I         1.00000               
![                                   1  I         1.00000               
![                                   2  I         2.00000               
![                                   3  I         6.00000               
![                                   4  I        24.00000               
![                                   5  I       120.00000               
![                                                                      
![ RESULTS:                                                             
![     I        NODES X(I)        WEIGHTS (I)                           
![ --------------------------------------------                         
![     1  I     .4157745568  I     .7110930099                          
![     2  I    2.2942803603  I     .2785177336                          
![     3  I    6.2899450829  I     .0103892565                          
![                                                                      
![ NO ERROR                                                             
![                                                                      
!                                                                       
!     Other tests with other data are possible.                         
!                                                                       
!     Remark: The given values of the integral are derived from Laguerre
!             polynomials.                                              
!                                                                       
!***********************************************************************
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      PARAMETER (N = 3) 
      DIMENSION AINT (0:2 * N - 1), X (N), W (N), WK (N * N + 4 * N + 1)&
      , IWK (N)                                                         
      DIMENSION WKD (N + 1) 
!                                                                       
!     initialize; switch data if desired                                
!                                                                       
      DATA AINT / 1.D0, 1.D0, 2.D0, 6.D0, 24.D0, 120.D0 / 
!                                                                       
!     Output of test example                                            
!                                                                       
      WRITE ( *, 2000) 
      WRITE ( *, 2010) (I, AINT (I), I = 0, 2 * N - 1) 
!                                                                       
!                                                                       
      CALL ORTOGP (N, AINT, X, W, IERR, WK, WKD, IWK) 
!                                                                       
!                                                                       
      IF (IERR.NE.0) THEN 
!                                                                       
!     Output of results                                                 
!                                                                       
      WRITE ( * ,  * ) 'ABANDONED COMPUTATIONS DUE TO ILL CONDITIONED MA&
     &TRIX!'                                                            
      ELSE 
         WRITE ( *, 2020) 
         WRITE ( *, 2030) (I, X (I), W (I), I = 1, N) 
         WRITE ( *, 2040) 
      ENDIF 
      STOP 
!                                                                       
!                                                                       
 2000 FORMAT (1X, 'C[ ', T74, ']*', /,                                  &
     &        1X, 'C[ TEST EXAMPLE:' ,T74, ']*', /,                     &
     &        1X, 'C[ ', 13('='), T74, ']*', /,                         &
     &        1X, 'C[', T74, ']*', /,                                   &
     &        1X, 'C[ GIVEN INTEGRAL VALUES:', 12X, 'I', 8X,            &
     &            'AINT(I)', T74, ']*', /,                              &
     &        1X, 'C[ ', 30X, 25('-'), T74, ']*')                       
 2010 FORMAT (1X, 'C[ ', 30X, I5, '  I ', F15.5, T74, ']*') 
 2020 FORMAT (1X, 'C[', T74, ']*', /,                                   &
     &        1X, 'C[ RESULTS:', T74, ']*', /,                          &
     &        1X, 'C[ ', 4X, 'I', 8X, 'NODES X(I)', 8X,                 &
     &            'WEIGHTS (I)', T74, ']*', /,                          &
     &        1X, 'C[ ', 44('-'), T74, ']*')                            
 2030 FORMAT (1X, 'C[ ', I5,'  I ', F15.10, '  I ', F15.10, T74, ']*') 
 2040 FORMAT (1X, 'C[', T74, ']*', /,                                   &
     &        1X, 'C[ ', 'NO ERROR', T74, ']*', /,                      &
     &        1X, 'C[ ', T74, ']*')                                     
      END PROGRAM TEST                              
