      PROGRAM TEST 
!                                       9.19.1992  (Dubois Guido)       
!*****************************************************************      
!                                                                *      
!     Testprogram for the subroutine LAGUER.                     *      
!     Compute the zeros of the polynomial                        *      
!     PN(X) = A(N) * X**N +...+ A(1) * X + A(0)  with real roots *      
!     using the method of Laguerre.                              *      
!                                                                *      
!     Test polynomial:                                           *      
!     PN(X) = X**4 - 10 X**3 + 35 X**2 - 50 X + 24               *      
!     Test results:                                              *      
!                                                                *      
![    EXAMPLE OF SUBROUTINE LAGUER                              ]*      
![    ============================                              ]*      
![                                                              ]*      
![    COEFFICIENTS:                                             ]*      
![                                                              ]*      
![     I     A(I))                                              ]*      
![    --------------                                            ]*      
![     0     24.000                                             ]*      
![     1    -50.000                                             ]*      
![     2     35.000                                             ]*      
![     3    -10.000                                             ]*      
![     4      1.000                                             ]*      
![                                                              ]*      
![    ZEROS:                                                    ]*      
![                                                              ]*      
![     I     XI(I)    NITER                                     ]*      
![    ---------------------                                     ]*      
![     1      1.000     4                                       ]*      
![     2      2.000     4                                       ]*      
![     3      3.000     0                                       ]*      
![     4      4.000     0                                       ]*      
!                                                                *      
!*****************************************************************      
!                                                                       
      PARAMETER (N = 4) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION A (0:N), XI (1:N), WORK (0:N) 
      INTEGER NITER (1:N) 
!                                                                       
!     Initialize; enter other polynomial if desired                     
!                                                                       
      DATA A / 24.0D0, - 50.0D0, 35.0D0, - 10.0D0, 1.0D0 / 
      DATA ABSERR / 0.0D0 / 
      DATA RELERR / 1.D-10 / 
      DATA IMAX / 100 / 
!                                                                       
!     Output of test coefficients                                       
!                                                                       
      WRITE ( *, 900) 
      WRITE ( *, 910) (I, A (I), I = 0, N) 
      CALL LAGUER (A, N, ABSERR, RELERR, IMAX, XI, NITER, IANZ, WORK,   &
      IERR)                                                             
!                                                                       
!     Output of test results                                            
!                                                                       
      IF (IERR.EQ.1) THEN 
         WRITE ( *, 920) 
         WRITE ( *, 930) (I, XI (I), NITER (I), I = 1, N) 
      ELSE 
         WRITE ( *, 940) IERR 
      ENDIF 
      STOP 
  900 FORMAT(1X,'C[',4X,'EXAMPLE OF SUBROUTINE LAGUER',T66,']*',/,      &
     &       1X,'C[',4X,28('='),T66,']*',/,                             &
     &       1X,'C[',T66,']*',/,                                        &
     &       1X,'C[',4X,'COEFFICIENTS:',T66,']*',/,                     &
     &       1X,'C[',T66,']*',/,                                        &
     &       1X,'C[',4X,' I     A(I))',T66,']*',/,                      &
     &       1X,'C[',4X,14('-'),T66,']*')                               
  910 FORMAT(1X,'C[',4X,I2,1X,F10.3,T66,']*') 
  920 FORMAT(1X,'C[',T66,']*',/,                                        &
     &       1X,'C[',4X,'ZEROS:',T66,']*',/,                            &
     &       1X,'C[',T66,']*',/,                                        &
     &       1X,'C[',4X,' I     XI(I)    NITER',T66,']*',/,             &
     &       1X,'C[',4X,21('-'),T66,']*')                               
  930 FORMAT(1X,'C[',4X,I2,1X,F10.3,3X,I3,T66,']*') 
  940 FORMAT(1X,'C[',T66,']*',/,                                        &
     &       1X,'C[',4X,'ERROR WHEN EXECUTING SUBROUTINE LAGUER:',      &
     &       I2,T66,']*')                                               
      END PROGRAM TEST                              
