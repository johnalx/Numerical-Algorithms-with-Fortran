      PROGRAM TEST 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Testprogram for subroutine  ZEROIN.                           *      
!                                                                *      
!  F(X)=(X*X-2.)*X-5.                                            *      
!                                                                *      
!  1st start at : A=2.                                           *      
!  2nd start at : B=3.                                           *      
!  Absolute error : ABSERR=0                                     *      
!  Relative error : RELERR=5.E-13                                *      
!                                                                *      
![  EXAMPLE SUBROUTINE ZEROIN                                   ]*      
![  =========================                                   ]*      
![                                                              ]*      
![  FUNCTIONAL ZERO AT:               .20945514815423E+01       ]*      
![  FUNCTION VALUE AT COMPUTED ZERO:  .40453751459779E-14       ]*      
![  NUMBER OF EVALUATIONS:  8                                   ]*      
![  IFEHL = 2                                                   ]*      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      EXTERNAL F 
      A = 2. 
      B = 3. 
      NMAX = 50 
      ABSERR = 0. 
      RELERR = 5.0E-13 
      CALL ZEROIN (F, ABSERR, RELERR, NMAX, 6, A, B, FB, IANZ, IFEHL) 
      WRITE ( *, 900) B, FB, IANZ, IFEHL 
  900 FORMAT(/,                                                         &
     &       1X,'C[  EXAMPLE SUBROUTINE ZEROIN',T66,']*',/,             &
     &       1X,'C[  ',25('='),T66,']*',/,                              &
     &       1X,'C[',T66,']*',/,                                        &
     &       1X,'C[  FUNCTIONAL ZERO AT: ',T39,E20.14,T66,']*',/,       &
     &       1X,'C[  FUNCTION VALUE AT COMPUTED ZERO: ',T39,            &
     &          E20.14,T66,']*',/,                                      &
     &       1X,'C[  NUMBER OF EVALUATIONS:',I3,T66,']*',/,             &
     &       1X,'C[  IFEHL =',I2,T66,']*')                              
      STOP 
      END PROGRAM TEST                              
                                                                        
      DOUBLEPRECISION FUNCTION F (X) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      F = (X * X - 2.) * X - 5. 
      RETURN 
      END FUNCTION F                                
