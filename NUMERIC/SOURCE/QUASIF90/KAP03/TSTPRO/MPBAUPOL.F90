      PROGRAM TEST 
!                                       4.21.1992  (Dubois Guido)       
!*****************************************************************      
!                                                                *      
!     Testprogram for the subroutine BAUPOL.                     *      
!     We compute all real and complex roots of a polynomial      *      
!     P(Z) = A(N) * Z**N +...+ A(1) * Z + A(0) with complex      *      
!     coefficients.                                              *      
!                                                                *      
!     With the given polynomial                                  *      
!     P(Z) = (Z-3*I-1) * (Z-2*I+1)**2 * (Z+4*I-2)                *      
!     the test results are as follows:                           *      
!                                                                *      
![    EXAMPLE OF SUBROUTINE BAUPOL                              ]*      
![    ============================                              ]*      
![                                                              ]*      
![    COEFFICIENTS:                                             ]*      
![                                                              ]*      
![     I   RE(A(I))   IM(A(I))                                  ]*      
![    -------------------------                                 ]*      
![     0  1.000E+00  0.000E+00                                  ]*      
![     1 -1.000E+00 -3.000E+00                                  ]*      
![     2  9.000E+00  1.200E+01                                  ]*      
![     3  4.900E+01 -4.300E+01                                  ]*      
![     4 -3.400E+01 -6.200E+01                                  ]*      
![                                                              ]*      
![    ZEROS:                                                    ]*      
![                                                              ]*      
![     I   RE(Z(I))   IM(Z(I))                                  ]*      
![    -------------------------                                 ]*      
![     1 -2.000E-01 -4.000E-01                                  ]*      
![     2  1.000E-01 -3.000E-01                                  ]*      
![     3  1.000E-01  2.000E-01                                  ]*      
![     4 -2.000E-01 -4.000E-01                                  ]*      
![                                                              ]*      
![       43 ITERATIONS                                          ]*      
!                                                                *      
!*****************************************************************      
!                                                                       
      PARAMETER (N = 4) 
      DOUBLEPRECISION KOEFRE (0:100), KOEFIM (0:100), WURZRE (1:100),   &
      WURZIM (1:100)                                                    
!                                                                       
!     Initialize; other test problems can be inserted here for          
!     differing N                                                       
!                                                                       
      DATA KOEFRE / 1.0D0, - 1.0D0, 9.0D0, 49.0D0, - 34.0D0, 96 * 0.0D0 &
      /                                                                 
      DATA KOEFIM / 0.0D0, - 3.0D0, 12.0D0, - 43.0D0, - 62.0D0, 96 *    &
      0.0D0 /                                                           
!                                                                       
!     Output of Test coefficients                                       
!                                                                       
      WRITE ( *, 100) 
      WRITE ( *, 110) (I, KOEFRE (I), KOEFIM (I), I = 0, N) 
      CALL BAUPOL (KOEFRE, KOEFIM, N, .TRUE., WURZRE, WURZIM, ISHRIT) 
!                                                                       
!     Output of Test results                                            
!                                                                       
      WRITE ( *, 105) 
      WRITE ( *, 110) (I, WURZRE (I), WURZIM (I), I = 1, N) 
      WRITE ( *, 120) ISHRIT 
      STOP 
  100 FORMAT(1X,'C[',4X,'EXAMPLE OF SUBROUTINE BAUPOL',T66,']*',/,      &
     &       1X,'C[',4X,28('='),T66,']*',/,                             &
     &       1X,'C[',T66,']*',/,                                        &
     &       1X,'C[',4X,'COEFFICIENTS:',T66,']*',/,                     &
     &       1X,'C[',T66,']*',/,                                        &
     &       1X,'C[',4X,' I   RE(A(I))   IM(A(I))',T66,']*',/,          &
     &       1X,'C[',4X,25('-'),T66,']*')                               
  110 FORMAT(1X,'C[',4X,I2,1X,1PE10.3,1X,E10.3,T66,']*') 
  105 FORMAT(1X,'C[',T66,']*',/,                                        &
     &       1X,'C[',4X,'ZEROS:',T66,']*',/,                            &
     &       1X,'C[',T66,']*',/,                                        &
     &       1X,'C[',4X,' I   RE(Z(I))   IM(Z(I))',T66,']*',/,          &
     &       1X,'C[',4X,25('-'),T66,']*')                               
  120 FORMAT(1X,'C[',T66,']*',/,                                        &
     &       1X,'C[',4X,I5,' ITERATIONS',T66,']*')                      
      END PROGRAM TEST                              
