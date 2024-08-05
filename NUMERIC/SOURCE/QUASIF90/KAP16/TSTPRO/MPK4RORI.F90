      PROGRAM TEST 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Test program for cubature over rectangular regions using the  *      
!  Romberg-Richardson method                                     *      
!                                                                *      
!  We test the subroutine K4RORI                                 *      
!                                                                *      
!  Our test example produces the output:                         *      
!                                                                *      
![                                                              ]*      
![  EXACT SOLUTION: .38682227139506E+00                         ]*      
![                                                              ]*      
![                                                              ]*      
![  IERR: 0  APPROXIMATE VALUE:  .38682227139290E+00            ]*      
![                                                              ]*      
![  ESTIMATED ERROR: .2468E-09  ACTUAL ERROR:  .2157E-11        ]*      
![                                                              ]*      
![  NUMBER OF FUNCTION CALLS:    404                            ]*      
![                                                              ]*      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  required subroutines : SXCY, K4RORI, K4ROST, RORIEX           *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author      : Volker KrÅger                                   *      
!  Date        : 6.12.1991                                       *      
!  Source code : FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
! declarations                                                          
!                                                                       
!    N = number of trapezoidal cubatures                                
!                                                                       
      PARAMETER (N = 4) 
      EXTERNAL SXCY 
      DOUBLEPRECISION WORK (0:N - 1, 2) 
      DOUBLEPRECISION A, B, C, D, CREC, DIVIAT, EXACT, SXCY 
!                                                                       
! corners of rectangle                                                  
!                                                                       
      A = 0.0D+00 
      B = 1.0D+00 
      C = 0.0D+00 
      D = 1.0D+00 
!                                                                       
! number of rectangles in X and Y directions for first                  
! trapezoidal cubature                                                  
!                                                                       
      IP = 2 
      IQ = 2 
!                                                                       
! exact solution                                                        
!                                                                       
      EXACT = (COS (B) - COS (A) ) * (SIN (C) - SIN (D) ) 
      WRITE ( *, 1000) EXACT 
!                                                                       
! approximate solution and error estimate                               
!                                                                       
      CALL K4RORI (SXCY, A, B, IP, C, D, IQ, N, CREC, DIVIAT, WORK,     &
      IERR, IUFCLL)                                                     
      WRITE ( *, 1100) IERR, CREC 
      WRITE ( *, 1200) DIVIAT, ABS (CREC - EXACT) 
      WRITE ( *, 1300) IUFCLL 
      STOP 
!                                                                       
! Format statements                                                     
!                                                                       
 1000 FORMAT (1X,'C[',T66,']*',/,                                       &
     &        1X,'C[',2X,'EXACT SOLUTION:',E20.14,T66,']*',/,           &
     &        1X,'C[',T66,']*')                                         
 1100 FORMAT (1X,'C[',T66,']*',/,                                       &
     &        1X,'C[',2X,'IERR: ',I1,2X,'APPROXIMATE VALUE: ',          &
     &        E20.14,T66,']*')                                          
 1200 FORMAT (1X,'C[',T66,']*',/,                                       &
     &        1X,'C[',2X,'ESTIMATED ERROR:',E10.4,2X,'ACTUAL ERROR: ',  &
     &        E10.4,T66,']*')                                           
 1300 FORMAT (1X,'C[',T66,']*',/,                                       &
     &        1X,'C[',2X,'NUMBER OF FUNCTION CALLS:',I7,T66,']*',/,     &
     &        1X,'C[',T66,']*')                                         
      END PROGRAM TEST                              
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION SXCY (X, Y) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Testfunktion zur Kubatur                                      *      
!                                                                *      
!     SXCY = SIN(X)*COS(Y)                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION X, Y 
!                                                                       
! Z-coordinate                                                          
!                                                                       
      SXCY = SIN (X) * COS (Y) 
      RETURN 
      END FUNCTION SXCY                             
