      PROGRAM TEST 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Test program for cubature over triangular regions using the   *      
!  NEWTON-COTES formulas:                                        *      
!                                                                *      
!  We test the subroutine K3NEC3 .                               *      
!                                                                *      
!  For F(X,Y)=SIN(X)*COS(Y) we approximate the double integral   *      
!  over a rectangular region by dividing the rectangle into two  *      
!  triangles. The sum of the two triangular integrals gives the  *      
!  value of the rectangular integration.                         *      
!                                                                *      
!  The test example gives rise to the following results:         *      
!                                                                *      
![                                                              ]*      
![  EXACT SOLUTION: .38682227139506E+00                         ]*      
![                                                              ]*      
![                                                              ]*      
![  IERR: 0  APPROXIMATE VALUE:  .15058428623848E+00            ]*      
![                                                              ]*      
![  IERR: 0  APPROXIMATE VALUE:  .23623794484260E+00            ]*      
![                                                              ]*      
![  APPROXIMATE VALUE:  .38682223108107E+00  ERROR: -.4031E-07  ]*      
![                                                              ]*      
![  NUMBER OF FUNCTION CALLS:    330                            ]*      
![                                                              ]*      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  required subroutines :    SXCY, K3NEC3                        *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author      : Volker KrÅger                                   *      
!  Date        : 6.12.1991                                       *      
!  Source code : FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
! DEclarations                                                          
!                                                                       
      EXTERNAL SXCY 
      DOUBLEPRECISION A, B, C, D, PX, PY, QX, QY, RX, RY, CTRI, CHELP,  &
      EXACT, SXCY                                                       
!                                                                       
! corners of rectangle                                                  
!                                                                       
      A = 0.0D+00 
      B = 1.0D+00 
      C = 0.0D+00 
      D = 1.0D+00 
!                                                                       
! number of triangles along one edge                                    
!                                                                       
      N = 10 
!                                                                       
! find and print exact solution                                         
!                                                                       
      EXACT = (COS (B) - COS (A) ) * (SIN (C) - SIN (D) ) 
      WRITE ( *, 1000) EXACT 
!                                                                       
! set up outer corners of first triangle                                
!                                                                       
      PX = A 
      PY = C 
      QX = B 
      QY = C 
      RX = A 
      RY = D 
!                                                                       
! compute integral                                                      
!                                                                       
      CALL K3NEC3 (SXCY, PX, PY, QX, QY, RX, RY, N, CHELP, IERR, IUFCLL) 
      WRITE ( *, 1100) IERR, CHELP 
      CTRI = CHELP 
!                                                                       
! set up outer corners of second triangle                               
!                                                                       
      PX = RX 
      PY = RY 
      RX = B 
      RY = D 
!                                                                       
! compute integral                                                      
!                                                                       
      CALL K3NEC3 (SXCY, PX, PY, QX, QY, RX, RY, N, CHELP, IERR, IHELP) 
      WRITE ( *, 1100) IERR, CHELP 
!                                                                       
! add results for approximation of rectangular integral                 
!                                                                       
      CTRI = CTRI + CHELP 
      IUFCLL = IHELP + IUFCLL 
      WRITE ( *, 1200) CTRI, CTRI - EXACT 
      WRITE ( *, 1300) IUFCLL 
      STOP 
!                                                                       
! Format statements                                                     
!                                                                       
 1000 FORMAT (1X,'C[',T66,']*',/,                                       &
     &        1X,'C[',2X,'EXACT SOLUTION:',E20.14,T66,']*',/,           &
     &        1X,'C[',T66,']*')                                         
 1100 FORMAT (1X,'C[',T66,']*',/,                                       &
     &        1X,'C[',2X,'IERR: ',I1,2X,                                &
     &        'APPROXIMATE VALUE: ',E20.14,T66,']*')                    
 1200 FORMAT (1X,'C[',T66,']*',/,                                       &
     &        1X,'C[',2X,'APPROXIMATE VALUE: ',E20.14,2X,'ERROR: ',     &
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
!  Test function for cubature                                    *      
!                                                                *      
!     SXCY = SIN(X)*COS(Y)                                       *      
!                                                                *      
!                                                                *      
!*****************************************************************      
!                                                                       
!                                                                       
      DOUBLEPRECISION X, Y 
!                                                                       
! Z-coordinate                                                          
!                                                                       
      SXCY = SIN (X) * COS (Y) 
      RETURN 
      END FUNCTION SXCY                             
