      PROGRAM TEST 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Testprogramm zur Kubatur Åber Dreieck-Gebiete nach Newton-    *      
!  Cotes-Formeln und                                                    
!  Test program for cubature over triangular regions using the   *      
!  Newton-Cotes formulas and Romberg-Richardson extrapolation:   *      
!                                                                *      
!  We test the subroutine K3RORI .                               *      
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
![  IERR: 0  APPROXIMATE VALUE:  .15058433969224E+00            ]*      
![           ESTIMATE ERROR: .6282E-08                          ]*      
![                                                              ]*      
![  IERR: 0  APPROXIMATE VALUE:  .23623793189023E+00            ]*      
![           ESTIMATE ERROR: .1155E-08                          ]*      
![                                                              ]*      
![  APPROXIMATE VALUE:  .38682227158247E+00  ERROR:  .1874E-09  ]*      
![                                                              ]*      
![  NUMBER OF FUNCTION CALLS:    300                            ]*      
![                                                              ]*      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines used:   SXCY, K3RORI, RORIEX, K3NEC3              *      
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
!   number of cubatures                                                 
!                                                                       
      PARAMETER (N = 4) 
      EXTERNAL SXCY 
      DOUBLEPRECISION A, B, C, D, PX, PY, QX, QY, RX, RY, CTRI, CHILF,  &
      EXACT, DIVIAT, SXCY                                               
      DOUBLEPRECISION WORK (0:N - 1, 2) 
!                                                                       
! corners of rectangle                                                  
!                                                                       
      A = 0.0D+00 
      B = 1.0D+00 
      C = 0.0D+00 
      D = 1.0D+00 
!                                                                       
! exact solution                                                        
!                                                                       
      EXACT = (COS (B) - COS (A) ) * (SIN (C) - SIN (D) ) 
      WRITE ( *, 1000) EXACT 
!                                                                       
! coordinates of the outer corners of first triangle                    
!                                                                       
      PX = A 
      PY = C 
      QX = B 
      QY = C 
      RX = A 
      RY = D 
!                                                                       
! compute triangular integral                                           
!                                                                       
      CALL K3RORI (SXCY, PX, PY, QX, QY, RX, RY, N, WORK, CHILF, DIVIAT,&
      IERR, IUFCLL)                                                     
      WRITE ( *, 1100) IERR, CHILF, DIVIAT 
      CTRI = CHILF 
!                                                                       
! coordinates of the outer corners of second triangle                   
!                                                                       
      PX = RX 
      PY = RY 
      RX = B 
      RY = D 
!                                                                       
! compute second triangular integral                                    
!                                                                       
      CALL K3RORI (SXCY, PX, PY, QX, QY, RX, RY, N, WORK, CHILF, DIVIAT,&
      IERR, IHELP)                                                      
      WRITE ( *, 1100) IERR, CHILF, DIVIAT 
!                                                                       
! add for rectangular integral approximation                            
!                                                                       
      CTRI = CTRI + CHILF 
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
     &        1X,'C[',2X,'IERR: ',I1,2X,'APPROXIMATE VALUE: ',          &
     &        E20.14,T66,']*',/,                                        &
     &        1X,'C[',11X,'ESTIMATE ERROR:',E10.4,T66,']*')             
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
