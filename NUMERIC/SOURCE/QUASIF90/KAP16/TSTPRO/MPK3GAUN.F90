      PROGRAM TEST 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Test program for cubature over triangular regions using the   *      
!  N point GAUSS formula:                                        *      
!                                                                *      
!  We test the subroutine K3GAUN .                               *      
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
![  METHOD: 1  IERR: 0  APPROXIMATE VALUE:  .15062615253663E+00 ]*      
![                                                              ]*      
![  METHOD: 1  IERR: 0  APPROXIMATE VALUE:  .23630357146945E+00 ]*      
![                                                              ]*      
![  APPROXIMATE VALUE:  .38692972400608E+00  ERROR:  .1075E-03  ]*      
![                                                              ]*      
![  NUMBER OF FUNCTION CALLS:    200                            ]*      
![                                                              ]*      
![                                                              ]*      
![  METHOD: 2  IERR: 0  APPROXIMATE VALUE:  .15062615253663E+00 ]*      
![                                                              ]*      
![  METHOD: 2  IERR: 0  APPROXIMATE VALUE:  .23627075228868E+00 ]*      
![                                                              ]*      
![  APPROXIMATE VALUE:  .38689690482531E+00  ERROR:  .7463E-04  ]*      
![                                                              ]*      
![  NUMBER OF FUNCTION CALLS:    400                            ]*      
![                                                              ]*      
![                                                              ]*      
![  METHOD: 3  IERR: 0  APPROXIMATE VALUE:  .15058434344993E+00 ]*      
![                                                              ]*      
![  METHOD: 3  IERR: 0  APPROXIMATE VALUE:  .23623792745072E+00 ]*      
![                                                              ]*      
![  APPROXIMATE VALUE:  .38682227090064E+00  ERROR: -.4944E-09  ]*      
![                                                              ]*      
![  NUMBER OF FUNCTION CALLS:    600                            ]*      
![                                                              ]*      
![                                                              ]*      
![  METHOD: 7  IERR: 0  APPROXIMATE VALUE:  .15058433946991E+00 ]*      
![                                                              ]*      
![  METHOD: 7  IERR: 0  APPROXIMATE VALUE:  .23623793192523E+00 ]*      
![                                                              ]*      
![  APPROXIMATE VALUE:  .38682227139514E+00  ERROR:  .8815E-13  ]*      
![                                                              ]*      
![  NUMBER OF FUNCTION CALLS:   1400                            ]*      
![                                                              ]*      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  required subroutines :    SXCY, K3GAUN                        *      
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
      EXTERNAL SXCY 
      DOUBLEPRECISION A, B, C, D, PX, PY, QX, QY, RX, RY, CTRI, CHELP,  &
      EXACT, WORK (3, 0:6), SXCY                                        
!                                                                       
! corners of rectangle                                                  
!                                                                       
      A = 0.0D+00 
      B = 1.0D+00 
      C = 0.0D+00 
      D = 1.0D+00 
!                                                                       
! initialize MOLD for first call                                        
!                                                                       
      MOLD = - 1 
!                                                                       
! number of triangles along one rectangular edge                        
!                                                                       
      N = 10 
!                                                                       
! compute exact solution and put out                                    
!                                                                       
      EXACT = (COS (B) - COS (A) ) * (SIN (C) - SIN (D) ) 
      WRITE ( *, 1000) EXACT 
!                                                                       
! Compute approximations                                                
!                                                                       
      DO 10 METHOD = 1, 7 
         IF (METHOD.LE.3.OR.METHOD.EQ.7) THEN 
!                                                                       
! set up coordinates of outer corner of first triangle                  
!                                                                       
            PX = A 
            PY = C 
            QX = B 
            QY = C 
            RX = A 
            RY = D 
!                                                                       
! compute approximations and print out                                  
!                                                                       
            CALL K3GAUN (SXCY, PX, PY, QX, QY, RX, RY, N, METHOD, MOLD, &
            CHELP, WORK, IERR, IUFCLL)                                  
            WRITE ( *, 1100) METHOD, IERR, CHELP 
            CTRI = CHELP 
!                                                                       
! set up coordinates of outer corner of second triangle                 
!                                                                       
            PX = RX 
            PY = RY 
            RX = B 
            RY = D 
!                                                                       
! compute approximate value and put out                                 
!                                                                       
            CALL K3GAUN (SXCY, PX, PY, QX, QY, RX, RY, N, METHOD, MOLD, &
            CHELP, WORK, IERR, IHELP)                                   
            WRITE ( *, 1100) METHOD, IERR, CHELP 
!                                                                       
! compute sum of two triangular integrations and put out                
!                                                                       
            CTRI = CTRI + CHELP 
            IUFCLL = IHELP + IUFCLL 
            WRITE ( *, 1200) CTRI, CTRI - EXACT 
            WRITE ( *, 1300) IUFCLL 
         ENDIF 
   10 END DO 
      STOP 
!                                                                       
! Format statements                                                     
!                                                                       
 1000 FORMAT (1X,'C[',T66,']*',/,                                       &
     &        1X,'C[',2X,'EXACT SOLUTION:',E20.14,T66,']*',/,           &
     &        1X,'C[',T66,']*')                                         
 1100 FORMAT (1X,'C[',T66,']*',/,                                       &
     &        1X,'C[',2X,'METHOD: ',I1,2X,'IERR: ',I1,2X,               &
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
!  Test function for cubature                                           
!                                                                *      
!     SXCY = SIN(X)*COS(Y)                                       *      
!                                                                *      
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
