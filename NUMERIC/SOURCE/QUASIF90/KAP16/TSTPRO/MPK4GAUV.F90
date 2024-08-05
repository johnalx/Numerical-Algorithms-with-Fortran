      PROGRAM TEST 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Test program for cubature over rectangular region using the   *      
!  Gaussian formulas:                                            *      
!                                                                *      
!  We test the subroutine K4GAUV                                 *      
!                                                                *      
!  Our test example produces the results:                        *      
!                                                                *      
![                                                              ]*      
![  EXACT SOLUTION: .38682227139506E+00                         ]*      
![                                                              ]*      
![                                                              ]*      
![  METHOD: 0  IERR: 0  APPROXIMATE VALUE:  .38732633996768E+00 ]*      
![                                                              ]*      
![  ESTIMATED ERROR:-.5056E-03  ACTUAL ERROR:  .5041E-03        ]*      
![                                                              ]*      
![  NUMBER OF FUNCTION CALLS:     80                            ]*      
![                                                              ]*      
![                                                              ]*      
![  METHOD: 1  IERR: 0  APPROXIMATE VALUE:  .38682222765159E+00 ]*      
![                                                              ]*      
![  ESTIMATED ERROR: .2191E-06  ACTUAL ERROR: -.4374E-07        ]*      
![                                                              ]*      
![  NUMBER OF FUNCTION CALLS:    320                            ]*      
![                                                              ]*      
![                                                              ]*      
![  METHOD: 2  IERR: 0  APPROXIMATE VALUE:  .38682227139652E+00 ]*      
![                                                              ]*      
![  ESTIMATED ERROR:-.3081E-10  ACTUAL ERROR:  .1465E-11        ]*      
![                                                              ]*      
![  NUMBER OF FUNCTION CALLS:    720                            ]*      
![                                                              ]*      
![                                                              ]*      
![  METHOD: 3  IERR: 0  APPROXIMATE VALUE:  .38682227139506E+00 ]*      
![                                                              ]*      
![  ESTIMATED ERROR: .2239E-14  ACTUAL ERROR:  .5551E-16        ]*      
![                                                              ]*      
![  NUMBER OF FUNCTION CALLS:   1280                            ]*      
![                                                              ]*      
![                                                              ]*      
![  METHOD: 4  IERR: 0  APPROXIMATE VALUE:  .38682227139506E+00 ]*      
![                                                              ]*      
![  ESTIMATED ERROR:-.1110E-15  ACTUAL ERROR: -.4996E-15        ]*      
![                                                              ]*      
![  NUMBER OF FUNCTION CALLS:   2000                            ]*      
![                                                              ]*      
![                                                              ]*      
![  METHOD: 5  IERR: 0  APPROXIMATE VALUE:  .38682227139506E+00 ]*      
![                                                              ]*      
![  ESTIMATED ERROR: .2035E-15  ACTUAL ERROR:  .4996E-15        ]*      
![                                                              ]*      
![  NUMBER OF FUNCTION CALLS:   2880                            ]*      
![                                                              ]*      
![                                                              ]*      
![  METHOD: 6  IERR: 0  APPROXIMATE VALUE:  .38682227139506E+00 ]*      
![                                                              ]*      
![  ESTIMATED ERROR:-.3331E-15  ACTUAL ERROR: -.3886E-15        ]*      
![                                                              ]*      
![  NUMBER OF FUNCTION CALLS:   3920                            ]*      
![                                                              ]*      
![                                                              ]*      
![  METHOD: 7  IERR: 0  APPROXIMATE VALUE:  .38682227139506E+00 ]*      
![                                                              ]*      
![  ESTIMATED ERROR: .2776E-15  ACTUAL ERROR: -.1110E-15        ]*      
![                                                              ]*      
![  NUMBER OF FUNCTION CALLS:   5120                            ]*      
![                                                              ]*      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: SXCY, K4GAUV, K4GINI                    *      
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
!                                                                       
!   number of rectangles in X nd Y directions                           
!                                                                       
      PARAMETER (NX = 4, NY = 4) 
      EXTERNAL SXCY 
      DOUBLEPRECISION WORK (2, 0:7), X (0:NX), Y (0:NY) 
      DOUBLEPRECISION A, B, C, D, DX, DY, CREC, DIVIAT, EXACT, SXCY 
      LOGICAL ESTDIV 
!                                                                       
! outer corners                                                         
!                                                                       
      A = 0.0D+00 
      B = 1.0D+00 
      C = 0.0D+00 
      D = 1.0D+00 
!                                                                       
! define X vector                                                       
!                                                                       
      X (0) = A 
      X (NX) = B 
      DX = (B - A) / DBLE (NX) 
      DO 10 I = 1, NX - 1 
         X (I) = X (I - 1) + DX 
   10 END DO 
!                                                                       
! define Y vector                                                       
!                                                                       
      Y (0) = C 
      Y (NY) = D 
      DY = (D-C) / DBLE (NY) 
      DO 20 I = 1, NY - 1 
         Y (I) = Y (I - 1) + DY 
   20 END DO 
!                                                                       
! initialize MOLD for first call of K4GAUV                              
!                                                                       
      MOLD = - 1 
!                                                                       
! swith error estimation                                                
!                                                                       
      ESTDIV = .TRUE. 
!                                                                       
! exact solution                                                        
!                                                                       
      EXACT = (COS (B) - COS (A) ) * (SIN (C) - SIN (D) ) 
      WRITE ( *, 1000) EXACT 
!                                                                       
! compute approximate solutions and error estimates                     
!                                                                       
      DO 30 METHOD = 0, 7 
         CALL K4GAUV (SXCY, X, NX, Y, NY, METHOD, MOLD, CREC, ESTDIV,   &
         DIVIAT, WORK, IERR, IUFCLL)                                    
         WRITE ( *, 1100) METHOD, IERR, CREC 
         WRITE ( *, 1200) DIVIAT, CREC - EXACT 
         WRITE ( *, 1300) IUFCLL 
   30 END DO 
      STOP 
!                                                                       
! Formatangaben                                                         
!                                                                       
 1000 FORMAT (1X,'C[',T66,']*',/,                                       &
     &        1X,'C[',2X,'EXACT SOLUTION:',E20.14,T66,']*',/,           &
     &        1X,'C[',T66,']*')                                         
 1100 FORMAT (1X,'C[',T66,']*',/,                                       &
     &        1X,'C[',2X,'METHOD: ',I1,2X,'IERR: ',I1,2X,               &
     &        'APPROXIMATE VALUE: ',E20.14,T66,']*')                    
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
