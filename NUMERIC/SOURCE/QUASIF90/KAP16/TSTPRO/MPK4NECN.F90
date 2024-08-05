      PROGRAM TEST 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Test program for cubature over rectangular region using       *      
!  Newton-Cotes formulas:                                        *      
!                                                                *      
!  We test subroutine  K4NECN                                    *      
!                                                                *      
!  The test example produces the output:                         *      
!                                                                *      
![                                                              ]*      
![  EXACT SOLUTION: .38682227139506E+00                         ]*      
![                                                              ]*      
![                                                              ]*      
![  METHOD: 1  IERR: 0  APPROXIMATE VALUE:  .38581531547004E+00 ]*      
![                                                              ]*      
![  ESTIMATED ERROR: .1005E-02  ACTUAL ERROR: -.1007E-02        ]*      
![                                                              ]*      
![  NUMBER OF FUNCTION CALLS:    106                            ]*      
![                                                              ]*      
![                                                              ]*      
![  METHOD: 2  IERR: 0  APPROXIMATE VALUE:  .38682233700823E+00 ]*      
![                                                              ]*      
![  ESTIMATED ERROR:-.6571E-07  ACTUAL ERROR:  .6561E-07        ]*      
![                                                              ]*      
![  NUMBER OF FUNCTION CALLS:    370                            ]*      
![                                                              ]*      
![                                                              ]*      
![  METHOD: 3  IERR: 0  APPROXIMATE VALUE:  .38682230055496E+00 ]*      
![                                                              ]*      
![  ESTIMATED ERROR:-.2920E-07  ACTUAL ERROR:  .2916E-07        ]*      
![                                                              ]*      
![  NUMBER OF FUNCTION CALLS:    794                            ]*      
![                                                              ]*      
![                                                              ]*      
![  METHOD: 4  IERR: 0  APPROXIMATE VALUE:  .38682227139353E+00 ]*      
![                                                              ]*      
![  ESTIMATED ERROR: .1528E-11  ACTUAL ERROR: -.1525E-11        ]*      
![                                                              ]*      
![  NUMBER OF FUNCTION CALLS:   1378                            ]*      
![                                                              ]*      
![                                                              ]*      
![  METHOD: 5  IERR: 0  APPROXIMATE VALUE:  .38682227139420E+00 ]*      
![                                                              ]*      
![  ESTIMATED ERROR: .8606E-12  ACTUAL ERROR: -.8588E-12        ]*      
![                                                              ]*      
![  NUMBER OF FUNCTION CALLS:   2122                            ]*      
![                                                              ]*      
![                                                              ]*      
![  METHOD: 6  IERR: 0  APPROXIMATE VALUE:  .38682227139506E+00 ]*      
![                                                              ]*      
![  ESTIMATED ERROR:-.3026E-16  ACTUAL ERROR: -.3331E-15        ]*      
![                                                              ]*      
![  NUMBER OF FUNCTION CALLS:   3026                            ]*      
![                                                              ]*      
![                                                              ]*      
![  METHOD: 7  IERR: 0  APPROXIMATE VALUE:  .38682227139506E+00 ]*      
![                                                              ]*      
![  ESTIMATED ERROR:-.1611E-16  ACTUAL ERROR:  .3331E-15        ]*      
![                                                              ]*      
![  NUMBER OF FUNCTION CALLS:   4090                            ]*      
![                                                              ]*      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  required subroutines : SXCY, K4NECN, K4INIT, GRIDOT           *      
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
      DOUBLEPRECISION WORK (0:9) 
      DOUBLEPRECISION A, B, C, D, CREC, DIVIAT, EXACT, SXCY 
      LOGICAL ESTDIV 
!                                                                       
! number of rectangles in X and Y direction                             
!                                                                       
      IP = 4 
      IQ = 4 
!                                                                       
! outer corners of rectangle                                            
!                                                                       
      A = 0.0D+00 
      B = 1.0D+00 
      C = 0.0D+00 
      D = 1.0D+00 
!                                                                       
! set MOLD for first call of  K4NECN                                    
!                                                                       
      MOLD = 8 
!                                                                       
! with error estimation                                                 
!                                                                       
      ESTDIV = .TRUE. 
!                                                                       
! exact solution                                                        
!                                                                       
      EXACT = (COS (B) - COS (A) ) * (SIN (C) - SIN (D) ) 
      WRITE ( *, 1000) EXACT 
!                                                                       
! compute approximate solutions for all methods and                     
! estimate the error                                                    
!                                                                       
      DO 10 METHOD = 1, 7 
         CALL K4NECN (SXCY, A, B, IP, C, D, IQ, METHOD, MOLD, CREC,     &
         ESTDIV, DIVIAT, WORK, IERR, IUFCLL)                            
         WRITE ( *, 1100) METHOD, IERR, CREC 
         WRITE ( *, 1200) DIVIAT, CREC - EXACT 
         WRITE ( *, 1300) IUFCLL 
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
