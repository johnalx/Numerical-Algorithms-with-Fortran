      PROGRAM TEST 
!                                                 12.12.1991  ( Dubois G
!***********************************************************************
!                                                                       
!   Testprogram for the subroutine  ZERORF  for finding the roots of a  
!   continuous real valued funstion.                                    
!                                                                       
![  FUNCTION:  F(X) = EXP(X)-2.*X-3.                                    
![                                                                      
![  STARTING VALUES:  X1 =  -2.00  AND  X2 =  -2.00          IMAX =  99 
![                                                                      
![  FUNCTIONAL VALUES:  F(X1) =  .11353352832366D+01                    
![                      F(X2) =  .11353352832366D+01     IEIN = 0       
![                                                                      
![                                                                      
![  IEXTRA = 0   ABSERR =  .5000D-14   IMETH = 1   IZAEHL =  7   IFEHL =
![  X(1) = -.13733746911009D+01                                         
![  X(2) = -.13733745453522D+01     F(2) =  .45097519468795D-12         
![  X(3) = -.13733745453519D+01     F(3) = -.10169816377914D-15         
![                                                                      
![  IEXTRA = 0   ABSERR =  .5000D-14   IMETH = 2   IZAEHL =  3   IFEHL =
![  X(1) = -.13733746911009D+01                                         
![  X(2) = -.13733745453522D+01     F(2) =  .45097519468795D-12         
![  X(3) = -.13733745453519D+01     F(3) = -.10169816377914D-15         
![                                                                      
![  IEXTRA = 0   ABSERR =  .5000D-14   IMETH = 3   IZAEHL =  3   IFEHL =
![  X(1) = -.13733746911009D+01                                         
![  X(2) = -.13733745453522D+01     F(2) =  .45097519468795D-12         
![  X(3) = -.13733745453519D+01     F(3) = -.10169816377914D-15         
![                                                                      
![  IEXTRA = 0   ABSERR =  .5000D-14   IMETH = 4   IZAEHL =  3   IFEHL =
![  X(1) = -.13733746911009D+01                                         
![  X(2) = -.13733745453522D+01     F(2) =  .45097519468795D-12         
![  X(3) = -.13733745453519D+01     F(3) = -.10169816377914D-15         
![                                                                      
![                                                                      
![  IEXTRA = 1   ABSERR =  .5000D-14   IMETH = 1   IZAEHL =  3   IFEHL =
![  X(1) = -.13733746911009D+01                                         
![  X(2) = -.13733745453522D+01     F(2) =  .45097519468795D-12         
![  X(3) = -.13733745453519D+01     F(3) = -.10169816377914D-15         
![                                                                      
![  IEXTRA = 1   ABSERR =  .5000D-14   IMETH = 2   IZAEHL =  3   IFEHL =
![  X(1) = -.13733746911009D+01                                         
![  X(2) = -.13733745453522D+01     F(2) =  .45097519468795D-12         
![  X(3) = -.13733745453519D+01     F(3) = -.10169816377914D-15         
![                                                                      
![  IEXTRA = 1   ABSERR =  .5000D-14   IMETH = 3   IZAEHL =  3   IFEHL =
![  X(1) = -.13733746911009D+01                                         
![  X(2) = -.13733745453522D+01     F(2) =  .45097519468795D-12         
![  X(3) = -.13733745453519D+01     F(3) = -.10169816377914D-15         
![                                                                      
![  IEXTRA = 1   ABSERR =  .5000D-14   IMETH = 4   IZAEHL =  3   IFEHL =
![  X(1) = -.13733746911009D+01                                         
![  X(2) = -.13733745453522D+01     F(2) =  .45097519468795D-12         
![  X(3) = -.13733745453519D+01     F(3) = -.10169816377914D-15         
![                                                                      
![                                                                      
![  FOR  IEXTRA = 0  AND  METHOD 1:  SUM(IZAEH0) =   7                  
![  FOR  IEXTRA = 0  AND  METHOD 2:  SUM(IZAEH0) =   3                  
![  FOR  IEXTRA = 0  AND  METHOD 3:  SUM(IZAEH0) =   3                  
![  FOR  IEXTRA = 0  AND  METHOD 4:  SUM(IZAEH0) =   3                  
![                                                                      
![  FOR  IEXTRA = 1  AND  METHOD 1:  SUM(IZAEH1) =   3                  
![  FOR  IEXTRA = 1  AND  METHOD 2:  SUM(IZAEH1) =   3                  
![  FOR  IEXTRA = 1  AND  METHOD 3:  SUM(IZAEH1) =   3                  
![  FOR  IEXTRA = 1  AND  METHOD 4:  SUM(IZAEH1) =   3                  
!                                                                       
!   The above results were computed on a PC and are displayed on screen;
!   can also be stored in a file.                                       
!                                                                       
!***********************************************************************
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION FKT, DELX, ABSERR, RELERR, FAB, X1, X2 
      INTEGER IMAX, IZAEH0 (4), IZAEH1 (4), IEIN 
      CHARACTER(14) FVONX 
      EXTERNAL FKT 
!                                                                       
!  Initialize starting values for the function                          
!                                                                       
      DATA DELX / 1.D-8 / 
      DATA IMAX / 99 / 
      DATA ABSERR / .5D-14 / 
      DATA RELERR / 0.D0 / 
      DATA FAB / .5D-14 / 
      DATA X1 / - 2.D0 / 
      DATA X2 / - 2.D0 / 
      DATA FVONX / 'EXP(X)-2.*X-3.' / 
      DATA IZAEH0 / 4 * 0 / 
      DATA IZAEH1 / 4 * 0 / 
!                                                                       
!  Put out and check actual starting values                             
!                                                                       
      WRITE ( *, 900) FVONX 
      WRITE ( *, 910) X1, X2, IMAX 
      IEIN = 0 
      IF (FKT (X1) * FKT (X2) .LT.0.D0) IEIN = 1 
      WRITE ( *, 920) FKT (X1), FKT (X2), IEIN 
!                                                                       
!  Call ZERORF                                                          
!                                                                       
      DO 10 IEXTRA = 0, 1 
         DO 20 IMETH = 1, 4 
!                                                                       
            CALL ZERORF (FKT, ABSERR, RELERR, FAB, IMAX, IMETH, IEXTRA, &
            DELX, X1, X2, X3, F1, F2, F3, IZAEHL, IHILF1, IHILF2,       &
            IEINST, IFEHL)                                              
            IF (IEXTRA.EQ.0) THEN 
               IZAEH0 (IMETH) = IZAEH0 (IMETH) + IZAEHL 
            ELSE 
               IZAEH1 (IMETH) = IZAEH1 (IMETH) + IZAEHL 
            ENDIF 
!                                                                       
!  Output of test results                                               
!                                                                       
            WRITE ( *, 930) IEXTRA, ABSERR, IMETH, IZAEHL, IFEHL 
            WRITE ( *, 940) X1, X2, F2, X3, F3 
   20    END DO 
         WRITE ( *, 980) 
   10 END DO 
      WRITE ( *, 950) (I, IZAEH0 (I), I = 1, 4) 
      WRITE ( *, 980) 
      WRITE ( *, 960) (I, IZAEH1 (I), I = 1, 4) 
      STOP 
!                                                                       
  900 FORMAT(1X,'C[  FUNCTION:  F(X) = ',A,T78,']*',/,                  &
     &       1X,'C[',T78,']*')                                          
  910 FORMAT(1X,'C[  STARTING VALUES:  X1 = ',F6.2,'  AND  X2 = ',F6.2, &
     &       10X,'IMAX = ',I3,T78,']*',/,                               &
     &       1X,'C[',T78,']*')                                          
  920 FORMAT(1X,'C[  FUNCTIONAL VALUES:',T26,'F(X1) = ',D20.14,         &
     &          T78,']*',/,                                             &
     &       1X,'C[',T26,'F(X2) = ',D20.14,5X,'IEIN = ',I1,T78,']*',/,  &
     &       1X,'C[',T78,']*',/,                                        &
     &       1X,'C[',T78,']*')                                          
  930 FORMAT(1X,'C[  IEXTRA = ',I1,3X,'ABSERR = ',D10.4,3X,'IMETH = ',  &
     &       I1,3X,'IZAEHL = ',I2,3X,'IFEHL = ',I1,T78,']*')            
  940 FORMAT(1X,'C[  X(1) = ',D20.14,T78,']*',/,                        &
     &       1X,'C[  X(2) = ',D20.14,5X,'F(2) = ',D20.14,T78,']*',/,    &
     &       1X,'C[  X(3) = ',D20.14,5X,'F(3) = ',D20.14,T78,']*',/,    &
     &       1X,'C[',T78,']*')                                          
  950 FORMAT(1X,'C[  FOR  IEXTRA = 0  AND  METHOD ',                    &
     &       I1,':  SUM(IZAEH0) = ',I3,T78,']*')                        
  960 FORMAT(1X,'C[  FOR  IEXTRA = 1  AND  METHOD ',                    &
     &       I1,':  SUM(IZAEH1) = ',I3,T78,']*')                        
  980 FORMAT(1X,'C[',T78,']*') 
      END PROGRAM TEST                              
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION FKT (X) 
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      FKT = DEXP (X) - 2.D0 * X - 3.D0 
      RETURN 
      END FUNCTION FKT                              
