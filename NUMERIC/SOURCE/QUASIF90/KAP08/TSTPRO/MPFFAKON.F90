      PROGRAM TESFKN 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This program tests the subroutine  FFAKON  with values of the *      
!  functions                                                     *      
!       SIN(X)       for  X in [0, PI]                           *      
!        0           otherwise                                   *      
!  and                                                           *      
!        0.5         for  X = 0                                  *      
!       (PI-X)/PI    for  X in (0, PI]                           *      
!        0           otherwise,                                  *      
!  whose (nonperiodic) convolution is given by                   *      
!       1 - X/PI - COS(X) + SIN(X)/PI   for  X  in [0, PI]       *      
!       2 - X/PI + SIN(X)/PI            for  X  in (PI, 2*PI]    *      
!       0                               otherwise,               *      
!  while their (nonperiodic) correlation is                      *      
!       1 + X/PI + SIN(X)/PI            for  X  in [-PI, 0]      *      
!       X/PI + COS(X) + SIN(X)/PI       for  X  in (0, PI]       *      
!       0                               otherwise.               *      
!  We use equidistant nodes. An increase in the number of nodes  *      
!  must increase the accuracy of the approximation.              *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      INTEGER M (8), INDEX, L (8), IFEHL 
      DIMENSION F1RE (0:1023), F1IM (0:1023), F2RE (0:1023), F2IM (0:   &
      1023)                                                             
      DATA M / 10, 10, 41, 41, 128, 128, 500, 500 / 
      DATA L / 32, 32, 128, 128, 512, 512, 1024, 1024 / 
      INDEX = 1 
!                                                                       
!     Assign function values                                            
!                                                                       
      PI = 4.0D0 * ATAN (1.0D0) 
      DO 50 I = 1, 8 
         DX = PI / DBLE (M (I) ) 
         DO 10 J = 0, M (I) 
            XJ = DBLE (J) * DX 
            F1RE (J) = SIN (XJ) 
            F1IM (J) = 0.0D0 
            F2RE (J) = (PI - XJ) / PI 
            F2IM (J) = 0.0D0 
   10    END DO 
         F2RE (0) = 0.5D0 
!                                                                       
!     Compute convolution (INDEX=0) or correlation (INDEX=1)            
!                                                                       
         INDEX = 1 - INDEX 
         CALL FFAKON (M (I), F1RE, F1IM, M (I), F2RE, F2IM, L (I),      &
         DX, INDEX, IFEHL)                                              
!                                                                       
!     Compare with theoretical output                                   
!                                                                       
         WRITE (6, * ) 'Error parameter :', IFEHL 
         FFALT = 0.0D0 
         FFALTI = 0.0D0 
         FKORR = 0.0D0 
         FKORRI = 0.0D0 
         IF (INDEX.EQ.0) THEN 
            DO 20 J = 0, 2 * M (I) 
               XJ = DBLE (J) * DX 
               IF (J.LE.M (I) ) THEN 
                  FUFAL = 1.0D0 - XJ / PI - COS (XJ) + SIN (XJ) / PI 
               ELSE 
                  FUFAL = 2.0D0 - XJ / PI + SIN (XJ) / PI 
               ENDIF 
               FFALT = MAX (FFALT, ABS (F1RE (J) - FUFAL) ) 
               FFALTI = MAX (FFALTI, ABS (F1IM (J) ) ) 
   20       END DO 
      WRITE (6, 900) M (I) , 'Convolution    ', FFALT, FFALTI 
         ELSE 
            DO 30 J = 0, 2 * M (I) 
               XJ = - PI + DBLE (J) * DX 
               IF (J.LE.M (I) ) THEN 
                  FUKOR = 1.0D0 + XJ / PI + SIN (XJ) / PI 
               ELSE 
                  FUKOR = XJ / PI + COS (XJ) + SIN (XJ) / PI 
               ENDIF 
               FKORR = MAX (FKORR, ABS (F1RE (J) - FUKOR) ) 
               FKORRI = MAX (FKORRI, ABS (F1IM (J) ) ) 
   30       END DO 
            WRITE (6, 900) M (I) , 'Correlation', FKORR, FKORRI 
         ENDIF 
   50 END DO 
      STOP 
  900 FORMAT(1X,'M = N = ',I3,'  Error at ',A,' : ',2(1PD10.3,2X)) 
      END PROGRAM TESFKN                            
