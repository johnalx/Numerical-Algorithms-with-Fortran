      PROGRAM TESFFK 
!                                                                       
!****************************************************************       
!                                                               *       
!  This program performs a test of the subroutine FFAKO  using  *       
!  values for the PI-periodic functions                         *       
!       | SIN X |                                               *       
!  and                                                          *       
!        0.5         for  X = 0                                 *       
!        (PI-X)/PI   for  0 < X < PI  ,                         *       
!  whose cyclic convolution is                                  *       
!       2*(PI-X)/PI**2 - COS(X)/PI ,                            *       
!  and whose cyclic correlation is                              *       
!       2*X/PI**2 + COS(X)/PI .                                 *       
!  We use equidistant nodes in the period interval [0, PI).     *       
!  Note that the matching accuracy must improve for increasing  *       
!  numbers of nodes.                                            *       
!                                                               *       
!****************************************************************       
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      INTEGER M (5) 
      DIMENSION F1RE (0:511), F1IM (0:511), F2RE (0:511), F2IM (0:511) 
      DATA M / 16, 41, 64, 188, 512 / 
!                                                                       
!     Set up real and imaginary parts of function values                
!                                                                       
      PI = 4.0D0 * ATAN (1.0D0) 
      DO 50 I = 1, 5 
         FAKT = PI / DBLE (M (I) ) 
         DO 10 J = 0, M (I) - 1 
            XJ = DBLE (J) * FAKT 
            F1RE (J) = ABS (SIN (XJ) ) 
            F1IM (J) = 0.0D0 
            F2RE (J) = (PI - XJ) / PI 
            F2IM (J) = 0.0D0 
   10    END DO 
         F2RE (0) = 0.5D0 
!                                                                       
!     Compute discrete cyclic convolution and correlation               
!                                                                       
         CALL FFAKO (M (I), F1RE, F1IM, F2RE, F2IM) 
!                                                                       
!     Compare with theoretically known exact results                    
!                                                                       
         FFALT = 0.0D0 
         FFALTI = 0.0D0 
         FKORR = 0.0D0 
         FKORRI = 0.0D0 
         DO 20 J = 0, M (I) - 1 
            XJ = DBLE (J) * FAKT 
            FUFAL = 2.0D0 * (PI - XJ) / PI**2 - COS (XJ) / PI 
            FUKOR = 2.0D0 * XJ / PI**2 + COS (XJ) / PI 
            FFALT = MAX (FFALT, ABS (F1RE (J) - FUFAL) ) 
            FFALTI = MAX (FFALTI, ABS (F1IM (J) ) ) 
            FKORR = MAX (FKORR, ABS (F2RE (J) - FUKOR) ) 
            FKORRI = MAX (FKORRI, ABS (F2IM (J) ) ) 
   20    END DO 
         WRITE (6, * ) 
      WRITE (6, 900) M (I) , 'Convolution  ', FFALT, FFALTI 
      WRITE (6, 900) M (I) , 'Correlation  ', FKORR, FKORRI 
   50 END DO 
      STOP 
  900 FORMAT(1X,'M = ',I3,'  Error at   ',A,' : ',2(1PD10.3,2X)) 
      END PROGRAM TESFFK                            
