      SUBROUTINE SHPGLO (X, Y, FX, FY, F, W, R, N, DMUE, PHI) 
!********************************************************************   
!                                                                   *   
!   Program name: SHPGLO                                            *   
!                                                                   *   
!********************************************************************   
!                                                                   *   
!   This subroutine computes one functional value at (X,Y) for given*   
!   nodes using the global Shepard method.                          *   
!   The exponent dmue must be specified externally.                 *   
!   DMUE should be chosen to lie betwen 2 and 6.                    *   
!                                                                   *   
!********************************************************************   
!                                                                   *   
!   Input parameters:                                               *   
!   =================                                               *   
!   X    :  X value for which we want to interpolate the Z value    *   
!   Y    :  Y value for which we want to interpolate the Z value    *   
!   FX, FY, F   :  vectors ..(0:N) with X and Y coordinates of nodes*   
!                  (FX,FY) and corresponding functional value F.    *   
!   N    :  Index of last node                                      *   
!   DMUE :  Exponent, 0 < DMUE < infinity, reasonable results can   *   
!           be achieved for 2 < DMUE < 6. If on input DMUE <= 0, we *   
!           set DMUE = 2 internally.                                *   
!                                                                   *   
!   AUX VECTORS:                                                    *   
!   ============                                                    *   
!   W, R :  vectors ..(0:N)                                         *   
!                                                                   *   
!                                                                   *   
!   Output parameters:                                              *   
!   ==================                                              *   
!   PHI  : Interpolated Z value at (X,Y)                            *   
!                                                                   *   
!********************************************************************   
!                                                                   *   
!   Required subroutines: none                                      *   
!                                                                   *   
!********************************************************************   
!                                                                   *   
!   Author      : Bjoern Terwege                                    *   
!   Date        : 6.12.1995                                         *   
!   Source code : FORTRAN 77                                        *   
!                                                                       
!********************************************************************   
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION R (0:N), W (0:N), FX (0:N), FY (0:N), F (0:N) 
!                                                                       
!  Check that DMUE > 0, otherwise set DMUE = 2                          
!                                                                       
      IF (DMUE.LE.0.) DMUE = 2 
!                                                                       
!  Compute the R(I)                                                     
!                                                                       
      DO 20 I = 0, N 
         R (I) = DSQRT ( ( (X - FX (I) ) * (X - FX (I) ) ) + ( (Y - FY (&
         I) ) * (Y - FY (I) ) ) )                                       
         IF (R (I) .EQ.0) THEN 
            PHI = F (I) 
            GOTO 333 
         ENDIF 
   20 END DO 
!                                                                       
!  Compute the sum in the denominator of W(I)                           
!                                                                       
      DO 30 J = 0, N 
         SUM = 0 
         DO 40 I = 0, N 
            SUM = SUM + 1 / (R (I) **DMUE) 
   40    END DO 
!       Berechnung der Gewichte W(I)                                    
!                                                                       
         W (J) = 1 / (R (J) **DMUE) * 1 / SUM 
   30 END DO 
      PHI = 0 
!                                                                       
!  Compute the approximate function value at (X,Y)                      
!                                                                       
      DO 50 J = 0, N 
         PHI = PHI + W (J) * F (J) 
   50 END DO 
  333 CONTINUE 
      RETURN 
      END SUBROUTINE SHPGLO                         
