      SUBROUTINE SHPLOK (X, Y, FX, FY, F, PPHI, W, R, N, DMUE, RR, PHI) 
!********************************************************************   
!                                                                   *   
!   Program name: SHPLOK                                            *   
!                                                                   *   
!********************************************************************   
!                                                                   *   
!   This subroutine computes one functional value at (X,Y) for given*   
!   nodes using the local Shepard method.                           *   
!                                                                   *   
!********************************************************************   
!                                                                   *   
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
!   RR   :  Radius around (X,Y) inside which all nodes are used to  *   
!           interpolate at (X,Y).                                   *   
!                                                                   *   
!                                                                   *   
!   AUX VECTORS:                                                    *   
!   ============                                                    *   
!   W, R, PPHI :  vectors ..(0:N)                                   *   
!                                                                   *   
!                                                                   *   
!   Output parameters:                                              *   
!   ==================                                              *   
!   PHI  : Interpolated Z value at (X,Y)                            *   
!                                                                   *   
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
!                                                                   *   
!********************************************************************   
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION R (0:N), W (0:N), FX (0:N), FY (0:N), F (0:N), PPHI (0: &
      N)                                                                
      IF (DMUE.LT.0.) DMUE = 2. 
      DO 20 I = 0, N 
         R (I) = DSQRT ( ( (X - FX (I) ) * (X - FX (I) ) ) + ( (Y - FY (&
         I) ) * (Y - FY (I) ) ) )                                       
         IF (R (I) .EQ.0) THEN 
            PHI = F (I) 
            GOTO 111 
         ENDIF 
   20 END DO 
!                                                                       
!  Compute the PPHI(I)                                                  
!                                                                       
      DO 40 I = 0, N 
         IF (R (I) .GE.RR) THEN 
            PPHI (I) = 0 
         ELSE 
            PPHI (I) = (RR / R (I) ) - 1 
         ENDIF 
   40 END DO 
!                                                                       
!  Compute the denominators of the weights                              
!                                                                       
      SUM = 0 
      DO 50 J = 0, N 
         IF (PPHI (J) .EQ.0) GOTO 50 
         SUM = SUM + 1 / (PPHI (J) **DMUE) 
   50 END DO 
!                                                                       
!  Compute the weights                                                  
!                                                                       
      DO 60 J = 0, N 
         IF (PPHI (J) .EQ.0) THEN 
            W (J) = 0 
            GOTO 60 
         ELSE 
            W (J) = 1 / (PPHI (J) **DMUE) * 1 / SUM 
         ENDIF 
   60 END DO 
!                                                                       
!  Compute the approximate function value Z at (X,Y)                    
!                                                                       
      PHI = 0 
      DO 70 I = 0, N 
         PHI = PHI + W (I) * F (I) 
   70 END DO 
  111 CONTINUE 
      RETURN 
      END SUBROUTINE SHPLOK                         
