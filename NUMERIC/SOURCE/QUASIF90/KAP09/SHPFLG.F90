      SUBROUTINE SHPFLG (X, Y, FX, FY, F, W, R, EPS, N, DMUE, RR, PHI,  &
      IERR)                                                             
!********************************************************************   
!                                                                   *   
!   Program name: SHPFLG                                            *   
!                                                                   *   
!********************************************************************   
!                                                                   *   
!   This subroutine computes one functional value at (X,Y) for given*   
!   nodes using the local Shepard method and Franke-Little weights. *   
!   The exponent dmue and the radius rr must be specified externally*   
!   DMUE should be chosen to lie betwen 2 and 6.                    *   
!   The radius RR should be so that the circle of radius RR around  *   
!   (X,Y) contains some nodes.                                      *   
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
!   RR   :  Radius around (X,Y) inside which all nodes are used to  *   
!           interpolate at (X,Y).                                   *   
!                                                                   *   
!                                                                   *   
!   AUX VECTORS:                                                    *   
!   ============                                                    *   
!   W, R, EPS   : vectors ..(0:N)                                   *   
!                                                                   *   
!                                                                   *   
!   Output parameters:                                              *   
!   ==================                                              *   
!   PHI  : Interpolated Z value at (X,Y)                            *   
!   IERR : Error parameter:                                         *   
!          = 0 : all correct                                        *   
!          = 1 : Sum of Franke-Little weights is zero               *   
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
      DIMENSION R (0:N), W (0:N), FX (0:N), FY (0:N), F (0:N), EPS (0:N) 
      IERR = 0 
!                                                                       
!  Check that DMUE > 0, otherwise set DMUE = 2                          
!                                                                       
      IF (DMUE.LE.0) DMUE = 2 
!                                                                       
!  Check that  RR > 0, otherwise set RR = 0.1                           
!                                                                       
      IF (RR.LE.0) RR = 0.1 
!                                                                       
!  Compute the R(I)                                                     
!                                                                       
      DO 20 I = 0, N 
         R (I) = DSQRT ( ( (X - FX (I) ) * (X - FX (I) ) ) + ( (Y - FY (&
         I) ) * (Y - FY (I) ) ) )                                       
         IF (R (I) .EQ.0) THEN 
            PHI = F (I) 
            GOTO 111 
         ENDIF 
   20 END DO 
!                                                                       
!  Compute the EPS(I)                                                   
!                                                                       
      DO 40 I = 0, N 
         IF (R (I) .GE.RR) THEN 
            EPS (I) = 0 
         ELSE 
            EPS (I) = 1 - (R (I) / RR) 
         ENDIF 
   40 END DO 
!                                                                       
!  Compute the numerators needed for the weights                        
!                                                                       
      SUM = 0 
      DO 50 I = 0, N 
         SUM = SUM + EPS (I) **DMUE 
   50 END DO 
      IF (SUM.EQ.0) THEN 
         IERR = 1 
         GOTO 111 
      ENDIF 
!                                                                       
!  Compute the weights                                                  
!                                                                       
      DO 60 J = 0, N 
         W (J) = (EPS (J) **DMUE) / SUM 
   60 END DO 
!                                                                       
!  Compute the approximate function value at (X,Y)                      
!                                                                       
      PHI = 0 
      DO 70 I = 0, N 
         PHI = PHI + W (I) * F (I) 
   70 END DO 
  111 CONTINUE 
      RETURN 
      END SUBROUTINE SHPFLG                         
