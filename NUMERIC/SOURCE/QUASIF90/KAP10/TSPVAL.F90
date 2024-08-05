      DOUBLEPRECISION FUNCTION TSPVAL (PHI, N, PHIN, A, B, C, D) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  The FUNCTION TSPVAL evaluates the value of a transformed      *      
!  parametric cubic spline given in the form:                    *      
!                                                                *      
!  S(PHI) = A(I) + B(I)(PHI-PHIN(I)) + C(I)(PHI-PHIN(I))**2 +    *      
!                                    + D(I)(PHI-PHIN(I))**3      *      
!                                                                *      
!  for PHI in the interval [PHIN(I),PHIN(I+1)], I=0,1,...,N-1.   *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  PHI  :  value where the spline shall be evaluated (in radians)*      
!  N    :  index of the lst node                                 *      
!  PHIN :  vector PHIN(0:N); the nodes PHIN(I), I=0,1,...,N      *      
!  A    :  ] N+1-vectors ..(0:N);                                *      
!  B    :  ] the elements in positions 0 to N-1 denote the       *      
!  C    :  ] coefficients of the spline.                         *      
!  D    :  ]                                                     *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Subroutines required: SPVAL                                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : GÅnter Palm                                        *      
!  Date     : 06.01.1991                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
!-----Declarations                                                      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION PHIN (0:N), A (0:N), B (0:N), C (0:N), D (0:N) 
!                                                                       
!-----Initialize                                                        
!                                                                       
      ZWOPI = 8.0D0 * DATAN (1.0D0) 
!                                                                       
!-----Assign PHI to PHIX,                                               
!     if necessary, recompute PHIX so that it lies                      
!     in the interval [0,2*PI]                                          
!                                                                       
      IF (PHI.LT.0.0D0) THEN 
         L = ABS (INT (PHI / ZWOPI) ) + 1 
         PHIX = L * ZWOPI - PHI 
      ELSEIF (PHI.GT.ZWOPI) THEN 
         L = INT (PHI / ZWOPI) 
         PHIX = PHI - L * ZWOPI 
      ELSE 
         PHIX = PHI 
      ENDIF 
!                                                                       
!-----Compute the functional value and that of the derivatives          
!     at PHIX in SUBROUTINE SPVAL                                       
!                                                                       
      TSPVAL = SPVAL (PHIX, N, PHIN, A, B, C, D) 
!                                                                       
      RETURN 
      END FUNCTION TSPVAL                           
