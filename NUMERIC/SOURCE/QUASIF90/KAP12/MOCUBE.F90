![  {Modified (Interpolating) Cubic B\'ezier Splines}                   
![  {Modified (Interpolating) Cubic B\'ezier Splines}*)                 
      SUBROUTINE MOCUBE (D, B, M, EPS) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     PROGRAM OBJECTIVE:                                         *      
!     ==================                                         *      
!     SUBROUTINE MOCUBE determines the coefficients of a         *      
!     modified BEZIER spline from the weight points so that the  *      
!     weight points will be located on the computed curve, except*      
!     for an accuracy of EPS.                                    *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     M                : number of curve segments of the spline  *      
!     EPS              : accuracy bound for the interpolation    *      
!     D(1,J), J=0,...,M: X-coordinates of the weight points      *      
!     D(2,J), M=0,...,M: Y-coordinates of the weight points      *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     B(1,J), J=0,1,...,3*M: X-coordinates of the BEZIER points  *      
!     B(2,J), J=0.1,...,3*M: Y-coordinates of the BEZIER points  *      
!                                                                *      
!                                                                *      
!     LOCAL VARIABLES:                                           *      
!     ================                                           *      
!     DINT:  difference of a weight points and the corresponding *      
!            BEZIER point                                        *      
!     J,K:   control variables                                   *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Juergen Dietel                                     *      
!  date     : 04.23.1987                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      INTEGER M, J, K 
      DOUBLEPRECISION D (2, 0:M), B (2, 0:3 * M), EPS, DINT 
!                                                                       
!*****************************************************************      
!     determine the BEZIER points                                *      
!     B(1,J), J=0,1,...,3*M (X-direction) and                    *      
!     B(2,J), J=0,1,...,3*M (Y-direction)                        *      
!*****************************************************************      
!                                                                       
      DO 40 K = 1, 2 
         DO 30 J = 1, M - 1 
            B (K, 3 * J - 2) = (2.0D0 * D (K, J - 1) + D (K, J) )       &
            / 3.0D0                                                     
            B (K, 3 * J) = (D (K, J - 1) + 4.0D0 * D (K, J) + D (K, J + &
            1) ) / 6.0D0                                                
            B (K, 3 * J + 2) = (D (K, J) + 2.0D0 * D (K, J + 1) )       &
            / 3.0D0                                                     
   30    END DO 
         B (K, 0) = D (K, 0) 
         B (K, 2) = (D (K, 0) + 2.0D0 * D (K, 1) ) / 3.0D0 
         B (K, 3 * M - 2) = (2.0D0 * D (K, M - 1) + D (K, M) ) / 3.0D0 
         B (K, 3 * M) = D (K, M) 
   40 END DO 
!                                                                       
!*****************************************************************      
!     correction of the BEZIER points                            *      
!*****************************************************************      
!                                                                       
   50 DO 90 K = 1, 2 
         DO 80 J = 3, 3 * M - 3, 3 
            DINT = D (K, J / 3) - B (K, J) 
            IF (J.NE.3) B (K, J - 3) = B (K, J - 3) + DINT / 4.0D0 
            B (K, J - 2) = B (K, J - 2) + DINT / 2.0D0 
            B (K, J - 1) = B (K, J - 1) + DINT 
            B (K, J) = B (K, J) + DINT 
            B (K, J + 1) = B (K, J + 1) + DINT 
            B (K, J + 2) = B (K, J + 2) + DINT / 2.0D0 
            IF (J.NE.3 * M - 3) B (K, J + 3) = B (K, J + 3) + DINT /    &
            4.0D0                                                       
   80    END DO 
   90 END DO 
!                                                                       
!*****************************************************************      
!     check whether interpolation has reached an accuracy of EPS,*      
!     otherwise repeat correcting the BEZIER points              *      
!*****************************************************************      
!                                                                       
      DO 100 J = 1, M - 1 
         IF (DABS (D (1, J) - B (1, 3 * J) ) + DABS (D (2, J) - B (2, 3 &
         * J) ) .GT.EPS) GOTO 50                                        
  100 END DO 
      RETURN 
      END SUBROUTINE MOCUBE                         
