      SUBROUTINE CUBBEZ (B, D, M) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     This subroutine determines BEZIER points of a curve by     *      
!     the cubic BEZIER-method.                                   *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     M                  number of curve segments                *      
!     D(J,K)             coordinates of the weight points        *      
!                        J=1, ..., 3, K=0, ..., M                *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETER S                                         *      
!     ==================                                         *      
!     B(J,K)             coordinates of the BEZIER points        *      
!                        J=1, ..., 3, K=0, ..., 3*M              *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Michael Radermacher                                *      
!  date     : 04.30.1985                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION B (3, 0:3 * M), D (3, 0:M) 
!                                                                       
!*****************************************************************      
!     loop over the X- , Y- and Z-coordinates                    *      
!*****************************************************************      
!                                                                       
      DO 20 K = 1, 3 
         DO 10 J = 1, M - 1 
            B (K, 3 * J - 2) = (2.0D0 * D (K, J - 1) + D (K, J) )       &
            / 3.0D0                                                     
            B (K, 3 * J) = (D (K, J - 1) + 4.0D0 * D (K, J) + D (K, J + &
            1) ) / 6.0D0                                                
   10    B (K, 3 * J + 2) = (D (K, J) + 2.0D0 * D (K, J + 1) ) / 3.0D0 
         B (K, 2) = (D (K, 0) + 2.0D0 * D (K, 1) ) / 3.0D0 
         B (K, 3 * M - 2) = (D (K, M) + 2.0D0 * D (K, M - 1) ) / 3.0D0 
!                                                                       
!*****************************************************************      
!        the boundary points B(K,0) and B(K,3*M), K=1,...,3, are *      
!        preset so that a natural cubic BEZIER spline will result*      
!*****************************************************************      
!                                                                       
         B (K, 0) = D (K, 0) 
   20 B (K, 3 * M) = D (K, M) 
      RETURN 
      END SUBROUTINE CUBBEZ                         
