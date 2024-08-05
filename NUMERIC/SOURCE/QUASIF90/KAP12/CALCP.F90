      SUBROUTINE CALCP (B, M, N, VP, WP, POINT) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     This subroutine computes a point on the spline surface     *      
!     at the intersection (VV,WW) of two parameter lines that    *      
!     are defined by WP and VP.                                  *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     B(3,0:3M,0:3N)  Double Precision coordinates of the        *      
!                     BEZIER-points                              *      
!     M               INTEGER  number of patches in 1st direction*      
!     N               INTEGER  number of patches in 2nd direction*      
!     VP, WP          Double precision parameter lines at whose  *      
!                     point of intersection a point of the       *      
!                     BEZIER-plain is to be determined.          *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETER:                                          *      
!     =================                                          *      
!     POINT(3)        Double Precision coordinates of the point  *      
!                     on the BEZIER surface                      *      
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
      DOUBLEPRECISION POINT (3), B (3, 0:3 * M, 0:3 * N) 
      VV = DBLE (VP * 3 * N) 
      I = INT (VV / 3) * 3 
      IF (I.GE.3 * N) I = 3 * (N - 1) 
      V = (VV - I) / 3.0D0 
      WW = DBLE (WP * 3 * M) 
      J = INT (WW / 3) * 3 
      IF (J.GE.3 * M) J = 3 * (M - 1) 
      W = (WW - J) / 3.0D0 
      F1 = (1 - V) **3 
      F2 = 3.0D0 * (1 - V) **2 * V 
      F3 = 3.0D0 * (1 - V) * V**2 
      F4 = V**3 
      F5 = (1 - W) **3 
      F6 = 3.0D0 * (1 - W) **2 * W 
      F7 = 3.0D0 * (1 - W) * W**2 
      F8 = W**3 
      DO 10 K = 1, 3 
         POINT (K) = (B (K, J, I) * F1 + B (K, J, I + 1) * F2 + B (K, J,&
         I + 2) * F3 + B (K, J, I + 3) * F4) * F5 + (B (K, J + 1, I)    &
         * F1 + B (K, J + 1, I + 1) * F2 + B (K, J + 1, I + 2) * F3 + B &
         (K, J + 1, I + 3) * F4) * F6 + (B (K, J + 2, I) * F1 + B (K, J &
         + 2, I + 1) * F2 + B (K, J + 2, I + 2) * F3 + B (K, J + 2, I + &
         3) * F4) * F7 + (B (K, J + 3, I) * F1 + B (K, J + 3, I + 1)    &
         * F2 + B (K, J + 3, I + 2) * F3 + B (K, J + 3, I + 3) * F4)    &
         * F8                                                           
   10 END DO 
      RETURN 
      END SUBROUTINE CALCP                          
