      SUBROUTINE CALCWP (B, M, N, WP, ICALC, POINTS) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     This subroutine determines points, ICALC in number, on the *      
!     parameter line defined by WP.                              *      
!     (If WP=0.: I=0 , if WP=1.: I=3*N ; i.e., WP scales the     *      
!     (MxN) - patches in the first direction M)                  *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     B(3,0:3M,0:3N)  Double precision coordinates of the        *      
!                     BEZIER-points                              *      
!     M               INTEGER  number of patches in 1st direction*      
!     N               INTEGER  number of patches in 2nd direction*      
!     WP              Double precision parameter line on which   *      
!                     intermediate points of the BEZIER surface  *      
!                     are to be determined                       *      
!     ICALC           INTEGER  number of points to be determined *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETER:                                          *      
!     =================                                          *      
!     POINTS(3,ICALC) Double precision coordinates of the com-   *      
!                     puted intermediate points                  *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: CALCP                                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : Michael Radermacher                                *      
!  Date     : 04.30.1985                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION POINTS (3, ICALC), B (3, 0:3 * M, 0:3 * N) 
      FCALC = DBLE (ICALC - 1) 
      DO 10 I = 1, ICALC 
!                                                                       
!*****************************************************************      
!        VP covers the interval [0,1] and defines the step size  *      
!        with which points on the surface shall be computed on   *      
!        the parameter curve defined by WP.                      *      
!*****************************************************************      
!                                                                       
         VP = DBLE (I - 1) / FCALC 
!                                                                       
!*****************************************************************      
!        call of the SUBROUTINE CALCP for determining a point    *      
!        of the surface.                                         *      
!*****************************************************************      
!                                                                       
   10 CALL CALCP (B, M, N, VP, WP, POINTS (1, I) ) 
      RETURN 
      END SUBROUTINE CALCWP                         
