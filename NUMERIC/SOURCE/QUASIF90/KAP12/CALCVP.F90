      SUBROUTINE CALCVP (B, M, N, VP, ICALC, POINTS) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     This subroutine determines points, ICALC in number, on the *      
!     parameter line defined by VP.                              *      
!     (If VP=0.: I=0 , if VP=1.: I=3*N ; i.e., VP scales the     *      
!     (MxN) - patches in the second direction N)                 *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     B(3,0:3M,0:3N)  Double precision coordinates of the        *      
!                     BEZIER-points                              *      
!     M               INTEGER  number of patches in 1st direction*      
!     N               INTEGER  number of patches in 2nd direction*      
!     VP              Double precision parameter line on which   *      
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
!  author   : Michael Radermacher                                *      
!  date     : 04.30.1985                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION POINTS (3, ICALC), B (3, 0:3 * M, 0:3 * N) 
      FCALC = DBLE (ICALC - 1) 
      DO 10 I = 1, ICALC 
!                                                                       
!*****************************************************************      
!        WP covers the interval [0,1] and defines the step size  *      
!        with which points on the surface shall be computed on   *      
!        parameter curve defined by VP.                          *      
!*****************************************************************      
!                                                                       
         WP = DBLE (I - 1) / FCALC 
!                                                                       
!*****************************************************************      
!        call of SUBROUTINE CALCP for determining a point        *      
!        of the surface.                                         *      
!*****************************************************************      
!                                                                       
   10 CALL CALCP (B, M, N, VP, WP, POINTS (1, I) ) 
      RETURN 
      END SUBROUTINE CALCVP                         
