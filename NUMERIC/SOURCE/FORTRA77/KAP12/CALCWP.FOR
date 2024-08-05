C[BA*)
C[LE*)
      SUBROUTINE CALCWP (B,M,N,WP,ICALC,POINTS)
C[IX{CALCWP}*)
C
C*****************************************************************
C                                                                *
C     This subroutine determines points, ICALC in number, on the *
C     parameter line defined by WP.                              *
C[BE*)
C     (If WP=0.: I=0 , if WP=1.: I=3*N ; i.e., WP scales the     *
C     (MxN) - patches in the first direction M)                  *
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     B(3,0:3M,0:3N)  Double precision coordinates of the        *
C                     BEZIER-points                              *
C     M               INTEGER  number of patches in 1st direction*
C     N               INTEGER  number of patches in 2nd direction*
C     WP              Double precision parameter line on which   *
C                     intermediate points of the BEZIER surface  *
C                     are to be determined                       *
C     ICALC           INTEGER  number of points to be determined *
C                                                                *
C                                                                *
C     OUTPUT PARAMETER:                                          *
C     =================                                          *
C     POINTS(3,ICALC) Double precision coordinates of the com-   *
C                     puted intermediate points                  *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: CALCP                                   *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : Michael Radermacher                                *
C  Date     : 04.30.1985                                         *
C  Source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION POINTS(3,ICALC) , B(3,0:3*M,0:3*N)
      FCALC=DBLE(ICALC-1)
      DO 10   I=1,ICALC
C
C*****************************************************************
C        VP covers the interval [0,1] and defines the step size  *
C        with which points on the surface shall be computed on   *
C        the parameter curve defined by WP.                      *
C*****************************************************************
C
         VP=DBLE(I-1)/FCALC
C
C*****************************************************************
C        call of the SUBROUTINE CALCP for determining a point    *
C        of the surface.                                         *
C*****************************************************************
C
   10    CALL CALCP (B,M,N,VP,WP,POINTS(1,I))
      RETURN
      END
