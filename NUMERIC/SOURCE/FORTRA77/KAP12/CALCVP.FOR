C[BA*)
C[LE*)
      SUBROUTINE CALCVP (B,M,N,VP,ICALC,POINTS)
C[IX{CALCVP}*)
C
C*****************************************************************
C                                                                *
C     This subroutine determines points, ICALC in number, on the *
C     parameter line defined by VP.                              *
C[BE*)
C     (If VP=0.: I=0 , if VP=1.: I=3*N ; i.e., VP scales the     *
C     (MxN) - patches in the second direction N)                 *
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     B(3,0:3M,0:3N)  Double precision coordinates of the        *
C                     BEZIER-points                              *
C     M               INTEGER  number of patches in 1st direction*
C     N               INTEGER  number of patches in 2nd direction*
C     VP              Double precision parameter line on which   *
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
C  author   : Michael Radermacher                                *
C  date     : 04.30.1985                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION POINTS(3,ICALC), B(3,0:3*M,0:3*N)
      FCALC=DBLE(ICALC-1)
      DO 10   I=1,ICALC
C
C*****************************************************************
C        WP covers the interval [0,1] and defines the step size  *
C        with which points on the surface shall be computed on   *
C        parameter curve defined by VP.                          *
C*****************************************************************
C
         WP=DBLE(I-1)/FCALC
C
C*****************************************************************
C        call of SUBROUTINE CALCP for determining a point        *
C        of the surface.                                         *
C*****************************************************************
C
   10    CALL CALCP (B,M,N,VP,WP,POINTS(1,I))
      RETURN
      END
