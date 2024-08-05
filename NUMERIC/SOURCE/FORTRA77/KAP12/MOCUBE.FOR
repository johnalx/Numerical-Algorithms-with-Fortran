C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 12.3.3}
C[  {Modified (Interpolating) Cubic B\'ezier Splines}
C[  {Modified (Interpolating) Cubic B\'ezier Splines}*)
C[LE*)
      SUBROUTINE MOCUBE (D,B,M,EPS)
C[IX{MOCUBE}*)
C
C*****************************************************************
C                                                                *
C     SUBROUTINE MOCUBE determines the coefficients of a         *
C     modified BEZIER spline from the weight points so that the  *
C     weight points will be located on the computed curve, except*
C     for an accuracy of EPS.                                    *
C[BE*)
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     M                : number of curve segments of the spline  *
C     EPS              : accuracy bound for the interpolation    *
C     D(1,J), J=0,...,M: X-coordinates of the weight points      *
C     D(2,J), M=0,...,M: Y-coordinates of the weight points      *
C                                                                *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     B(1,J), J=0,1,...,3*M: X-coordinates of the BEZIER points  *
C     B(2,J), J=0.1,...,3*M: Y-coordinates of the BEZIER points  *
C                                                                *
C                                                                *
C     LOCAL VARIABLES:                                           *
C     ================                                           *
C     DINT:  difference of a weight points and the corresponding *
C            BEZIER point                                        *
C     J,K:   control variables                                   *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Juergen Dietel                                     *
C  date     : 04.23.1987                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER M,J,K
      DOUBLE PRECISION D(2,0:M),B(2,0:3*M),EPS,DINT
C
C*****************************************************************
C     determine the BEZIER points                                *
C     B(1,J), J=0,1,...,3*M (X-direction) and                    *
C     B(2,J), J=0,1,...,3*M (Y-direction)                        *
C*****************************************************************
C
      DO 40 K=1,2
         DO 30 J=1,M-1
            B(K,3*J-2) = (2.0D0*D(K,J-1)+D(K,J))/3.0D0
            B(K,3*J  ) = (D(K,J-1)+4.0D0*D(K,J)+D(K,J+1))/6.0D0
            B(K,3*J+2) = (D(K,J)+2.0D0*D(K,J+1))/3.0D0
   30       CONTINUE
         B(K,0)     =        D(K,0)
         B(K,2)     = (      D(K,0)+2.0D0*D(K,1))/3.0D0
         B(K,3*M-2) = (2.0D0*D(K,M-1)+    D(K,M))/3.0D0
         B(K,3*M)   =        D(K,M)
   40    CONTINUE
C
C*****************************************************************
C     correction of the BEZIER points                            *
C*****************************************************************
C
   50 DO 90 K=1,2
         DO 80 J=3,3*M-3,3
            DINT=D(K,J/3)-B(K,J)
            IF (J .NE. 3)     B(K,J-3) = B(K,J-3)+DINT/4.0D0
                              B(K,J-2) = B(K,J-2)+DINT/2.0D0
                              B(K,J-1) = B(K,J-1)+DINT
                              B(K,J  ) = B(K,J  )+DINT
                              B(K,J+1) = B(K,J+1)+DINT
                              B(K,J+2) = B(K,J+2)+DINT/2.0D0
            IF (J .NE. 3*M-3) B(K,J+3) = B(K,J+3)+DINT/4.0D0
   80       CONTINUE
   90    CONTINUE
C
C*****************************************************************
C     check whether interpolation has reached an accuracy of EPS,*
C     otherwise repeat correcting the BEZIER points              *
C*****************************************************************
C
      DO 100 J=1,M-1
         IF (DABS(D(1,J)-B(1,3*J))+DABS(D(2,J)-B(2,3*J)) .GT. EPS)
     *      GOTO 50
  100    CONTINUE
      RETURN
      END
