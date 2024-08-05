C[BA*)
C[LE*)
      SUBROUTINE CALCP (B,M,N,VP,WP,POINT)
C[IX{CALCP}*)
C
C*****************************************************************
C                                                                *
C     This subroutine computes a point on the spline surface     *
C     at the intersection (VV,WW) of two parameter lines that    *
C     are defined by WP and VP.                                  *
C[BE*)
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     B(3,0:3M,0:3N)  Double Precision coordinates of the        *
C                     BEZIER-points                              *
C     M               INTEGER  number of patches in 1st direction*
C     N               INTEGER  number of patches in 2nd direction*
C     VP, WP          Double precision parameter lines at whose  *
C                     point of intersection a point of the       *
C                     BEZIER-plain is to be determined.          *
C                                                                *
C                                                                *
C     OUTPUT PARAMETER:                                          *
C     =================                                          *
C     POINT(3)        Double Precision coordinates of the point  *
C                     on the BEZIER surface                      *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
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
      DOUBLE PRECISION POINT(3), B(3,0:3*M,0:3*N)
      VV=DBLE(VP*3*N)
      I=INT(VV/3)*3
      IF ( I .GE. 3*N) I=3*(N -1)
      V=(VV-I)/3.0D0
      WW=DBLE(WP*3*M)
      J=INT(WW/3)*3
      IF ( J .GE. 3*M) J=3*(M -1)
      W=(WW-J)/3.0D0
         F1=(1-V)**3
         F2=3.0D0*(1-V)**2*V
         F3=3.0D0*(1-V)*V**2
         F4=V**3
         F5=(1-W)**3
         F6=3.0D0*(1-W)**2*W
         F7=3.0D0*(1-W)*W**2
         F8=W**3
         DO 10   K=1,3
            POINT(K)=(B(K,J,I)*F1+B(K,J,I+1)*F2+
     F                B(K,J,I+2)*F3+B(K,J,I+3)*F4)*F5
     F              +(B(K,J+1,I)*F1+B(K,J+1,I+1)*F2+
     F                B(K,J+1,I+2)*F3+B(K,J+1,I+3)*F4)*F6
     F              +(B(K,J+2,I)*F1+B(K,J+2,I+1)*F2+
     F                B(K,J+2,I+2)*F3+B(K,J+2,I+3)*F4)*F7
     F              +(B(K,J+3,I)*F1+B(K,J+3,I+1)*F2+
     F                B(K,J+3,I+2)*F3+B(K,J+3,I+3)*F4)*F8
   10       CONTINUE
      RETURN
      END
