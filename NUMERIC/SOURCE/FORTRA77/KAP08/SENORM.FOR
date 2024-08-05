C[BA*)
C[LE*)
      DOUBLE PRECISION FUNCTION SENORM (X, N)
C[IX{SENORM}*)
C
C*****************************************************************
C                                                                *
C  The FUNCTION-subroutine SENORM computes the square of the     *
C  euclidean norm of a real vector X of length N+1               *
C[BE*)
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C                                                                *
C  X    (N+1)-vector X(0:N) containing the vector whose norm     *
C       shall be computed                                        *
C  N    number of components of the vector X                     *
C                                                                *
C  OUTPUT PARAMETER:                                             *
C  =================                                             *
C                                                                *
C  none                                                          *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Ilona Westermann                                   *
C  date     : 09.01.1987                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N
      DOUBLE PRECISION X(0:N)
      SENORM = 0.0D0
      DO 10 I=0,N
         SENORM = SENORM + X(I) * X(I)
   10 CONTINUE
      RETURN
      END
