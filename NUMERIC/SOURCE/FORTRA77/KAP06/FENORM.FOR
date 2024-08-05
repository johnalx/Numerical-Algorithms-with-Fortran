C[BA*)
C[LE*)
      DOUBLE PRECISION FUNCTION FENORM (N,X)
C[IX{FENORM}*)
C
C*****************************************************************
C                                                                *
C Determining the euclidean norm of a vector X.                  *
C[BE*)
C                                                                *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C X       : N-vector X(1:N); vector for which the euclidean norm *
C           is to be determined                                  *
C N       : size of the vector X                                 *
C                                                                *
C                                                                *
C OUTPUT PARAMETER:                                              *
C =================                                              *
C FENORM  : euclidean norm of the vector X                       *
C                                                                *
C                                                                *
C LOCAL VARIABLES:                                               *
C ================                                               *
C SCPROD  : dot product of the vector X with itself              *
C I       : iteration variable                                   *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Thomas Eul                                         *
C  date     : 05.03.1985                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C   Input declarations
C
      INTEGER N
      DOUBLE PRECISION X(N)
C
C   local variables
C
      DOUBLE PRECISION SCPROD
      INTEGER I
C
      SCPROD = 0.0D0
      FENORM = 0.0D0
      DO 10 I = 1,N
         SCPROD = SCPROD + X(I)*X(I)
   10 CONTINUE
      FENORM = DSQRT(SCPROD)
C
      RETURN
      END
