C[BA*)
C[LE*)
      DOUBLE PRECISION FUNCTION VMNORM (VECTOR,N)
C[IX{VMNORM}*)
C
C*****************************************************************
C                                                                *
C This FUNCTION routine determines the maximum norm of a vector. *
C[BE*)
C                                                                *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C VECTOR  : DOUBLE PRECISION vector VECTOR(1:N)                  *
C N       : size of the vector                                   *
C                                                                *
C                                                                *
C LOCAL VARIABLE:                                                *
C ===============                                                *
C I       : loop variable                                        *
C                                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  Author  : Volker KrÅger                                     *
C  Date    : 07.08.1990                                        *
C  Source  : FORTRAN 77                                        *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C Declaration
C
      DOUBLE PRECISION VECTOR(N)
C
C Compute maximum norm
C
      VMNORM=DABS(VECTOR(1))
      DO 10 I=2,N
         VMNORM=DMAX1(VMNORM,DABS(VECTOR(I)))
10    CONTINUE
      RETURN
      END
