C[BA*)
C[LE*)
      SUBROUTINE PSPVAL (TV,N,T,AX,BX,CX,DX,AY,BY,CY,DY,SX,SY)
C[IX{PSPVAL}*)
C
C*****************************************************************
C                                                                *
C  This SUBROUTINE evaluates a parametric cubic spline           *
C[BE*)
C  given by its component functions SX(T), SY(T) in the form:    *
C                                                                *
C  SX := SX(T) = AX(I) + BX(I)(T-T(I)) + CX(I)(T-T(I))**2 +      *
C                                      + DX(I)(T-T(I))**3        *
C                                                                *
C  SY := SY(T) = AY(I) + BY(I)(T-T(I)) + CY(I)(T-T(I))**2 +      *
C                                      + DY(I)(T-T(I))**3        *
C                                                                *
C  for T in the interval [T(I),T(I+1)], I=0,1,...,N-1.           *
C                                                                *
C  The program computes the values of SX(T) and SY(T) at T=TV.   *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  TV :  value where the component functions SX(T) and SY(T)     *
C        are to be evaluated                                     *
C  N  :  index of the last node                                  *
C  T  :  vector T(0:N); the nodes T(I), I=0,1,...,N              *
C  AX :  ] N+1-vectors ..(0:N);                                  *
C  BX :  ] their elements 0, ..., N-1 describe the component     *
C  CX :  ] function SX(T).                                       *
C  DX :  ]                                                       *
C                                                                *
C  AY :  ] N+1-vectors ..(0:N);                                  *
C  BY :  ] their elements 0, ..., N-1 describe the component     *
C  CY :  ] function SY(T).                                       *
C  DY :  ]                                                       *
C          The elements BX(N), CX(N), DX(N), BY(N), CY(N) and    *
C          DY(N) are used for auxiliary purposes.                *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  SX :  Functional value for SX(T) at T=TV                      *
C  SY :  Functional value for SY(T) at T=TV                      *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Subroutines required: SPVAL                                   *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : GÅnter Palm                                        *
C  Date     : 06.01.1991                                         *
C  Source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C-----Declarations
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION T(0:N), AX(0:N), BX(0:N), CX(0:N), DX(0:N),
     +                         AY(0:N), BY(0:N), CY(0:N), DY(0:N)
C
C-----Compute the functional values of the component splines
C     SX(T) and SY(T)
C
      SX = SPVAL (TV,N,T,AX,BX,CX,DX)
      SY = SPVAL (TV,N,T,AY,BY,CY,DY)
C
      RETURN
      END
