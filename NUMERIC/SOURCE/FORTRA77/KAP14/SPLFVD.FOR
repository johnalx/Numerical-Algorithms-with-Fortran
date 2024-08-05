C[BA*)
C[KA{F 14}
C[  {Numerical Differentiation}
C[  {Numerical Differentiation}*)
C[FE{F 14.3}
C[  {Differentiation by Using Interpolating Polynomials}
C[  {Differentiation by Using Interpolating Polynomials}*)
C[LE*)
      SUBROUTINE SPLFVD (X,N,XN,A,B,C,D,S,S1,S2,S3)
C[IX{SPLFVD}*)
C
C*****************************************************************
C                                                                *
C  Subroutine SPLFVD determines the functional value as well as  *
C  the 1st, 2nd and 3rd derivatives of a cubic spline function S *
C[BE*)
C  given in the following form                                   *
C                                                                *
C  S(X) = P(I)(X) = A(I) + B(I)*(X-XN(I)) + C(I)*(X-XN(I))**2 +  *
C                        + D(I)*(X-XN(I))**3,                    *
C                                                                *
C  where X lies in the interval [XN(I), XN(I+1)], I=0, ..., N-1. *
C                                                                *
C  For X < XN(0) the boundary polynomial P(0) is evaluated,      *
C  for X > XN(N) we use the boundary polynomial P(N-1).          *
C  No plausibility check of the input values is performed.       *
C  This program may be used in connection with other spline      *
C  routines.                                                     *
C                                                                *
C                                                                *
C  INPUT PARAMETER:                                              *
C  ================                                              *
C  X  : value where S(X), S'(X), S''(X), S'''(X) are to be       *
C       determined.                                              *
C  N  : index of the last node.                                  *
C  XN : (N+1)-vector XN(0:N) containing the nodes XN(I),         *
C       I=0, ..., N.                                             *
C  A  : ] (N+1)-vectors A(0:N), ..., D(0:N), containing the      *
C  B  : ] spline coefficients A(I), B(I), C(I), D(I) for         *
C  C  : ] I=0, ..., N-1 and an auxiliary variable in position    *
C  D  : ] I=N.                                                   *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  S  : functional value S(X) of the spline function S at X.     *
C  S1 : 1st derivative S'(X) of the spline function S at X.      *
C  S2 : 2nd derivative S''(X) of the spline function S at X.     *
C  S3 : 3rd derivative S'''(X) of the spline function S at X.    *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Gisela Engeln-Muellges                             *
C  date     : 04.12.1988                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  declarations
C
      DOUBLE PRECISION XN(0:N),A(0:N),B(0:N),C(0:N),D(0:N)
      SAVE I
C
C  initializing
C
      DATA I /0/
      IF (I .GE. N) I=0
C
C  If this is a repeated call of the subroutine, the
C  loop for determining the interval [XN(I), XN(I+1)],
C  that contains X, is only executed if X is not
C  inside the same interval as during the previous call.
C
      IF (X .LT. XN(I) .OR. X .GE. XN(I+1)) THEN
         I = 0
         K = N
   10    M = (I+K)/2
         IF (X .LT. XN(M)) THEN
            K = M
         ELSE
            I = M
         ENDIF
         IF (K .GT. I+1) GOTO 10
      ENDIF
C
C  determine the linear factor X-XN(I) needed for the polynomial evaluation
C
      XL = X - XN(I)
C
C  determine the functional and derivative values via a HORNER
C  scheme using certain auxiliary variables
C
      DUMMY1 = 3.0D0*D(I)
      DUMMY2 = 2.0D0*C(I)
      DUMMY3 = 2.0D0*DUMMY1
      S  = ((D(I)*XL + C(I))*XL + B(I))*XL + A(I)
      S1 = (DUMMY1*XL + DUMMY2)*XL + B(I)
      S2 = DUMMY3*XL + DUMMY2
      S3 = DUMMY3
      RETURN
      END
