C[BA*)
C[LE*)
      DOUBLE PRECISION FUNCTION SPVAL (X,N,XN,A,B,C,D)
C[IX{SPVAL}*)
C
C*****************************************************************
C                                                                *
C  The FUNCTION SPVAL determines the functional value of a cubic *
C  spline S                                                      *
C[BE*)
C  of the form:                                                  *
C                                                                *
C  S(X) = P(I)(X) = A(I) + B(I)*(X-XN(I)) + C(I)*(X-XN(I))**2 +  *
C                                         + D(I)*(X-XN(I))**3    *
C                                                                *
C  for X in the interval [XN(I),XN(I+1)], where I=0,1,...,N-1.   *
C                                                                *
C  If X < XN(0) we evaluate the first polynom P(0), if X > XN(N) *
C  we evaluate the polynom P(N-1).                               *
C  We do not check the input X in any way.                       *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  X  :  value for which we want to find the value of the spline *
C  N  :  index of the final node                                 *
C  XN :  vector of nodes XN(0:N)                                 *
C  A  :  ] N+1-vectors ..(0:N);                                  *
C  B  :  ] the spline coefficients A(I), B(I), C(I), D(I) for    *
C  C  :  ] I = 0, ..., N-1; used for storage for I = N.          *
C  D  :  ]                                                       *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : GÅnter Palm                                        *
C  Date     : 01.06.1991                                         *
C  Source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  Declarations
C
      DOUBLE PRECISION XN(0:N), A(0:N), B(0:N), C(0:N), D(0:N)
      SAVE I
C
C  Initializing
C
      DATA I /0/
      IF (I .GE. N) I=0
C
C  In a repeat call we determine the interval [XN(I), XN(I+1)]
C  that contains X only if X does not lie in the same interval
C  as last used.
C
      IF (X .LT. XN(I)  .OR.  X .GE. XN(I+1)) THEN
        I = 0
        K = N
   10   M = (I+K)/2
        IF (X .LT. XN(M)) THEN
          K = M
        ELSE
          I = M
        ENDIF
        IF (K .GT. I+1) GOTO 10
      ENDIF
C
C  Compute X-XN(I) in order to evaluate the polynomial
C
      XL = X - XN(I)
C
C  Compute the value of the spline via a Horner scheme
C
      SPVAL = ((D(I)*XL + C(I))*XL + B(I))*XL + A(I)
      RETURN
      END
