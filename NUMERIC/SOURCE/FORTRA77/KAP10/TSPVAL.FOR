C[BA*)
C[LE*)
      DOUBLE PRECISION FUNCTION TSPVAL (PHI,N,PHIN,A,B,C,D)
C[IX{TSPVAL}*)
C
C*****************************************************************
C                                                                *
C  The FUNCTION TSPVAL evaluates the value of a transformed      *
C  parametric cubic spline                                       *
C[BE*)
C  given in the form:                                            *
C                                                                *
C  S(PHI) = A(I) + B(I)(PHI-PHIN(I)) + C(I)(PHI-PHIN(I))**2 +    *
C                                    + D(I)(PHI-PHIN(I))**3      *
C                                                                *
C  for PHI in the interval [PHIN(I),PHIN(I+1)], I=0,1,...,N-1.   *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  PHI  :  value where the spline shall be evaluated (in radians)*
C  N    :  index of the lst node                                 *
C  PHIN :  vector PHIN(0:N); the nodes PHIN(I), I=0,1,...,N      *
C  A    :  ] N+1-vectors ..(0:N);                                *
C  B    :  ] the elements in positions 0 to N-1 denote the       *
C  C    :  ] coefficients of the spline.                         *
C  D    :  ]                                                     *
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
      DOUBLE PRECISION PHIN(0:N), A(0:N), B(0:N), C(0:N), D(0:N)
C
C-----Initialize
C
      ZWOPI = 8.0D0*DATAN(1.0D0)
C
C-----Assign PHI to PHIX,
C     if necessary, recompute PHIX so that it lies
C     in the interval [0,2*PI]
C
      IF (PHI .LT. 0.0D0) THEN
        L = ABS(INT(PHI/ZWOPI)) + 1
        PHIX = L*ZWOPI - PHI
      ELSEIF (PHI .GT. ZWOPI) THEN
        L = INT(PHI/ZWOPI)
        PHIX = PHI - L*ZWOPI
      ELSE
        PHIX = PHI
      ENDIF
C
C-----Compute the functional value and that of the derivatives
C     at PHIX in SUBROUTINE SPVAL
C
      TSPVAL = SPVAL (PHIX,N,PHIN,A,B,C,D)
C
      RETURN
      END
