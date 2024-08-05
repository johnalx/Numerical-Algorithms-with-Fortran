C[BA*)
C[LE*)
      SUBROUTINE ISPL1D (N,XN,FN,ALPHA,BETA,MREP,B,C,D,
C[IX{ISPL1D}*)
     +                   H,DU,DM,RS,IERR)
C
C*******************************************************************
C                                                                  *
C  ISPL1D computes the coefficients B(I), C(I), D(I) for I=0,1,.., *
C  N-1 of a cubic interpolating spline with prescribed first end   *
C  point derivatives.                                              *
C[BE*)
C                                                                  *
C  The spline has the form:                                        *
C                                                                  *
C  S(X) = FN(I) + B(I)(X-XN(I)) + C(I)(X-XN(I))**2 +               *
C                               + D(I)(X-XN(I))**3                 *
C                                                                  *
C  for X in the interval [XN(I),XN(I+1)] for I=0,1,...,N-1.        *
C                                                                  *
C                                                                  *
C  ASSUMPTIONS:   1.         N > 2                                 *
C  ============   2.     XN(I) < XN(I+1), I=0,1,...,N-1            *
C                                                                  *
C                                                                  *
C  NOTE:  ISPL1D should not be used by itself, but rather via the  *
C  =====  SUBROUTINE ISPLNP, or ISPLPA for parametric splines.     *
C         These subroutines also check the input data.             *
C                                                                  *
C                                                                  *
C  INPUT PARAMETERS:                                               *
C  =================                                               *
C  N  :  index of the final node                                   *
C  XN :  vector XN(0:N); the nodes XN(I), I = 0,1,...,N            *
C  FN :  vector FN(0:N); the functional values FN(I) = FN( XN(I) ) *
C                                                                  *
C  ALPHA :  first derivative at XN(0)                              *
C  BETA  :  first derivative at XN(N)                              *
C                                                                  *
C  MREP  : index for repeated call of the SUBROUTINE:              *
C          MREP = 1: In order to compute C(I) the system matrix and*
C                    its factorization must be evaluated in SUB-   *
C                    ROUTINE TRDSY anew.                           *
C          MREP = 2: Only the right hand side of the system needs  *
C                    to be recomputed. The vectors DU and DM as    *
C                    computed during the first call in SUBROUTINE  *
C                    TRDSYS help us find the solution.             *
C                    (This prevents a duplicate decomposition in   *
C                    case of parametric splines).                  *
C                    The entries of H, DU and DM must not be       *
C                    altered after the first call.                 *
C                                                                  *
C                                                                  *
C  AUXILIARY VARIABLES:                                            *
C  ====================                                            *
C  H  :   N-vector H(0:N-1)                                        *
C  DU : ]                                                          *
C  DM : ] N-1-vectors ..(1:N-1)                                    *
C  RS : ]                                                          *
C                                                                  *
C                                                                  *
C  OUTPUT PARAMETERS:                                              *
C  ==================                                              *
C  FN :  ]  N+1-vectors ..(0:N);                                   *
C  B  :  ]  The first N entries of B, C and D are the spline       *
C  C  :  ]  coefficients for S. B(N), C(N), D(N) are auxiliary     *
C  D  :  ]  variables.                                             *
C  IERR :  error parameter                                         *
C          =  0 :  All is ok                                       *
C          = -1 :  N < 3                                           *
C          = -5 :  erroneous value for MREP                        *
C          =  1 :  crash in SUBROUTINE TRDSY, system matrix        *
C                  numerically singular                            *
C                                                                  *
C------------------------------------------------------------------*
C                                                                  *
C  Subroutines required: TRDSY, TRDSYS                             *
C                                                                  *
C                                                                  *
C                                                                  *
C*******************************************************************
C                                                                  *
C  Author   : Gnter Palm                                          *
C  Date     : 04.15.1988                                           *
C  Source   : FORTRAN 77                                           *
C                                                                  *
C[BA*)
C*******************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION XN(0:N), FN(0:N), B(0:N), C(0:N), D(0:N),
     +                 H(0:N-1), DU(1:N-1), DM(1:N-1), RS(1:N-1)
C
C-----Check MREP for repeated call
C
      IERR = -5
      IF (MREP .NE. 1  .AND.  MREP .NE. 2) RETURN
C
C-----Compute auxiliary variables and matrix elements
C     on and above the diagonal in case of first call
C
      IF (MREP .EQ. 1) THEN
C
C       auxiliary variables
C
        DO 10 I=0,N-1,1
          H(I) = XN(I+1) - XN(I)
   10   CONTINUE
C
C       co diagonal
C
        DO 20 I=1,N-2,1
          DU(I) = H(I)
   20   CONTINUE
C
C       main diagonal
C
        DM(1) = 1.5D0*H(0) + 2.D0*H(1)
        DO 30 I=2,N-2,1
          DM(I) = 2.0D0*(H(I-1)+H(I))
   30   CONTINUE
        DM(N-1) = 2.0D0*H(N-2) + 1.5D0*H(N-1)
      ENDIF
C
C-----Compute right hand side
C
      DUMMY1 = (FN(2) - FN(1))/H(1)
      RS(1)  = 3.0D0*(DUMMY1 - 0.5D0*(3.0D0*(FN(1)-FN(0))/H(0)
     +         - ALPHA))
      DO 40 I = 2,N-2,1
        DUMMY2 = (FN(I+1) - FN(I))/H(I)
        RS(I)  = 3.0D0*(DUMMY2 - DUMMY1)
        DUMMY1 = DUMMY2
   40 CONTINUE
      RS(N-1) = 3.0D0*(0.5D0*(3.0D0*(FN(N)-FN(N-1))/H(N-1)-BETA)
     +          - DUMMY1)
C
C-----Solve the linear system to obtain C(1), ... ,C(N-1)
C
      IF (MREP .EQ. 1) THEN
C
C       ... factor the system matrix
C
        CALL TRDSY (N-1,DM,DU,RS,C(1),IFLAG)
        IF (IFLAG .NE. 1) THEN
          IF (IFLAG .EQ. -2) THEN
            IERR = -1
          ELSE
            IERR =  1
          ENDIF
          RETURN
        ENDIF
      ELSE
C
C       ... without factoring (repeated call)
C
        CALL TRDSYS (N-1,DM,DU,RS,C(1))
      ENDIF
      IERR = 0
C
C-----Compute the remaining spline coefficients
C
      C(0) = 1.5D0*((FN(1)-FN(0))/H(0) - ALPHA)/H(0) - 0.5D0*C(1)
      C(N) = 1.5D0*(BETA-(FN(N)-FN(N-1))/H(N-1))/H(N-1) - 0.5D0*C(N-1)
C
      DO 50 I=0,N-1,1
        B(I) = (FN(I+1)-FN(I))/H(I) - H(I)/3.0D0*(C(I+1)+2.0D0*C(I))
        D(I) = (C(I+1)-C(I))/(3.0D0*H(I))
   50 CONTINUE
      RETURN
      END
