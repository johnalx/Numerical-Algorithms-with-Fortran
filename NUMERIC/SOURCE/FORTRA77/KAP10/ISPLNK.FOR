C[BA*)
C[LE*)
      SUBROUTINE ISPLNK (N,XN,FN,B,C,D,H,DM,RS,IERR)
C[IX{ISPLNK}*)
C
C*******************************************************************
C                                                                  *
C  ISPLNK computes the coefficients B(I), C(I), D(I) for I=0,1,.., *
C  N-1 of a cubic interpolating spline with 'not a node' end       *
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
C  NOTE:  ISPLNK should not be used by itself, but rather via the  *
C  =====  SUBROUTINE ISPLNP, which also checks the input data.     *
C                                                                  *
C                                                                  *
C  INPUT PARAMETERS:                                               *
C  =================                                               *
C  N  :  index of the final node                                   *
C  XN :  vector XN(0:N); the nodes XN(I), I = 0,1,...,N            *
C  FN :  vector FN(0:N); the functional values FN(I) = FN( XN(I) ) *
C                                                                  *
C                                                                  *
C  AUXILIARY VARIABLES:                                            *
C  ====================                                            *
C  H  :   N-vector H(0:N-1)                                        *
C  DM :   N-1-vector DM(1:N-1)                                     *
C  RS :   N-1-vector RS(1:N-1)                                     *
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
C          =  1 :  crash in SUBROUTINE TRDIG, system matrix        *
C                  numerically singular                            *
C                                                                  *
C------------------------------------------------------------------*
C                                                                  *
C  Subroutines required: TRDIG                                     *
C                                                                  *
C                                                                  *
C                                                                  *
C*******************************************************************
C                                                                  *
C  Author   : G�nter Palm                                          *
C  Date     : 04.15.1988                                           *
C  Source   : FORTRAN 77                                           *
C                                                                  *
C[BA*)
C*******************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION XN(0:N), FN(0:N), B(0:N), C(0:N), D(0:N),
     +                 H(0:N-1), DM(1:N-1), RS(1:N-1)
C
C-----Computing auxiliary variables
C
      DO 10 I=0,N-1,1
        H(I) = XN(I+1) - XN(I)
   10 CONTINUE
C
C-----Compute the main and both co diagonals of the system
C     matrix and the right hand side of A*C=RS with a
C     tridiagonal matrix A
C
C     co diagonals
C
      D(1) = H(1) - H(0)
      B(2) = H(0)
      DO 20 I=2,N-3,1
        D(I)   = H(I)
        B(I+1) = H(I)
   20 CONTINUE
      D(N-2) = H(N-2)
      B(N-1) = H(N-2) - H(N-1)
C
C     main diagonal
C
      DM(1) = H(0) + 2.0D0*H(1)
      DO 30 I=2,N-2,1
        DM(I) = 2.0D0*(H(I-1) + H(I))
   30 CONTINUE
      DM(N-1) = 2.0D0*H(N-2) + H(N-1)
C
C     right hand side
C
      DUMMY1 = (FN(2) - FN(1))/H(1)
      RS(1)  = 3.0D0*H(1)/(H(1)+H(0))*(DUMMY1 - (FN(1)-FN(0))/H(0))
      DO 40 I = 2,N-2,1
        DUMMY2 = (FN(I+1) - FN(I))/H(I)
        RS(I)  = 3.0D0*(DUMMY2 - DUMMY1)
        DUMMY1 = DUMMY2
   40 CONTINUE
      RS(N-1) = 3.0D0*H(N-2)/(H(N-2)+H(N-1))*((FN(N)-FN(N-1))/H(N-1)
     +          - DUMMY1)
C
C-----Solve the linear system to find C(1), ..., C(N-1)
C
      CALL TRDIG (N-1,B(1),DM,D(1),RS,C(1),IFLAG)
      IF (IFLAG .NE. 1) THEN
        IF (IFLAG .EQ. 0) THEN
          IERR = 1
        ELSE
          IERR = -1
        ENDIF
      ENDIF
      IERR = 0
C
C-----Compute the remaining spline coefficients
C
      C(0) = C(1) + H(0)/H(1)*(C(1) - C(2))
      C(N) = C(N-1) + H(N-1)/H(N-2)*(C(N-1) - C(N-2))
C
      DO 50 I=0,N-1,1
        B(I) = (FN(I+1)-FN(I))/H(I) - H(I)/3.0D0*(C(I+1)+2.0D0*C(I))
        D(I) = (C(I+1)-C(I))/(3.0D0*H(I))
   50 CONTINUE
      RETURN
      END
