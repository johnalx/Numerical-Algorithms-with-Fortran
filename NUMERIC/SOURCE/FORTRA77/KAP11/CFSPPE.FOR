C[BA*)
C[LE*)
      SUBROUTINE CFSPPE (N,XN,FN,W,MREP,A,B,C,D,H,H1,H2,H3,RS,
C[IX{CFSPPE}*)
     +                   WORK,IERR)
C
C*****************************************************************
C                                                                *
C  'CFSPPE' computes the coefficients A(I), B(I), C(I), D(I) for *
C  I = 0, 1, ..., N-1 of a periodic cubic least squares spline.  *
C[BE*)
C  The spline function is represented in the form:               *
C                                                                *
C  S(X) = A(I) + B(I)(X-XN(I)) + C(I)(X-XN(I))**2 +              *
C                              + D(I)(X-XN(I))**3                *
C                                                                *
C  for X in the interval [XN(I),XN(I+1)] for I=0,1,...,N-1.      *
C                                                                *
C                                                                *
C  ASSUMPTIONS:        1.         N > 5                          *
C  ================    2.     XN(I) < XN(I+1), I=0,1,...,N-1     *
C                      3.      W(I) > 0.0    , I=0,1,...,N       *
C                      4.      W(0) = W(N)                       *
C                      5.     FN(0) = FN(N)                      *
C                                                                *
C                                                                *
C  BEMERKUNG:  'CFSPPE' should not be called by itself, but      *
C  ==========  rather using the subroutine 'CFSPNP' - or in      *
C              case of parametric or transformed parametric      *
C              splines via the subroutines 'CFSPPA' or 'CFSPTR'. *
C              These subroutines check the assumptions 1 to 3.   *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  N  :  Index of final node                                     *
C  XN :  (N+1)-vector XN(0:N); XN(I) = nodes, I = 0,1,...,N      *
C  FN :  (N+1)-vector FN(0:N); FN(I) = function values at XN(I)  *
C  W  :  (N+1)-vector W(0:N);   W(I) = weights for FN(I)         *
C                                                                *
C  MREP :  Marker for repeated call of this subroutine:          *
C          MREP = 1: We must form the system matrix and its      *
C                    complete LR decomposition in order to       *
C                    compute C(I) in the subroutine NCYFSY.      *
C          MREP = 2: We only have to compute the right hand side.*
C                    We can use the vector WORK from our first   *
C                    call of NCYFSS to find the solution.        *
C                    (This avoids a duplicate LR decomposition   *
C                    for parametric splines). If MREP = 2, the   *
C                    entries of the vectors H, H1, H2, H3 and    *
C                    WORK must not be altered after the first    *
C                    call of NCYFSS !                            *
C  REMARK:    For parametric splines with differing weights      *
C             (WX(I) different from WY(I) for at least one index *
C             I) put MREP = 1 always!                            *
C                                                                *
C                                                                *
C  AUX VARIABLES:                                                *
C  ==============                                                *
C  H    :]                                                       *
C  H1   :]  (N+1)-vectors ..(0:N)                                *
C  H2   :]                                                       *
C  H3   :]                                                       *
C                                                                *
C  RS   :  N-vector RS(1:N)                                      *
C  WORK :  vector WORK(1:9*N-14)                                 *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  A    :  ] (N+1)-vectors    The elements in positions 0 to N-1 *
C  B    :  ] ..(0:N)          are the coefficients of the spline *
C  C    :  ]                  function S.                        *
C  D    :  ]                  A(N), B(N), C(N), D(N) are aux     *
C                             variables.                         *
C  IERR :  error parameter                                       *
C          =  0 :  All ok                                        *
C          = -1 :  N < 6                                         *
C          = -5 :  FN(0) not equal to FN(N),  or                 *
C                   W(0) not equal to W(N)                       *
C          = -6 :  wrong value assigned to MREP                  *
C          =  1 :  fatal error in NCYFSY                         *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Subroutines required : NCYFSY, NCYFSS                         *
C                                                                *
C                                                                *
C  Reference : Engeln-Mllges, G.; Reutter, F., see [ENGE87].    *
C                                                                *
C*****************************************************************
C                                                                *
C  Author      : Gnter Palm                                     *
C  Date        : 18.04.1988                                      *
C  Source code : FORTRAN 77                                      *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XN(0:N), FN(0:N), W(0:N), A(0:N), B(0:N),
     +                 C(0:N), D(0:N), H(0:N), H1(0:N), H2(0:N),
     +                 H3(0:N), RS(1:N), WORK(1:9*N-14)
C
C-----Check input for periodicity---------------
C
      IERR = -5
      IF (FN(N) .NE. FN(0)) RETURN
      IF (W(N)  .NE. W(0))  RETURN
C
C-----Check marker for a repeated call---------------
C
      IERR = -6
      IF (MREP .NE. 1  .AND.  MREP .NE. 2) RETURN
C
C-----Compute aux variables and matrix elements-----------
C     Main diagonal, lower, and upper co-diagonals of system matrix
C     for linear system  (this matrix is symmetric, almost cyclic
C     and five-diagonal),  in case of a first call.
C
      IF (MREP .EQ. 1) THEN
C
C       Aux variables
C
        DO 10 I=0,N-1,1
          H(I)  = XN(I+1) - XN(I)
          H1(I) = 1.0D0/H(I)
          C(I)  = H1(I)*H1(I)
          H2(I) = 6.0D0/W(I)
   10   CONTINUE
        H(N)  = H(0)
        H1(N) = H1(0)
        C(N)  = C(0)
        H2(N) = H2(0)
C
        DO 20 I=0,N-1,1
          H3(I) = H1(I) + H1(I+1)
   20   CONTINUE
C
C       Upper co-diagonal
C
        DO 30 I=1,N-1,1
          D(I) = H2(I+1)*H1(I)*H1(I+1)
   30   CONTINUE
        D(N) = H2(1)*H1(0)*H1(1)
C
C       Lower co-diagonal
C
        DO 40 I=1,N-1,1
          B(I) = H(I) - H2(I)*H1(I)*H3(I-1) -H2(I+1)*H1(I)*H3(I)
   40   CONTINUE
        B(N) = H(0) - H2(0)*H1(0)*H3(N-1) - H2(1)*H1(0)*H3(0)
C
C       Main diagonal
C
        DO 50 I=1,N-1,1
          K = I-1
          A(I) = 2.0D0*(H(K)+H(I)) + H2(K)*C(K) +
     +           H2(I)*H3(K)*H3(K) + H2(I+1)*C(I)
   50   CONTINUE
        A(N) = 2.0D0*(H(N-1)+H(N)) + H2(N-1)*C(N-1) +
     +         H2(N)*H3(N-1)*H3(N-1) + H2(1)*C(0)
      ENDIF
C
C-----Compute right hand side--------------------------------
C
      DUMMY1 = (FN(1) - FN(0)) * H1(0)
      DO 60 I=1,N-1,1
        DUMMY2 = (FN(I+1) - FN(I)) * H1(I)
        RS(I)  = 3.0D0*(DUMMY2 - DUMMY1)
        DUMMY1 = DUMMY2
   60 CONTINUE
      RS(N) = 3.0D0*((FN(1)-FN(0))*H1(0) - DUMMY1)
C
C-----Find coefficients C(1) to C(N) by------------
C     solving the linear system ...
C
      IF (MREP .EQ. 1) THEN
C
C       ... form LR decomposition
C
        CALL NCYFSY (N,A(1),B(1),D(1),RS,C(1),WORK(1),WORK(N+1),
     +               WORK(2*N+1),WORK(3*N+1),WORK(4*N-3),
     +               WORK(5*N-6),WORK(6*N-6),WORK(7*N-6),
     +               WORK(8*N-10),IERR)
        IF (IERR .NE. 0) RETURN
      ELSE
C
C       ... without decomposing anew
C
        IERR = 0
        CALL NCYFSS (N,RS,C(1),WORK(1),WORK(N+1),WORK(2*N+1),
     +               WORK(3*N+1),WORK(4*N-3),WORK(5*N-6),
     +               WORK(6*N-6),WORK(7*N-6),WORK(8*N-10))
      ENDIF
C
C-----Compute remaining spline coefficients---------------
C
      C(0) = C(N)
C
      A(0) = FN(0) - H2(0)/3.0D0*H1(0)*(C(1)-C(0)) +
     +       H2(N)/3.0D0*H1(N-1)*(C(N)-C(N-1))
      DO 70 I=1,N-1,1
        A(I) = FN(I) - H2(I)/3.0D0*(C(I-1)*H1(I-1) -
     +         H3(I-1)*C(I) + C(I+1)*H1(I))
   70 CONTINUE
      A(N) = A(0)
C
      DO 80 I=0,N-1,1
        B(I) = H1(I)*(A(I+1)-A(I)) -
     +         H(I)/3.0D0*(C(I+1)+2.0D0*C(I))
        D(I) = H1(I)/3.0D0*(C(I+1)-C(I))
   80 CONTINUE
      RETURN
      END
