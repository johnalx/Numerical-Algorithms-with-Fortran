C[BA*)
C[LE*)
      SUBROUTINE TRDNPD (N,DM,DU,RS,X,MARK)
C[IX{TRDNPD}*)
C
C*****************************************************************
C                                                                *
C     Solving a linear system of equations                       *
C                  A * X = RS                                    *
C     for a tridiagonal, symmetric, strongly nonsingular         *
C     matrix A.                                                  *
C[BE*)
C     The matrix A is defined by the two N-vectors               *
C     DM and DU. The system of equations is given as follows:    *
C                                                                *
C     DM(1) * X(1) + DU(1) * X(2)                      = RS(1)   *
C                                                                *
C     DU(I-1) * X(I-1) + DM(I) * X(I) + DU(I) * X(I+1) = RS(I)   *
C            for I = 2, ... ,N-1, and                            *
C                                                                *
C     DU(N-1) * X(N-1) + DM(N) * X(N)                  = RS(N)   *
C                                                                *
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     N    : number of equations, N > 2                          *
C     DM   : N-vector DM(1:N); main diagonal of A                *
C            DM(1), DM(2), ... , DM(N)                           *
C     DU   : N-vector DU(1:N); co-diagonal of A                  *
C            DU(1), DU(2), ... , DU(N-1)                         *
C     RS   : N-vector X(1:N); the right hand side                *
C                                                                *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     DM   :)                                                    *
C     DU   :) overwritten with auxiliary vectors                 *
C     RS   :)                                                    *
C     X    : N-vector X(1:N) containing the solution of the      *
C            system of equations                                 *
C     MARK : error parameter                                     *
C            MARK= 1 : A is positive definite.                   *
C            MARK= 0 : numerically the matrix A is not strongly  *
C                      nonsingular                               *
C            MARK=-1 : A is strongly nonsingular, but not        *
C                      positive definite.                        *
C            MARK=-2 : condition N > 2 is not satisfied          *
C                                                                *
C     NOTE: If MARK = 1 or MARK = -1, then the determinant of A  *
C           can be calculated as:                                *
C              DET A = DM(1) * DM(2) * ... * DM(N),              *
C           while the inertia of A is given by the number of     *
C           positive as well as negative entries among the DM(I).*
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: TRDSYP, TRDSYS, MACHPD                  *
C                                                                *
C*****************************************************************
C                                                                *
C  authors    : Gisela Engeln-Muellges                           *
C  date       : 01.03.1992                                       *
C  source     : FORTRAN 77                                       *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION DM(1:N),DU(1:N),RS(1:N),X(1:N)
      MARK = -2
      IF (N .LT. 3) RETURN
C
C  Factoring A
C
      CALL TNPSYP (N,DM,DU,MARK)
C
C  if MARK = 1 or MARK = -1, update and backsubstitute
C
      IF ((MARK .EQ. 1) .OR. (MARK .EQ. -1)) THEN
         CALL TRDSYS (N,DM,DU,RS,X)
      ENDIF
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE TNPSYP (N,DM,DU,MARK)
C[IX{TNPSYP}*)
C
C*****************************************************************
C                                                                *
C     Factoring a tridiagonal, symmetric, and strongly           *
C     nonsingular matrix A, that is given by the two N-vectors DM*
C     and DU, into the product A = R(TRANSP) * D * R  for a unit *
C     upper triangular matrix R by applying the root-free        *
C     Cholesky-method for tridiagonal matrices. The form of the  *
C     system matrix A is identical with the one desribed in      *
C     SUBROUTINE TRDSY.                                          *
C[BE*)
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     N    : number of equations, N > 2                          *
C     DM   : N-vector DM(1:N); main diagonal of A                *
C            DM(1), DM(2), ... , DM(N)                           *
C     DU   : N-vector DU(1:N); co-diagonal of A                  *
C            DU(1), DU(2), ... , DU(N-1);                        *
C            due to symmetry of A its lower and upper            *
C            co-diagonals coincide                               *
C                                                                *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     DM   :) overwritten with auxiliary vectors containing the  *
C     DU   :) factors of A. The co-diagonal of the unit upper    *
C             triangular bidiagonal matrix R is stored in DU,    *
C             while the diagonal matrix D is stored in DM.       *
C     MARK : error parameter                                     *
C            MARK= 1 : A is strongly nonsingular and positive    *
C                      definite.                                 *
C            MARK= 0 : numerically the matrix A is not strongly  *
C                      nonsingular.                              *
C            MARK=-1 : A is strongly nonsingular, but not        *
C                      positive definite.                        *
C            MARK=-2 : condition N > 2 is not met.               *
C     Note : If Mark = +/- 1, then the inertia of A, i.e., the   *
C            number of positive and negative eigenvalues of A,   *
C            is the same as the number of positive and negative  *
C            numbers among the components of DM.                 *
C                                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: MACHPD                                  *
C                                                                *
C*****************************************************************
C                                                                *
C  authors    : Gisela Engeln-Muellges                           *
C  date       : 25.04.1988                                       *
C  source     : FORTRAN 77                                       *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION DM(1:N),DU(1:N)
C
C   calculating the machine constant
C
      FMACHP = 1.0D0
   10 FMACHP = 0.5D0 * FMACHP
      IF (MACHPD(1.0D0+FMACHP) .EQ. 1) GOTO 10
      FMACHP = FMACHP * 2.0D0
C
C   determining the relative error bound
C
      EPS = 4.0D0 * FMACHP
C
C   checking whether N > 2
C
      MARK = -2
      IF (N .LT. 3) RETURN
      MARK = 1
      DU(N) = 0.0D0
C
C   testing for strong nonsingularity of the matrix A for N=1
C
      ROW = DABS(DM(1)) + DABS(DU(1))
      IF (ROW .EQ. 0.0D0) THEN
         MARK = 0
         RETURN
      ENDIF
      D = 1.0D0/ROW
      IF (DM(1) .LT. 0.0D0) MARK = -1
      IF (DABS(DM(1))*D .LE. EPS) THEN
         MARK = 0
         RETURN
      ENDIF
C
C   factoring A while checking for strong nonsingularity
C
      DUMMY = DU(1)
      DU(1) = DU(1)/DM(1)
      DO 20 I=2,N,1
         ROW = (DABS (DM(I)) + DABS(DU(I)) + DABS(DUMMY))
         IF (ROW .EQ. 0.0D0) THEN
            MARK = 0
            RETURN
         ENDIF
         D = 1.0D0/ROW
         DM(I) = DM(I) - DUMMY * DU(I-1)
         IF (DM(I) .LT. 0.0D0) MARK = -1
         IF (DABS(DM(I))*D .LE. EPS) THEN
            MARK = 0
            RETURN
         ENDIF
         IF (I .LT. N) THEN
            DUMMY = DU(I)
            DU(I) = DU(I)/DM(I)
         ENDIF
   20 CONTINUE
      RETURN
      END
