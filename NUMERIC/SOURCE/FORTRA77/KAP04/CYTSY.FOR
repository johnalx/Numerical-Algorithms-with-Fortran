C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 4.11.2}
C[  {Systems with Symm.\ Cycl.\ Trid.\ Str.\ Nonsing.\ Matrices}
C[  {Systems with Symmetric Cyclically
C[   Tridiagonal Strongly Nonsingular Matrices}*)
C[LE*)
      SUBROUTINE CYTSY (N,DM,DU,CR,RS,X,MARK)
C[IX{CYTSY}*)
C
C*****************************************************************
C                                                                *
C     Solving a system of linear  equations                      *
C                  A * X = RS                                    *
C     for a cyclically tridiagonal, symmetric, positive definite *
C     matrix A.                                                  *
C[BE*)
C     The matrix A is given by two N-vectors DM and DU.          *
C     The system of equations has the form:                      *
C                                                                *
C     DM(1) * X(1) + DU(1) * X(2) + DU(N) * X(N)       = RS(1)   *
C                                                                *
C     DU(I-1) * X(I-1) + DM(I) * X(I) + DU(I) * X(I+1) = RS(I)   *
C                for I = 2, ..., N - 1, and                      *
C                                                                *
C     DU(N) * X(1) + DU(N-1) * X(N-1) + DM(N) * X(N)   = RS(N)   *
C                                                                *
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     N    : number of equations, N > 2                          *
C     DM   : N-vector DM(1:N); main diagonal of A                *
C            DM(1), DM(2), ... , DM(N)                           *
C     DU   : N-vector DU(1:N); upper co-diagonal of A            *
C            DU(1), DU(2), ... , DU(N-1); the off-diagonal       *
C            element A(1,N) is stored in DU(N).                  *
C     RS   : N-vector RS(1:N); the right hand side               *
C                                                                *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     DM   :)                                                    *
C     DU   :) overwritten with intermediate vectors              *
C     CR   :)                                                    *
C     RS   :)                                                    *
C     X    : N-vector X(1:N), containing the solution            *
C     MARK : error parameter                                     *
C            MARK=-2 : condition N > 2 is not satisfied          *
C            MARK=-1 : A is not positive definite                *
C            MARK= 0 : numerically the matrix A is not strongly  *
C                      nonsingular                               *
C            MARK= 1 : ok                                        *
C                                                                *
C     NOTE : If MARK =  1, the determinant of A is given as:     *
C               DET A = DM(1) * DM(2) * ... * DM(N)              *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: CYTSYP, CYTSYS, MACHPD                  *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Gisela Engeln-Muellges                             *
C  date     : 27.04.1988                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION DM(1:N),DU(1:N),CR(1:N),RS(1:N),X(1:N)
      MARK = -2
      IF (N .LT. 3) RETURN
C
C  factorization of the matrix A
C
      CALL CYTSYP (N,DM,DU,CR,MARK)
C
C  if MARK = 1, update and backsubstitute
C
      IF (MARK .EQ. 1) THEN
         CALL CYTSYS (N,DM,DU,CR,RS,X)
      ENDIF
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE CYTSYP (N,DM,DU,CR,MARK)
C[IX{CYTSYP}*)
C
C*****************************************************************
C                                                                *
C     Factoring a cyclically tridiagonal, symmetric and positive *
C     definite matrix A, that is given by the two N-vectors DM   *
C     and DU, into its Cholesky factors                          *
C                    A = R(TRANSP) * D * R                       *
C     by applying the Cholesky-method for tridiagonal cyclic     *
C     matrices.                                                  *
C[BE*)
C     The form of the system of equations is identical           *
C     to the one described in SUBROUTINE CYTSY.                  *
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     N    : number of equations, N > 2                          *
C     DM   : N-vector DM(1:N); main diagonal of A                *
C            DM(1), DM(2), ... , DM(N)                           *
C     DU   : N-vector DU(1:N); upper co-diagonal of A            *
C            DU(1), DU(2), ... , DU(N-1); the off-diagonal       *
C            element A(1,N) is stored in DU(N).                  *
C            Due to symmetry the lower co-diagonal does not need *
C            to be stored separately.                            *
C                                                                *
C                                                                *
C     OUTPUT PARAMETER:                                          *
C     =================                                          *
C     DM   :) overwritten with auxiliary vectors from the        *
C     DU   :) factorization of A. The co-diagonal of the unit    *
C     CR   :) upper tridiagonal matrix R is stored in DU, the    *
C             diagonal matrix D appears in DM and the right hand *
C             side in CR.                                        *
C     MARK : error parameter                                     *
C            MARK=-2 : condition N > 2 is not satisfied          *
C            MARK=-1 : A is not positive definite                *
C            MARK= 0 : numerically A is not strongly             *
C                      nonsingular                               *
C            MARK= 1 : ok                                        *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutine required: MACHPD                                   *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Gisela Engeln-Muellges                             *
C  date     : 27.04.1988                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION DM(1:N),DU(1:N),CR(1:N)
C
C   calculating the machine constant
C
      FMACHP = 1.0D0
   10 FMACHP = 0.5D0 * FMACHP
      IF (MACHPD(1.0D0+FMACHP) .EQ. 1) GOTO 10
      FMACHP = FMACHP * 2.0D0
C
C   determinaing the relative error bound
C
      EPS = 4.0D0 * FMACHP
C
C   testing of condition N > 2
C
      MARK = -2
      IF (N .LT. 3) RETURN
C
C   checking for positive definite matrix A and for strong
C   nonsingularity of A for N=1 and if A is
C
      ROW = DABS(DM(1)) + DABS(DU(1)) + DABS(DU(N))
      IF (ROW .EQ. 0.0D0) THEN
         MARK = 0
         RETURN
      END IF
      D = 1.0D0/ROW
      IF (DM(1) .LT. 0.0D0) THEN
         MARK = -1
         RETURN
      ELSEIF (DABS(DM(1))*D .LE. EPS) THEN
         MARK = 0
         RETURN
      END IF
C
C   factoring A while checking for a positive definite and strong
C   nonsingular matrix A
C
      DUMMY = DU(1)
      DU(1) = DU(1)/DM(1)
      CR(1) = DU(N)/DM(1)
      DO 20 I=2,N-1,1
         ROW = DABS(DM(I)) + DABS(DU(I)) + DABS(DUMMY)
         IF (ROW .EQ. 0.0D0) THEN
            MARK = 0
            RETURN
         END IF
         D = 1.0D0/ROW
         DM(I) = DM(I) - DUMMY * DU(I-1)
         IF (DM(I) .LT. 0.0D0) THEN
            MARK = -1
            RETURN
         ELSEIF (DABS(DM(I))*D .LE. EPS) THEN
            MARK = 0
            RETURN
         ENDIF
         IF (I .LT. (N-1)) THEN
            CR(I) = -DUMMY * CR(I-1)/DM(I)
            DUMMY = DU(I)
            DU(I) = DU(I)/DM(I)
         ELSE
            DUMMY2 = DU(I)
            DU(I) = (DU(I) - DUMMY * CR(I-1))/DM(I)
         ENDIF
   20 CONTINUE
      ROW = DABS(DU(N)) + DABS(DM(N)) + DABS(DUMMY2)
      IF (ROW .EQ. 0.0D0) THEN
         MARK = 0
         RETURN
      END IF
      D = 1.0D0/ROW
      DM(N) = DM(N) - DM(N-1) * DU(N-1) * DU(N-1)
      DUMMY = 0.0D0
      DO 30 I=1,N-2,1
         DUMMY = DUMMY + DM(I) * CR(I) * CR(I)
   30 CONTINUE
      DM(N) = DM(N) - DUMMY
      IF (DM(N) .LT. 0) THEN
         MARK = -1
         RETURN
      ELSEIF (DABS(DM(N))*D .LE. EPS) THEN
         MARK = 0
         RETURN
      ENDIF
      MARK = 1
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE CYTSYS (N,DM,DU,CR,RS,X)
C[IX{CYTSYS}*)
C
C*****************************************************************
C                                                                *
C     Solving a system of linear equations                       *
C                  A * X = RS                                    *
C     for a cyclically tridiagonal, symmetric, positive definite *
C     matrix A, once its Cholesky factors have been calculated   *
C     in SUBROUTINE CYTSYP.                                      *
C[BE*)
C     Here the factoring matrices are used                       *
C     as input arrays. They are stored in the three N-vectors    *
C     DM, CR and DU.                                             *
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     N    : number of equations, N > 2                          *
C     DM   : N-vector DM(1:N); ) the vectors DU, DM and CR       *
C     DU   : N-vector DU(1:N); ) contain the factors of A        *
C     CR   : N-vector CR(1:N); )                                 *
C     RS   : N-vector RS(1:N); the right hand side               *
C                                                                *
C                                                                *
C     OUTPUT PARAMETER:                                          *
C     =================                                          *
C     X    : N-vector X(1:N) containing the solution of the      *
C            system of equations                                 *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Gisela Engeln-Muellges                             *
C  date     : 27.04.1988                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION DM(1:N),DU(1:N),CR(1:N),RS(1:N),X(1:N)
C
C  updating
C
      DUMMY = RS(1)
      RS(1) = DUMMY/DM(1)
      SUM = CR(1)*DUMMY
      DO 10 I=2,N-1,1
         DUMMY = RS(I)-DU(I-1)*DUMMY
         RS(I) = DUMMY/DM(I)
         IF (I .NE. (N-1)) SUM = SUM + CR(I)*DUMMY
   10 CONTINUE
      DUMMY = RS(N)-DU(N-1)*DUMMY
      DUMMY = DUMMY-SUM
      RS(N) = DUMMY/DM(N)
C
C  backsubstitution
C
      X(N) = RS(N)
      X(N-1) = RS(N-1)-DU(N-1)*X(N)
      DO 30 I=N-2,1,-1
         X(I) = RS(I)-DU(I)*X(I+1)-CR(I)*X(N)
   30 CONTINUE
      RETURN
      END
