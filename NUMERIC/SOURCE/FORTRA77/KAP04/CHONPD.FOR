C[BA*)
C[LE*)
      SUBROUTINE CHONPD (N,A,LDA,RS,X,Z,MARK)
C[IX{CHONPD}*)
C
C*****************************************************************
C                                                                *
C     Solving a linear system of equations                       *
C                    A * X = RS                                  *
C     for a symmetric, strongly nonsingular system matrix A      *
C     using the root-free Cholesky-method.                       *
C[BE*)
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     N    : order of the square matrix A                        *
C     A    : 2-dimensional array A(1:LDA,1:N) containing the     *
C            symmetric matrix A. It is sufficient to specify     *
C            only the elements of the upper triangle of A        *
C            (including the elements of the main diagonal);      *
C            only the elemente in A's upper triangle will be     *
C            processed.                                          *
C     LDA  : leading dimension of A as defined in the calling    *
C            program                                             *
C     RS   : N-vector RS(1:N) containing the right hand side     *
C                                                                *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     A    : 2-dimensional array A(1:LDA,1:N) which contains     *
C            the Cholesky factors of A = R(TRANSP) * D * R       *
C            for a diagonal D and a unit upper triangular matrix *
C            R. The elements of R, excluding the diagonal ones,  *
C            are stored in the upper triangle of A. The elements *
C            of D appear on the main diagonal of A.              *
C     X    : N-vector X(1:N) containing the solution of the      *
C            system of equations                                 *
C     MARK : error parameter                                     *
C            MARK= 1 : A is positive definite                    *
C            MARK= 0 : numerically the matrix A is not stongly   *
C                      nonsingular                               *
C            MARK=-1 : A is strongly nonsingular, but not        *
C                      positive definite                         *
C                                                                *
C     NOTE : If MARK = +/-1 , then the determinant of A can be   *
C            calculated as follows:                              *
C               DET A = A(1,1) * A(2,2) * ... * A(N,N)           *
C            If Mark = +/-1 , then the inertia of A, i.e., the   *
C            number of positive and negative eigenvalues of A,   *
C            is the same as the number of positive and negative  *
C            entries on the diagonal of the output matrix A.     *
C                                                                *
C                                                                *
C     AUXILIARY PARAMETER:                                       *
C     ====================                                       *
C     Z    : N-vector Z(1:N)                                     *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: CHOKYP, CHOKYS, MACHPD                  *
C                                                                *
C*****************************************************************
C                                                                *
C  authors   : Gisela Engeln-Muellges                            *
C  date      : 01.07.1992                                        *
C  source    : FORTRAN 77                                        *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A(LDA,N),RS(N),X(N),Z(N)
C
C  Cholesky factorization of the matrix A
C
      CALL CNPKYP (N,A,LDA,Z,MARK)
C
C  Updating and backsustitution
C
      IF ((MARK .EQ. 1) .OR. (MARK .EQ. -1)) THEN
         CALL CHOKYS (N,A,LDA,RS,X,Z)
      ENDIF
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE CNPKYP (N,A,LDA,Z,MARK)
C[IX{CNPKYP}*)
C
C*****************************************************************
C                                                                *
C     Factoring a strongly nonsingular, symmetric matrix A into  *
C             A = R(TRANSP) * D * R                              *
C     for a diagonal D and a unit upper triangular R using the   *
C     root-free Cholesky decomposition method.                   *
C[BE*)
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     N    : order of the matrix A                               *
C     A    : 2-dimensional array A(1:LDA,1:N) containing the     *
C            symmetric matrix A. It is sufficient to store only  *
C            the elements of the upper triangle of A including   *
C            the main diagonal due to symmetry.                  *
C     LDA  : leading dimension of A as defined in the calling    *
C            program                                             *
C                                                                *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     A    : 2-dimensional array A(1:LDA,1:N) containing the     *
C            Cholesky Factors R and D in the array A where       *
C            A = R(TRANSP) * D * R. The elements of R except for *
C            the unit diagonal are stored in the upper triangle  *
C            of A. The elements of D are on the main diagonal of *
C            A.                                                  *
C     MARK : error parameter                                     *
C            MARK= 1 : A is positive definite                    *
C            MARK= 0 : numerically the matrix A is not strongly  *
C                      nonsingular                               *
C            MARK=-1 : A is strongly nonsingular, but not        *
C                      positive definite                         *
C                                                                *
C                                                                *
C     AUXILIARY PARAMETER:                                       *
C     ====================                                       *
C     Z    : N-vector Z(1:N)                                     *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: MACHPD                                  *
C                                                                *
C*****************************************************************
C                                                                *
C  authors   : Gisela Engeln-Muellges                            *
C  date      : 01.07.1992                                        *
C  source    : FORTRAN 77                                        *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A(LDA,N),Z(N)
      MARK=1
C
C  calculating the machine constant
C
      FMACHP = 1.0D0
   10 FMACHP = 0.5D0 * FMACHP
      IF (MACHPD(1.0D0+FMACHP) .EQ. 1) GOTO 10
      FMACHP = FMACHP * 2.0D0
C
C  determining the relative error bound
C
      EPS = 4.0D0 * FMACHP
C
C  calculating the absolute row sums of the matrix A,
C  with storage in the auxiliary vector Z
C
      DO 20 I=1,N
         S = 0.0D0
         DO 30 K=I,N
            S = S + DABS(A(I,K))
   30    CONTINUE
         DO 40 K=1,I-1
            S = S + DABS(A(K,I))
   40    CONTINUE
         IF (S .EQ. 0.0D0) THEN
            MARK = 0
            RETURN
         ENDIF
         Z(I) = 1.0D0/S
   20 CONTINUE
C
C  Factoring the matrix A
C
      DO 50 J=1,N
         DO 60 I=1,J-1
            H = A(I,J)
            A(I,J) = H/A(I,I)
            DO 60 K=I+1,J
               A(K,J) = A(K,J) - H * A(I,K)
   60    CONTINUE
         IF (A(J,J) .LT. 0.0D0) MARK = -1
         IF (DABS(A(J,J))*Z(J) .LE. EPS) THEN
            MARK = 0
            RETURN
         END IF
   50 CONTINUE
      RETURN
      END
