C[BA*)
C[KA{F 4}
C[  {Direct Methods for Solving Linear Systems}
C[  {Direct Methods for Solving Systems of Linear Equations}*)
C[FE{F 4.5.1}
C[  {Gau"s Algorithm with Column Pivot Search}
C[  {Gau"s Algorithm with Column Pivot Search}*)
C[LE*)
      SUBROUTINE GAUSS(N,A,LDA,Y,X,MARK,D,IPIVOT)
C[IX{GAUSS}*)
C
C*****************************************************************
C                                                                *
C  Solving a linear system of equations  A * X = Y  by applying  *
C  the Gauss-elimination method with scaling and column pivot    *
C  search.                                                       *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  N        : order of the linear system.                        *
C  A        : 2-dimensional array A(1:LDA,1:N); the matrix A is  *
C             the system matrix of the equations, (A = A(ORG)).  *
C  LDA      : leading dimension of A as defined in the calling   *
C             program.                                           *
C  Y        : N-vector Y(1:N); the right hand side of the system *
C             of equations.                                      *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  A        : 2-dimensional array A(1:LDA,1:N), containing the   *
C             factors L and R with P * A(ORG) = L * R.           *
C             P = permutation matrix, L = unit lower triangular  *
C             matrix, and R = upper triangular matrix.           *
C  X        : N-vector X(1:N); the solution vector of the system *
C             of equations.                                      *
C  MARK     : = 1, even number of row permutations.              *
C             =-1, odd number of row permutations.               *
C             = 0, input array A is numerically singular.        *
C             The determinant of A can be computed as :          *
C                DET(A(ORG)) = MARK * A(1,1) * ... * A(N,N).     *
C  D        : N-vector D(1:N); the reciprocals of the row sum    *
C             norms of A(ORG) that serve as scaling factors:     *
C             D(I) = 1./(ABS(A(I,1)) + ... + ABS(A(I,N)))  for   *
C             I = 1, ..., N.                                     *
C  IPIVOT   : N-vector IPIVOT(1:N); it indicates the row         *
C             permutations for the scaled column pivot search    *
C             and thereby defines the permutation matrix P.      *
C             If e.g. IPIVOT(2) = 7, then the 7th row of A(ORG)  *
C             is permuted to become the 2nd row of P * A(ORG).   *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: GAUSSP, GAUSSS, MACHPD                  *
C                                                                *
C*****************************************************************
C                                                                *
C  authors   : Gisela Engeln-Muellges, Guido Dubois              *
C  date      : 04.25.88                                          *
C  source    : FORTRAN 77                                        *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1:LDA,1:N),Y(1:N),X(1:N),D(1:N)
      INTEGER IPIVOT(1:N)
C
C  Factor the matrix A by using SUBROUTINE GAUSSP.
C
      CALL GAUSSP(N,A,LDA,IPIVOT,MARK,D)
C
C  Updating and backsubstitution via SUBROUTINE GAUSSS
C  in order to find the solution of the system of equations.
C
      IF(MARK .NE. 0) CALL GAUSSS(N,A,LDA,IPIVOT,Y,X)
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE GAUSSP(N,A,LDA,IPIVOT,MARK,D)
C[IX{GAUSSP}*)
C
C*****************************************************************
C                                                                *
C  Factoring the matrix A into the product of two matrices L and *
C  R so that  P * A = L * R, where P = permutation matrix,       *
C  L = unit lower triangular matrix and R = upper triangular     *
C  matrix by applying the Gauss-elimination method with          *
C  scaling and column pivot search.                              *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  N        : order of the system of equations.                  *
C  A        : 2-dimensional array A(1:LDA,1:N); the system matrix*
C             of the system of equations, (A = A(ORG)).          *
C  LDA      : leading dimension of A as defined in the calling   *
C             program.                                           *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  A        : 2-dimensional array A(1:LDA,1:N), containing the   *
C             factors L and R with  P * A(ORG) = L * R  for      *
C             P a permutation matrix. The upper triangular R     *
C             is stored in the upper triangle of A. The unit     *
C             lower triangular matrix L, except for the diagonal *
C             ones, is stored in the lower triangle of A.        *
C  IPIVOT   : N-vector IPIVOT(1:N); it indicates the row         *
C             permutations of the scaled column pivot search     *
C             algorithm and thus defines the permutation matrix  *
C             P. If e.g. IPIVOT(2) = 7, then the 7th row of      *
C             A(ORG) has become the 2nd row of P * A(ORG).       *
C  MARK     : = 1, even number of row permutations.              *
C             =-1, odd number of row permutations.               *
C             = 0, system matrix A is numerically singular.      *
C             The determinant of A is :                          *
C                DET(A(ORG)) = MARK * A(1,1) * ... * A(N,N).     *
C  D        : N-vector D(1:N); the reciprocals of the row sum    *
C             norms of A(ORG) that serve as scaling factors:     *
C             D(I) = 1./(ABS(A(I,1)) + ... + ABS(A(I,N)))  for   *
C             I = 1, ..., N.                                     *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: MACHPD                                  *
C                                                                *
C*****************************************************************
C                                                                *
C  authors   : Gisela Engeln-Muellges, Guido Dubois              *
C  date      : 04.25.88                                          *
C  source    : FORTRAN 77                                        *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1:LDA,1:N),D(1:N)
      INTEGER IPIVOT(1:N)
C
C  Local storage of error parameter RELERR in case that the
C  SUBROUTINE is called repeatedly.
C
      SAVE RELERR,IFLAG
      DATA IFLAG /0/
      MARK=1
C
C  Calculation of the machine constant and initializing the relative
C  error.
C
      IF(IFLAG .EQ. 0) THEN
         IFLAG=1
         FMACHP=1.0D0
   10    FMACHP=0.5D0*FMACHP
         IF(MACHPD(1.0D0+FMACHP) .EQ. 1) GOTO 10
         RELERR=8.0D0*FMACHP
      END IF
C
C  Calculation of row sum norms of A and initializing
C  the PIVOT vector.
C
      DO 20 I=1,N
         IPIVOT(I)=I
         SUM=DABS(A(I,1))
         DO 30 K=2,N
            SUM=SUM+DABS(A(I,K))
   30    CONTINUE
         IF(SUM .EQ. 0.0D0) THEN
            MARK=0
            RETURN
         ELSE
            D(I)=1.0D0/SUM
         END IF
   20 CONTINUE
      IF(N .EQ. 1) RETURN
C
C  Triangular factorization.
C
      DO 40 I=1,N-1
C
C  Determine the pivot row.
C
         PIVOT=DABS(A(I,I))*D(I)
         IPVT=I
         DO 50 J=I+1,N
            DUMMY=DABS(A(J,I))*D(J)
            IF(DUMMY .GT. PIVOT) THEN
               PIVOT=DUMMY
               IPVT=J
            END IF
   50    CONTINUE
         IF(PIVOT .LT. RELERR) THEN
            MARK=0
            RETURN
         ELSE
            IF(IPVT .NE. I) THEN
C
C  Interchange the I-th and the IPVT-th row of A.
C
               MARK=-MARK
               J=IPIVOT(I)
               IPIVOT(I)=IPIVOT(IPVT)
               IPIVOT(IPVT)=J
               DUMMY=D(I)
               D(I)=D(IPVT)
               D(IPVT)=DUMMY
               DO 60 J=1,N
                  DUMMY=A(I,J)
                  A(I,J)=A(IPVT,J)
                  A(IPVT,J)=DUMMY
   60          CONTINUE
            END IF
C
C  Perform the elimination step.
C
            DO 70 J=I+1,N
               A(J,I)=A(J,I)/A(I,I)
               FAK=A(J,I)
               DO 80 K=I+1,N
                  A(J,K)=A(J,K)-FAK*A(I,K)
   80          CONTINUE
   70       CONTINUE
         END IF
   40 CONTINUE
      IF(DABS(A(N,N)) .LT. RELERR) MARK=0
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE GAUSSS(N,A,LDA,IPIVOT,Y,X)
C[IX{GAUSSS}*)
C
C*****************************************************************
C                                                                *
C  Calculating the solution X of a linear system of equations    *
C  A * X = Y, where A has been factored via Gauss-elimination    *
C  in SUBROUTINE GAUSSP.                                         *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  N        : order of the system of equations.                  *
C  A        : 2-dimensional array A(1:LDA,1:N) containing the    *
C             factors L and  R  with  P * A(ORG) = L * R  for    *
C             P a permutation matrix. This array is the output   *
C             matrix of SUBROUTINE GAUSSP.                       *
C  LDA      : leading dimension of A as defined in the calling   *
C             program.                                           *
C  IPIVOT   : N-vector IPIVOT(1:N); it indicates the row         *
C             interchanges in P * A relative to A(ORG). It is an *
C             output of SUBROUTINE GAUSSP.                       *
C  Y        : N-vector Y(1:N); the right hand side of the system *
C             of equations.                                      *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  X        : N-vector X(1:N); the solution vector for the       *
C             system of equations.                               *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  authors   : Gisela Engeln-Muellges, Guido Dubois              *
C  date      : 04.25.88                                          *
C  source    : FORTRAN 77                                        *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(1:LDA,1:N),Y(1:N),X(1:N)
      INTEGER IPIVOT(1:N)
      IF(N .EQ. 1) THEN
         X(1)=Y(1)/A(1,1)
         RETURN
      END IF
C
C  Updating the right hand side.
C
      IPVT=IPIVOT(1)
      X(1)=Y(IPVT)
      DO 10 I=2,N
         SUM=0.0D0
         DO 20 J=1,I-1
            SUM=SUM+A(I,J)*X(J)
   20    CONTINUE
         IPVT=IPIVOT(I)
         X(I)=Y(IPVT)-SUM
   10 CONTINUE
C
C  Compute the solution vector X by backsubstitution.
C
      X(N)=X(N)/A(N,N)
      DO 50 I=N-1,1,-1
         SUM=0.0D0
         DO 40 K=N,I+1,-1
            SUM=SUM+A(I,K)*X(K)
   40    CONTINUE
         X(I)=(X(I)-SUM)/A(I,I)
   50 CONTINUE
      RETURN
      END
