C[BA*)
C[KA{F 7}{Eigenvalues and Eigenvectors of Matrices}
C[       {Eigenvalues and Eigenvectors of Matrices}*)
C[FE{F 7.3.1}
C[  {Vector Iteration for the Dominant Eigenvalue}
C[  {Vector Iteration for the Dominant Eigenvalue and the
C[   Associated Eigenvector of a Matrix}*)
C[LE*)
      SUBROUTINE EIGVAL(A,N,LDA,M0,EPSI,X,Y,Z,EV,IERR)
C[IX{EIGVAL}*)
C
C*****************************************************************
C*                                                               *
C*  Determining the eigenvalue of largest magnitude of an n by n *
C*  matrix A with the corresponding eigenvector by vector        *
C*  iteration                                                    *
C[BE*)
C*                                                               *
C*                                                               *
C*   INPUT PARAMETERS:                                           *
C*   =================                                           *
C*   A   : 2-dimensional DOUBLE PRECISION array A(1:LDA,1:N);    *
C*         the input matrix                                      *
C*   N   : order of A                                            *
C*   LDA : leading dimension of A as defined in the calling      *
C*         program                                               *
C*   M0  : maximum iteration number                              *
C*   EPSI: desired relative accuracy                             *
C*         (larger than 1E-12) (DOUBLE PRECISION)                *
C*   X   : N-vector X(1:n) in DOUBLE PRECISION                   *
C*                                                               *
C*                                                               *
C*   OUTPUT PARAMETERS:                                          *
C*   ==================                                          *
C*   Y   : N-vector Y(1:n) in DOUBLE PRECISION; the eigenvector  *
C*   Z   : N-vector Z(1:n) in DOUBLE PRECISION; the residual     *
C*         vector A * Y - EV * Y                                 *
C*   EV  : the dominant eigenvalue in DOUBLE PRECISION           *
C*   IERR: error parameter:                                      *
C*         =0: run was successfully completed                    *
C*         =1: maximum number of iterations was reached, i.e.,   *
C*             eigenvalue/vector is complex or the problem is    *
C*             poorly conditioned                                *
C*                                                               *
C*---------------------------------------------------------------*
C*                                                               *
C*  subroutines required: DBNORM, MAVE, QUOT, QSCAL              *
C*                                                               *
C*                                                               *
C*                                                               *
C*****************************************************************
C*                                                               *
C*  author   : Juergen Beckmann                                  *
C*  date     : 10.24.1985                                        *
C*  source   : FORTRAN 77                                        *
C*                                                               *
C[BA*)
C*****************************************************************
C[BE*)
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(LDA,LDA),Z(LDA),X(LDA),Y(LDA)
        IERR = 0
        DO 100 I=1,N
           Y(I) = 1.0D0
 100    CONTINUE
        CALL DBNORM(Y,N)
        EV = 0.0D0
        ITER  = 0
 200    ITER  = ITER + 1
        EM    = EV
        DO 300 I=1,N
           X(I) = Y(I)
 300    CONTINUE
        CALL MAVE(A,X,Y,N,LDA)
        EV = QUOT(X,Y,N)
        CALL DBNORM(Y,N)
        IF (ITER .EQ. 1) THEN
           GOTO 200
        END IF
        DO 400 I = 1,N
           Z(I) = Y(I) - X(I)
 400    CONTINUE
        CALL QSCAL(Z,Z,S,N)
        S = DSQRT(S)
        IF (ITER .EQ. M0) THEN
           IERR = 1
           GOTO 500
        END IF
        IF(S .GT. EPSI .OR. DABS(EV-EM) .GT. EPSI) THEN
           GOTO 200
        END IF
 500    RETURN
        END
C
C
C[BA*)
C[LE*)
        SUBROUTINE MAVE(A,X,Y,N,LDA)
C[IX{MAVE}*)
C
C*****************************************************************
C*                                                               *
C*  compute  Y = A * X                                           *
C*                                                               *
C*****************************************************************
C[BE*)
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(LDA,LDA),X(LDA),Y(LDA)
        DO 200 I=1,N
           Y(I) = 0.0D0
           DO 100 J=1,N
              Y(I) = Y(I) + A(I,J) * X(J)
 100       CONTINUE
 200    CONTINUE
        RETURN
        END
C
C
C[BA*)
C[LE*)
        SUBROUTINE DBNORM(X,N)
C[IX{DBNORM}*)
C
C*****************************************************************
C*                                                               *
C*  Normalizes a DOUBLE PRECISION vector to euclidean length 1   *
C*                                                               *
C*****************************************************************
C[BE*)
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION X(N)
        CALL QSCAL(X,X,S,N)
        S = DSQRT(S)
        IF (S .NE. 0.0D0) THEN
           DO 100 I=1,N
              X(I) = X(I) / S
 100       CONTINUE
        END IF
        RETURN
        END
C
C
C[BA*)
C[LE*)
        DOUBLE PRECISION FUNCTION QUOT(X,Y,N)
C[IX{QUOT}*)
C[BE*)
C
C*****************************************************************
C*                                                               *
C*  auxiliary routine for EIGVAL                                 *
C*                                                               *
C*****************************************************************
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION X(N),Y(N)
        QUOT = 1.0D0
        S    = 0.0D0
        N1   = 0
        DO 100 I=1,N
           IF (X(I) .NE. 0.0D0) THEN
              S  = S + Y(I) / X(I)
              N1 = N1 + 1
           END IF
 100    CONTINUE
        IF (N1 .NE. 0) THEN
           QUOT = S / DBLE(N1)
        END IF
        RETURN
        END
C
C
C[BA*)
C[LE*)
        SUBROUTINE QSCAL(X,Y,R,N)
C[IX{QSCAL}*)
C
C*****************************************************************
C*                                                               *
C*  computes the dot product of two DOUBLE PRECISION vectors     *
C*                                                               *
C*****************************************************************
C[BE*)
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION X(N),Y(N)
        R = 0.0D0
        DO 100 I=1,N
           R = R + X(I)*Y(I)
 100    CONTINUE
        RETURN
        END
