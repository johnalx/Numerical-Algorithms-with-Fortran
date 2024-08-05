C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 4.10.1}
C[  {Systems with Tridiagonal Matrices}
C[  {Systems with Tridiagonal Matrices}*)
C[LE*)
      SUBROUTINE TRDIG (N,DL,DM,DU,RS,X,MARK)
C[IX{TRDIG}*)
C
C*****************************************************************
C                                                                *
C     This routine solves a linear system of equations           *
C                  A * X = RS                                    *
C     for a tridiagonal, strongly nonsingular matrix A.          *
C[BE*)
C     The matrix is given by the three vectors DL, DM and DU     *
C     which designate the lower co-diagonal, the diagonal and    *
C     the upper co-diagonal elements of A, respectively.         *
C     The system of equations has the form :                     *
C                                                                *
C     DM(1) * X(1) + DU(1) * X(2)                      = RS(1)   *
C                                                                *
C     DL(I) * X(I-1) + DM(I) * X(I) + DU(I) * X(I+1)   = RS(I)   *
C            for I = 2, ..., N-1, and                            *
C                                                                *
C     DL(N) * X(N-1) + DM(N) * X(N)                    = RS(N)   *
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     N    : number of equations, N > 2                          *
C     DL   : N-vector DL(1:N); lower co-diagonal of A            *
C            DL(2), DL(3), ... ,DL(N)                            *
C     DM   : N-vector DM(1:N); the diagonal of A                 *
C            DM(1), DM(2), ... , DM(N)                           *
C     DU   : N-vector DU(1:N); upper co-diagonal of A            *
C            DU(1), DU(2), ... , DU(N-1)                         *
C     RS   : N-vector RS(1:N); the right hand side               *
C                                                                *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     DL   :)                                                    *
C     DM   :)                                                    *
C     DU   :) these are overwritten with auxiliary vectors       *
C     RS   :)                                                    *
C     X    : N-vector X(1:N), containing the solution of the     *
C            system of equations                                 *
C     MARK : error parameter                                     *
C            MARK= 1 : everything is o.k.                        *
C            MARK= 0 : the matrix A is not strongly nonsingular  *
C            MARK=-1 : error on N: N <= 2                        *
C                                                                *
C     NOTE: if MARK = 1, the determinant of A can be calculated  *
C           after this subroutine has run as follows:            *
C              DET A = DM(1) * DM(2) * ... * DM(N)               *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: TRDIGP, TRDIGS, MACHPD                  *
C                                                                *
C*****************************************************************
C                                                                *
C  author     : Gisela Engeln-Muellges                           *
C  date       : 05.02.1988                                       *
C  source     : FORTRAN 77                                       *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION DL(1:N),DM(1:N),DU(1:N),RS(1:N),X(1:N)
      MARK = -1
      IF (N .LT. 3) RETURN
C
C  Factoring the matrix A
C
      CALL TRDIGP (N,DL,DM,DU,MARK)
C
C  If MARK = 1, update the right hand side and solve via backsubstitution
C
      IF (MARK .EQ. 1) THEN
         CALL TRDIGS (N,DL,DM,DU,RS,X)
      ENDIF
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE TRDIGP (N,DL,DM,DU,MARK)
C[IX{TRDIGP}*)
C
C*****************************************************************
C                                                                *
C     Factor the tridiagonal, strongly nonsingular matrix A,     *
C     that is given by the vectors DL, DM and DU of co-diagonal  *
C     and diagonal entries into the product  L * R  for a bi-    *
C     diagonal lower triangular matrix  L  and a unit bidiagonal *
C     upper triangular matrix R. The form of the system matrix A *
C     is identical to that given in SUBROUTINE TRDIG.            *
C[BE*)
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     N    : number of equations, N > 2                          *
C     DL   : N-vector DL(1:N); lower co-diagonal of A            *
C            DL(2), DL(3), ... ,DL(N)                            *
C     DM   : N-vector DM(1:N); the diagonal of A                 *
C            DM(1), DM(2), ... , DM(N)                           *
C     DU   : N-vector DU(1:N); upper co-diagonal of A            *
C            DU(1), DU(2), ... , DU(N-1)                         *
C                                                                *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     DL   : the lower co-diagonal remains unchanged, since it   *
C            forms the lower co-diagonal of L                    *
C     DM   :) overwritten with auxiliary vectors that contain    *
C     DU   :) the factors of A: The co-diagonal of the unit      *
C             bidiagonal upper tridiagonal matrix R is stored in *
C             DU. The diagonal elements of L are stored in DM    *
C     MARK : error parameter                                     *
C            MARK= 1 : everything is o.k.                        *
C            MARK= 0 : the matrix is not strongly nonsingular    *
C            MARK=-1 : condition N > 2 is not fulfilled          *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  required subroutine: MACHPD                                   *
C                                                                *
C*****************************************************************
C                                                                *
C  author     : Gisela Engeln-Muellges                           *
C  date       : 05.02.1988                                       *
C  source     : FORTRAN 77                                       *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION DL(1:N),DM(1:N),DU(1:N)
C
C   testing whether N > 2
C
      MARK = -1
      IF (N .LT. 3) RETURN
C
C   calculate the machine constant
C
      FMACHP = 1.0D0
   10 FMACHP = 0.5D0 * FMACHP
      IF (MACHPD(1.0D0+FMACHP) .EQ. 1) GOTO 10
      FMACHP = FMACHP * 2.0D0
C
C   determine the relative error bound
C
      EPS = 4.0D0 * FMACHP
C
C   checking for strong nonsingularity with N=1
C
      ROW = DABS(DM(1)) + DABS(DU(1))
      IF (ROW .EQ. 0.0D0) THEN
         MARK = 0
         RETURN
      ENDIF
      D = 1.0D0/ROW
      IF (DABS(DM(1))*D .LE. EPS) THEN
         MARK = 0
         RETURN
      ENDIF
C
C   factoring A while checking for strong nonsingularity
C
      DL(1) = 0.0D0
      DU(N) = 0.0D0
      DU(1) = DU(1)/DM(1)
      DO 20 I=2,N,1
         ROW = DABS(DL(I)) + DABS(DM(I)) + DABS(DU(I))
         IF (ROW .EQ. 0.0D0) THEN
            MARK = 0
            RETURN
         ENDIF
         D = 1.0D0/ROW
         DM(I) = DM(I) - DL(I) * DU(I-1)
         IF (DABS(DM(I))*D .LE. EPS) THEN
            MARK = 0
            RETURN
         ENDIF
         IF (I .LT. N) THEN
            DU(I) = DU(I)/DM(I)
         ENDIF
   20 CONTINUE
      MARK=1
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE TRDIGS (N,DL,DM,DU,RS,X)
C[IX{TRDIGS}*)
C
C*****************************************************************
C                                                                *
C     Solving a linear system of equations                       *
C                  A * X = RS                                    *
C     for a tridiagonal, strongly nonsingular matrix A, once the *
C     tringular factorization has been calculated by SUBROUTINE  *
C     TRDIGP. The factors are used as known input matrices and   *
C     are stored in the three N-vectors DL, DM and DU.           *
C[BE*)
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     N    : number of equations, N > 2                          *
C     DL   : N-vector DL(1:N); ) the vectors DL, DM and DU con-  *
C     DM   : N-vector DM(1:N); ) tain the factors of the matrix  *
C     DU   : N-vector DU(1:N); ) A obtained as output of TRDIGP  *
C     RS   : N-vector RS(1:N); right hand side                   *
C                                                                *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     X    : N-vector X(1:N); containing the solution of the     *
C            system of equations                                 *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  author     : Gisela Engeln-Muellges                           *
C  date       : 05.02.1988                                       *
C  source     : FORTRAN 77                                       *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION DL(1:N),DM(1:N),DU(1:N),RS(1:N),X(1:N)
C
C  updating
C
      RS(1) = RS(1)/DM(1)
      DO 10 I=2,N,1
         RS(I) = (RS(I) - DL(I) * RS(I-1)) / DM(I)
   10 CONTINUE
C
C  backsubstitution
C
      X(N) = RS(N)
      DO 20 I=N-1,1,-1
         X(I) = RS(I) - DU(I) * X(I+1)
   20 CONTINUE
      RETURN
      END
