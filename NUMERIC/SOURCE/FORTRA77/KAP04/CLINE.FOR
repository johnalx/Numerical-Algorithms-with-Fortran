C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 4.15.2}{Condition Estimates}{Condition Estimates}*)
C[LE*)
C[FE{}
C[  {Condition Estimate according to Cline}
C[  {Condition Estimate according to Cline}*)
C[LE*)
      SUBROUTINE CLINE (L,R,N,IA,COND,X,Y,Z,ZSUM,NA)
C[IX{CLINE}*)
C
C*****************************************************************
C                                                                *
C   This subroutine computes an estimate for the condition number*
C   of a matrix A whose  LR decoposition is known, where L is a  *
C   unit lower triangular and R is a nonsingular upper triangular*
C   matrix.                                                      *
C[BE*)
C                                                                *
C                                                                *
C   INPUT PARAMETERS:                                            *
C   =================                                            *
C   L    : array L(1:IA,1:N) containing the unit diagonal lower  *
C          triangular matrix of the LR decomposition of A        *
C   R    : array R(1:IA,1:N) containing the nonsingular upper    *
C          triangular matrix of the LR decomposition of A        *
C   N    : order of the matrices L and R                         *
C   IA   : leading dimension of the matrices L and R, as         *
C          stipulated in the calling program                     *
C                                                                *
C                                                                *
C   OUTPUT PARAMETER:                                            *
C   =================                                            *
C   COND : Estimate for the condition number of A                *
C                                                                *
C                                                                *
C   AUXILIARY PARAMETERS:                                        *
C   =====================                                        *
C   X    : ] N-vectors ..(1:N)                                   *
C   Y    : ]                                                     *
C   Z    : ]                                                     *
C   ZSUM : (N+1)-vector ZSUM(0:N)                                *
C   NA   : array NA(1:IA,1:N)                                    *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C   Required subroutines: BACK, MXNORM, ZNORM                    *
C                                                                *
C*****************************************************************
C                                                                *
C   Author    : Michaela Kisters                                 *
C   Date      : 12.09.1990                                       *
C   Source    : FORTRAN 77                                       *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      DOUBLE PRECISION  NA(IA,N), L(IA,N), R(IA,N), COND, X(N),
     F                  Y(N), Z(N), ZSUM(0:N), MXNORM, ZNORM,
     F                  KUNEND, SMI, SPL, V
C
C   For R transpose (=TRANS(R)) we determine X=X(I) with X(I)=+ or -1 and
C   Y=Y(I)=INV(TRANS(R))*X for I=1, ..., N, so that the 1-norm of Y becomes
C   as large as possible.
C
      X(1) = 1.0D0
      Y(1) = 1.0D0/R(1,1)
      DO 10 I=2,N
         Y(I) = -R(1,I)*Y(1)/R(I,I)
   10 CONTINUE
      DO 20 K=2,N
         V = 1.0D0/R(K,K)
         X(K) = Y(K) - V
         Y(K) = Y(K) + V
         SMI = DABS(X(K))
         SPL = DABS(Y(K))
         DO 30 I=K+1,N
            V = R(K,I)/R(I,I)
            X(I) = Y(I) - V*X(K)
            Y(I) = Y(I) - V*Y(K)
            SMI = SMI + DABS(X(I))
            SPL = SPL + DABS(Y(I))
   30    CONTINUE
         IF (SMI .GT. SPL) THEN
            DO 40 I=K,N
               Y(I) = X(I)
   40       CONTINUE
            X(K) = -1.0D0
         ELSE
            X(K) = 1.0D0
         ENDIF
   20 CONTINUE
C
C   Use backsubstitution to find Z with
C   TRANS(L)*Z=Y.
C
      CALL BACK (L,Y,N,IA,Z)
C
C   Estimate  KUNEND, the row sum norm of A
C
      KUNEND = MXNORM(Z,N)/MXNORM(X,N)
C
C   Estimate COND(A)
C
      COND = ZNORM(L,R,N,IA,ZSUM,NA)*KUNEND
c
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE BACK (L,B,N,IL,X)
C[IX{BACK}*)
C
C*****************************************************************
C                                                                *
C   Solving a triangular system TRANS(L)*X=B by backsubstitution *
C[BE*)
C                                                                *
C   INPUT PARAMETERS:                                            *
C   =================                                            *
C   L    : array L(1:IL,1:N) contaning the entries of the unit   *
C          diagonal lower tringular matrix L                     *
C   B    : N-vector B(1:N), the right hand side of the system    *
C   N    : order of the matrix L                                 *
C   IL   : leading dimension of L, as stipulated in the calling  *
C          program                                               *
C                                                                *
C                                                                *
C   OUTPUT PARAMETER:                                            *
C   =================                                            *
C   X    : N-vector X(1:N), the solution vector                  *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C   Required subroutines: none                                   *
C                                                                *
C*****************************************************************
C                                                                *
C   Author    : Michaela Kisters                                 *
C   Date      : 09.12.1990                                       *
C   Sourcee   : FORTRAN 77                                       *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      DOUBLE PRECISION  L(IL,N), B(N), X(N), SUM
c
      X(N) = B(N)/L(N,N)
      DO 10 K=N-1,1,-1
         SUM = 0.0D0
         DO 20 J=K+1,N
            SUM = SUM + L(J,K)*X(J)
   20    CONTINUE
         X(K) = 1.0D0/L(K,K)*(B(K) - SUM)
   10 CONTINUE
c
      RETURN
      END
C
C
C[BA*)
C[LE*)
      FUNCTION MXNORM(X,N)
C[IX{MXNORM}*)
C
C*****************************************************************
C                                                                *
C   Calculate the maximum norm MXNORM of a vector X              *
C[BE*)
C                                                                *
C   INPUT PARAMETERS:                                            *
C   =================                                            *
C   X    : N-vector X(1:N)                                       *
C   N    : dimension of X                                        *
C                                                                *
C                                                                *
C   OUTPUT PARAMETER:                                            *
C   =================                                            *
C   MXNORM : Maximum norm of X                                   *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C   Required subroutines: none                                   *
C                                                                *
C*****************************************************************
C                                                                *
C   Author    : Michaela Kisters                                 *
C   Date      : 12.09.1990                                       *
C   Source    : FORTRAN 77                                       *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      DOUBLE PRECISION  X(N), MXNORM
c
      MXNORM = DABS(X(1))
      DO 10 I=1,N
         IF (DABS(X(I)) .GT. MXNORM) MXNORM = DABS(X(I))
   10 CONTINUE
c
      RETURN
      END
C
C
C[BA*)
C[LE*)
      FUNCTION ZNORM(L,R,N,IA,ZSUM,NA)
C[IX{ZNORM}*)
C
C*****************************************************************
C                                                                *
C   Computes the row sum norm of a matrix A which is given by its*
C   LR decomposition by a unit diagonal lower triangular matrix L*
C   and a nonsingular upper tringular matrix R.                  *
C[BE*)
C                                                                *
C                                                                *
C   INPUT PARAMETERS:                                            *
C   =================                                            *
C   L    : array L(1:IA,1:N) containing the unit diagonal lower  *
C          triangular matrix of the LR decomposition of A        *
C   R    : array R(1:IA,1:N) containing the nonsingular upper    *
C          triangular matrix of the LR decomposition of A        *
C   N    : order of the matrices L and R                         *
C   IA   : leading dimension of the matrices L and R, as         *
C          stipulated in the calling program                     *
C                                                                *
C                                                                *
C   OUTPUT PARAMETER:                                            *
C   =================                                            *
C   ZNORM : Row sum norm of the matrix A                         *
C                                                                *
C                                                                *
C   AUXILIARY PARAMETERS:                                        *
C   =====================                                        *
C   ZSUM : (N+1)-vector ZSUM(0:N)                                *
C   NA   : array NA(1:IA,1:N)                                    *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C   Required subroutines: none                                   *
C                                                                *
C*****************************************************************
C                                                                *
C   Author    : Michaela Kisters                                 *
C   Date      : 12.09.1990                                       *
C   Source    : FORTRAN 77                                       *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      DOUBLE PRECISION  ZSUM(0:N), ZNORM, L(IA,N), R(IA,N), NA(IA,N)
c
      ZSUM(0) = -99.0D0
      DO 10 I=1,N
         ZSUM(I) = 0.0D0
         DO 20 J=1,N
            NA(I,J) = 0.0D0
            DO 30 K=1,I
               NA(I,J) = NA(I,J) + L(I,K)*R(K,J)
   30       CONTINUE
            ZSUM(I) = ZSUM(I) + DABS(NA(I,J))
   20    CONTINUE
         IF (ZSUM(I) .GT. ZSUM(I-1)) ZNORM = ZSUM(I)
   10  CONTINUE
c
       RETURN
       END
