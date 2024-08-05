C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 4.7.2}
C[  {The Conjugate Gradient Method}
C[  {The Conjugate Gradient Method}*)
C[LE*)
      SUBROUTINE CG (A, N, IA, Y, X, IERR, D, G, AMULD)
C[IX{CG}*)
C
C*****************************************************************
C                                                                *
C  This SUBROUTINE solves a linear system of equations AX = Y    *
C  using the conjugate gradient method.                          *
C                                                                *
C  ASSUMPTION:                                                   *
C  ===========                                                   *
C          A must be a symmetric and positive definite N by N    *
C          matrix                                                *
C[BE*)
C                                                                *
C  Input PARAMETERS:                                             *
C  =================                                             *
C  A     : 2-dim. array  A(1:IA, 1:N), containing the NxN system *
C          matrix A. Only the upper triangle of A shall be used  *
C          and we do not check whether A is indeed symmetric.    *
C  N     : order of the system                                   *
C  IA    : leading dimension of A, as specified in the calling   *
C          program                                               *
C  Y     : N-vector Y(1:N), the right hand side                  *
C                                                                *
C  HILFSPARAMETER:                                               *
C  ===============                                               *
C  D     : N-vector D(1:N)                                       *
C  G     : N-vector G(1:N)                                       *
C  AMULD : N-vector AMULD(1:N) containing  A*D                   *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  X     : N-vector X(1:N), the solution of the linear system    *
C  IERR  : error parameter:                                      *
C            = 0, if the denominator of ALPHA vanishes           *
C            = 1, all is ok, the solution has been found         *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Required subroutines: MACHPD                                  *
C                                                                *
C*****************************************************************
C                                                                *
C  Author    : Gisela Engeln-MÅllges                             *
C  Date      : 02.12.1991                                        *
C  Source    : FORTRAN  77                                       *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C  Declarations
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION A(1:IA, 1:N), Y(1:N), X(1:N), D(1:N), G(1:N),
     +                 AMULD(1:N)
C
C  Compute the machine constant
C
      FMACHP = 1.0D0
   10 FMACHP = FMACHP * 0.5D0
      IF(MACHPD(1.0D0+FMACHP) .EQ. 1) GOTO 10
      FMACHP = 8.0D0 * FMACHP
C
C  Initialize the auxiliary vectors  D and G and use
C  the zero vector for a start
C
      DO 20 I=1,N
        HELP = Y(I)
        D(I) = HELP
        G(I) = -HELP
        X(I) = 0.0D0
   20 CONTINUE
C
C  Perform N conjugate gradient steps
C
      DO 30 K=0,N-1
C
C       Initialize numerator and denominator for ALPHA
C
        XNUM = 0.0D0
        DENOM = 0.0D0
C
C       update ALPHA according to:
C       ALPHA = -(D(TRANSP)*G) /A*D(TRANSP)*(A*D))
C
        DO 40 I=1,N
          XNUM = XNUM + D(I) * G(I)
          HELP = 0.0D0
          DO 50 J=1,I-1
            HELP = HELP + A(J,I) * D(J)
   50     CONTINUE
          DO 60 J=I,N
            HELP = HELP + A(I,J) * D(J)
   60     CONTINUE
          AMULD(I) = HELP
          DENOM = DENOM + D(I) * HELP
   40   CONTINUE
C
C       check whether the denominator of  ALPHA  is zero
C
        IF(DABS(DENOM) .LT. FMACHP) THEN
          IERR = 0
          RETURN
        ENDIF
        ALPHA = -XNUM / DENOM
C
C       update  X := X + ALPHA * D
C
        DO 70 I=1,N
          X(I) = X(I) + ALPHA * D(I)
   70   CONTINUE

C
C       update  G := G + ALPHA * A * D
C       and find its norm: NORM;
C       we also check whether X is a good enough approximation
C       of the solution so that computations can be stopped
C       with less than  N  CG-steps.
C
        GNORM = 0.0D0
        DO 80 I=1,N
          G(I) = G(I) + ALPHA * AMULD(I)
          GNORM = GNORM + G(I)*G(I)
   80   CONTINUE
        IF(GNORM .LT. FMACHP) THEN
          IERR = 1
          RETURN
        ENDIF
C
C       Calculate a new  BETA :
C
C       BETA = (G(TRANSP)*(A*D)) / (D(TRANSP)*(A*D))
C
        XNUM = 0.0D0
        DO 90 I=1,N
          XNUM = XNUM + G(I) * AMULD(I)
   90   CONTINUE
        BETA = XNUM / DENOM
C
C       update  D := -G + BETA * D
C
        DO 100 I=1,N
          D(I) = -G(I) + BETA * D(I)
  100   CONTINUE
   30 CONTINUE
      IERR = 1
C
      RETURN
      END
