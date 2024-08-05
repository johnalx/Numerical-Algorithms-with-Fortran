C[BA*)
C[LE*)
      SUBROUTINE NLESYS (FX,M,N,JACMA,DFX,DF,LDDF,MAXIT,EPS,
C[IX{NLESYS}*)
     +                  KMAX,LUN,X,F,RNORM2,IERR,D,S,WORK)
C
C*****************************************************************
C                                                                *
C  SUBROUTINE NLESYS solves a nonlinear system of equations      *
C  consisting of M+1 equations in N+1 unknowns by applying the   *
C  damped Newton method provided that it converges for the       *
C  starting vector.                                              *
C[BE*)
C  The Jacobi matrix must be either determined by a subroutine   *
C  provided by the user or approximated by central difference    *
C  quotients.                                                    *
C  Three break-off criteria are used:                            *
C    1. maximum number of iterations has been reached            *
C    2. the relative difference between the euclidean norms of   *
C       the old and new approximations is smaller or equal to    *
C       EPS                                                      *
C    3. the euclidean norm of the function value for the newest  *
C       approximation is smaller than or equal to EPS            *
C                                                                *
C  this program was developed by using the SUBROUTINES SMNEWD    *
C  and SMNEWT by Thomas Eul and by SUBROUTINE SNLFIT by          *
C  Ilona Westermann.                                             *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  ================+                                             *
C  FX     : Function SUBROUTINE, that has to be provided by the  *
C           user. It defines the system of equations that are to *
C           be solved. In the calling program FX has to be       *
C           defined as EXTERNAL. It must have the form:          *
C               SUBROUTINE FX (X, N, F, M)                       *
C               INTEGER M, N                                     *
C               DOUBLE PRECISION X(0:N), F(0:M)                  *
C               ------------------                               *
C               F(0) = F0 (X(0), ... , X(N)                      *
C                .                                               *
C                .                                               *
C               F(M) = FM (X(0), ... , X(N)                      *
C               ------------------                               *
C               RETURN                                           *
C               END                                              *
C  M      : M+1 = number of equations                            *
C  N      : N+1 = number of unknowns                             *
C  JACMA  : logical parameter:                                   *
C           If JACMA = .TRUE., the user has provided his own     *
C           subroutine for determining the Jacobi matrix.        *
C           Otherwise the Jacobi matrix is computed intrinsically*
C           by central difference quotients.                     *
C  DFX    : A SUBROUTINE that has to be provided by the user if  *
C           JACMA = .TRUE.  DFX determines the Jacobi matrix.    *
C           In the calling program DFX has to be defined as      *
C           EXTERNAL and should have the following form:         *
C               SUBROUTINE DFX (X, M, N, DF, LDDF)               *
C               INTEGER M, N, LDDF                               *
C               DOUBLE PRECISION X(0:N), DF(0:LDDF,0:N)          *
C               ----------------------                           *
C               DF(I,K) = partial derivative of function I       *
C                         by X(K) at the position X              *
C                         K=0(1)N, I=0(1)M                       *
C               ----------------------                           *
C               RETURN                                           *
C               END                                              *
C  DF     : 2-dim. array DF(0:LDDF,0:N+1); used for storage      *
C           space for the Jacobi matrix                          *
C  LDDF   : leading dimension of DF as defined in the calling    *
C           program                                              *
C  X      : (N+1)-vector X(0:N) containing the starting vector   *
C  MAXIT  : maximum number of iterations to be executed          *
C  EPS    : relative error bound                                 *
C  KMAX   : damping bound, if KMAX=0 => standard Newton method   *
C  LUN    : > 0, file number onto which the iteration steps are  *
C                stored                                          *
C           = 0, no output                                       *
C                                                                *
C                                                                *
C  AUXILIARY PARAMETERS:                                         *
C  =====================                                         *
C  D      : (N+1)-vector D(0:N); auxiliary vector for SUBROUTINE *
C                                SHOUSE                          *
C  S      : (N+1)-vector S(0:N); for the Newton direction        *
C  WORK   : (N+1)-vector WORK(0:M); auxiliary vector for         *
C                                   SUBROUTINE JACOBI            *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  X      : (N+1)-vector X(0:N) containing the approximate       *
C           solution                                             *
C  F      : (N+1)-vector F(0:N) containing the function values   *
C           at the approximate solution                          *
C  IERR   : error parameter:                                     *
C           = 0  : everything o.k.                               *
C           = 1  : error in the input parameters                 *
C           = 2  : error while solving the linear minimization   *
C                  problem (matrix singular)                     *
C           = 3  : error bound not attained after MAXIT          *
C                  iterations                                    *
C  RNORM2 : an estimate for the achieved accuracy,               *
C           RNORM2 = MIN (XNORM2, FNORM2)                        *
C           (compare with the description of local variables)    *
C                                                                *
C                                                                *
C  LOCAL VARIABLES:                                              *
C  ================                                              *
C  FNORM2 : euclidean norm of F(X)                               *
C  FNORMN : euclidean norm of F(X + 1/2**K * DELTA X)            *
C  XNORM2 : relative accuracy                                    *
C  EPSMA  : machine constant                                     *
C  SMOOTH : damping factor 1/2**K                                *
C  MARK   : error parameter of SHOUSE                            *
C  K      : counter for the damping loop. The damping factor is  *
C           1/2**K                                               *
C  IT     : Newton iteration loop counter                        *
C  I      : control variable                                     *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: JACOBI, SHOUSE, SENORM, MACHPD          *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Ilona Westermann                                   *
C  date     : 09.01.1987                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER M, N, LDDF, MAXIT, KMAX, LUN, IERR
      LOGICAL JACMA
      DIMENSION DF(0:LDDF,0:N+1), X(0:N), F(0:M), D(0:N), S(0:N),
     +          WORK(0:M)
      EXTERNAL FX
C
C  test of the input parameters
C
      IF (M .GE. N .AND. LDDF .GE. M .AND. MAXIT .GT. 0 .AND.
     +    KMAX .GE. 0 .AND. EPS .GT. 0.0D0 .AND. N .GE. 0 .AND.
     +    M .GE. 0 .AND. LUN .GE. 0) THEN
         IERR = 3
         IT = 0
C
C  determine machine constant in case the Jacobi matrix is calculated in
C  this program by central difference quotients
C
         IF (.NOT. JACMA) THEN
            EPSMA = 1.0D0
   10       EPSMA = EPSMA * 0.5D0
            IF (MACHPD(EPSMA + 1.0D0) .EQ. 1) GOTO 10
            EPSMA = EPSMA * 2.0D0
         ENDIF
         IF (LUN .GT. 0) THEN
            I = 0
            WRITE (LUN,1000)
            WRITE (LUN,1100) IT,I,X(0)
            WRITE (LUN,1200) (I,X(I),I=1,N)
            WRITE (LUN,1300)
         ENDIF
C
C  determine functional value and its norm at the starting point
C
         CALL FX (X, N, F, M)
         FNORM2 = DSQRT(SENORM (F, M))
C
C  Newton- teration
C
  100    IT = IT + 1
C
C  determining the Jacobi matrix
C
         IF (JACMA) THEN
            CALL DFX (X, M, N, DF, LDDF)
         ELSE
            CALL JACOBI (FX, X, M, N, DF, LDDF, EPSMA, WORK)
         ENDIF
C
C  minimize (DF * DELTA X + F)      (DELTA X  here is called S)
C
         DO 15 I=0,M
            DF(I,N+1) = - F(I)
   15    CONTINUE
         CALL SHOUSE (DF, LDDF, M, N, D, S, MARK)
C
C  check for nonsingularity
C
         IF (MARK .NE. 1) THEN
C
C  damped Newton step
C
            DO 40 K=0,KMAX
               SMOOTH = 1.0D0 / 2.0D0 ** K
               DO 20 I=0,N
                  D(I) = X(I) + S(I) * SMOOTH
   20          CONTINUE
               CALL FX (D, N, F, M)
               FNORMN = DSQRT (SENORM (F, M))
C
C  if FNORMN < FNORM2 continue calculations with damping
C
               IF (FNORMN .LT. FNORM2) THEN
                  DO 30 I=0,N
                     S(I) = S(I) * SMOOTH
                     X(I) = D(I)
   30             CONTINUE
                  GOTO 60
               ENDIF
   40       CONTINUE
C
C  if the loop was completed, i.e., damping gave no improvement,
C  continue calculations without damping
C
            DO 50 I=0,N
               X(I) = X(I) + S(I)
   50       CONTINUE
            CALL FX (X, N, F, M)
            FNORMN = DSQRT (SENORM (F, M))
            K = 0
   60       CONTINUE
C
C  check accuracy and if necessary stop:
C
C  1. use relative error estimate for X
C
            XNORM2 = DSQRT (SENORM (X, N))
            IF (XNORM2 .EQ. 0.0D0) THEN
               XNORM2 = DSQRT (SENORM (S, N))
            ELSE
               XNORM2 = DSQRT (SENORM (S, N)) / XNORM2
            ENDIF
C
C  2. find the norm of the functional value for X
C
            FNORM2 = FNORMN
C
C     select the smaller one of both norms
C
            RNORM2 = DMIN1 (FNORM2, XNORM2)
            IF (LUN .GT. 0) THEN
               I = 0
               WRITE (LUN,1100) IT,I,X(0),RNORM2,K
               WRITE (LUN,1200) (I,X(I),I=1,N)
               WRITE (LUN,1300)
            ENDIF
            IF (RNORM2 .LE. EPS) THEN
               IERR = 0
               MAXIT = IT
            ENDIF
         ELSE
            IERR = 2
            MAXIT = IT
         ENDIF
         IF (IT .LT. MAXIT .AND. IERR .EQ. 3) GOTO 100
      ELSE
         IERR = 1
         MAXIT = 0
      ENDIF
      RETURN
 1000 FORMAT (1X,'ITERATION NUMBER',10X,'APROXIMATION',14X,
     1           'ACCURACY ESTIMATE   K')
 1100 FORMAT (1X,I6,2X,I6,8X,D22.15,6X,D22.15,2X,I3)
 1200 FORMAT (1X, 6X,2X,I6,8X,D22.15)
 1300 FORMAT (1X)
      END
