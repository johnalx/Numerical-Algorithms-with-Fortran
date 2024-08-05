C[BA*)
C[KA{F 6}{Systems of Nonlinear Equations}
C[       {Systems of Nonlinear Equations}*)
C[FE{F 6.2.1.2}{Damped Newton Method for Systems}
C[             {Damped Newton Method for Systems}*)
C[LE*)
      SUBROUTINE SMNEWD (FX,N,MAXIT,IERR,KMAX,LUN,IUPD,
C[IX{SMNEWD}*)
     1                   EPS,RNORM2,F,X,DF,LDDF,IWORK,WORK)
C
C*****************************************************************
C                                                                *
C     SMNEWD finds a solution of a nonlinear system of equations *
C                                                                *
C                   F1(X(1),...,X(N))=0                          *
C                   F2(X(1),...,X(N))=0                          *
C                   - - - - - - - - - -                          *
C                   FN(X(1),...,X(N))=0                          *
C                                                                *
C     by the damped Newton method, if this method converges for  *
C     the starting vector. Here the Jacobi matrix is replaced    *
C     by the forward difference quotients.                       *
C     Hence the user does not need to supply the partial         *
C     derivatives. The parameter IUPD can be used to predetermine*
C     after how many iterations the Jakobi matrix is to be re-   *
C     computed and LR decomposed.                                *
C[BE*)
C     IUPD = 1 specifies the damped Newton method.               *
C     IUPD > 1 specifies a damped version of the simplified      *
C              Newton method.                                    *
C     In general, more iteration steps are required if IUPD > 1. *
C     However, since the function then does not have to be       *
C     evaluated and the Jacobi matrix not LR factored as often,  *
C     this may save computational time.                          *
C                                                                *
C     Three break-off criteria are used:                         *
C     1.  maximum number of iterations has been reached          *
C     2.  euclidean norm of the difference between the old and   *
C         the new approximate solutions is smaller or equal to   *
C         the preset (relative) accuracy bound EPS               *
C     3.  euclidean norm of the function value at the new        *
C         approximate solution is smaller or equal to EPS        *
C                                                                *
C     If desired, there will be output of intermediate results   *
C     (via input parameter LUN).                                 *
C                                                                *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C FX      : SUBROUTINE, that has to be specified by the user. It *
C           describes the system of equations to be solved.      *
C           In the calling program FX has to be defined as       *
C           EXTERNAL. It has to be in the form:                  *
C               SUBROUTINE FX (N,X,F)                            *
C               DOUBLE PRECISION X(N),F(N)                       *
C               F(1) = F1 (X(1),...,X(N))                        *
C               F(2) = F2 (X(1),...,X(N))                        *
C               - - - - - - - - - - - - -                        *
C               F(N) = FN (X(1),...,X(N))                        *
C               RETURN                                           *
C               END                                              *
C N       : number of equations and number of unknowns           *
C LUN     : > 0, file number onto which the iterates are stored  *
C           = 0, no output                                       *
C MAXIT   : maximum number of iterations to be executed          *
C KMAX    : a bound for the damping factor >= 0;                 *
C           KMAX = 0  ==> standard Newton method is used,        *
C           usually a value between 4 and 6 is chosen for KMAX   *
C IUPD    : after each IUPD steps the Jacobi matrix is           *
C           reconfigured and LR decomposed anew. In general, the *
C           method will not converge if IUPD is chosen too large.*
C           IUPD between 1 and 4 usually will give meaningful    *
C           results.                                             *
C EPS     : error parameter                                      *
C X       : N-vector X(1:N); starting vector                     *
C DF      : 2-dim. array DF(1:LDDF,1:N); Jacobi matrix           *
C           (providing storage space)                            *
C LDDF    : leading dimension of DF as defined in the calling    *
C           program. LDDF >= N                                   *
C IWORK   : N-vector IWORK(1:N); auxiliary vector, the pivot     *
C           vector for solving the linear system of equations    *
C WORK    : (4N)-vector WORK(1:4*N); auxiliary vector for X,     *
C           -DELTA X, the old iterate X and the scaling factors  *
C           in GAUSSP or the functional values at X              *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C MAXIT   : number of iterations executed                        *
C IERR    : error parameter                                      *
C           IERR=0 : no error                                    *
C           IERR=1 : error bound was not reached after MAXIT     *
C                    steps                                       *
C           IERR=2 : error when solving the linear system of     *
C                    equations (matrix singular)                 *
C           IERR=3 : incorrect input parameter                   *
C RNORM2  : accuracy estimate                                    *
C               RNORM2 = MIN (XNORM2,FNORM2)                     *
C           (compare with the description of the local variables)*
C X       : N-vector X(1:N); approximate solution                *
C F       : N-vector F(1:N); function value at the approximate   *
C           solution                                             *
C                                                                *
C                                                                *
C LOCAL VARIABLES:                                               *
C ================                                               *
C MARK    : error parameter from GAUSSP                          *
C IX      : starting index for X in the auxiliary vector WORK    *
C IDELTA  : starting index for the stepsize DELTA X in WORK      *
C IXOLD   : starting index for the old approximate solution X    *
C           in WORK                                              *
C IFGAUS  : starting index for the scaling factors from GAUSSP   *
C           in WORK                                              *
C IF      : starting index for the vector of functional values   *
C           at X. In WORK this storage space is identical and    *
C           shared with the one for the scaling factors in       *
C           GAUSSP.                                              *
C XANRM2  : euclidean norm of F(X)                               *
C XNNRM2  : euclidean norm of F(X + 1/2**K * DELTA X)            *
C XNNRMH  : set equal to XNNRM2; however, this constant it is not*
C           altered during damping and thus it can be used later *
C           instead of XNNRM2 if damping will not be used. Thus  *
C           it equals XNNRM2 if K=0.                             *
C XNORM2  : relative accuracy                                    *
C FNORM2  : euclidean norm of the function value at the new      *
C           approximation                                        *
C K       : damping loop counter. The damping factor is 1/2**K   *
C IT      : Newton-iteration loop counter                        *
C I       : control variable                                     *
C IUP     : counter for estimating the Jacobi matrix             *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: GAUSSP, GAUSSS, FENORM, FDIFQU, MACHPD  *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Thomas Eul                                         *
C  date     : 07.02.1985                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER IWORK(N)
      DOUBLE PRECISION F(N), DF(LDDF,N), X(N), WORK(4*N)
      EXTERNAL FX
C
C*****************************************************************
C*               checking the input parameters                   *
C*****************************************************************
C
      IF (MAXIT .GT. 0 .AND. KMAX .GE. 0 .AND. EPS .GT. 0.0D0 .AND.
     1    N .GT. 0 .AND. LUN .GE. 0 .AND. IUPD .GT. 0 .AND.
     1    LDDF .GE. N) THEN
C
C*****************************************************************
C*                        initialization                         *
C*****************************************************************
C
         IERR   = 1
         IX     = 1
         IDELTA = N + 1
         IXOLD  = 2*N + 1
         IFGAUS = 3*N + 1
         IF     = 3*N + 1
         IUP    = IUPD - 1
         K      = 0
         IT     = 0
C
C*****************************************************************
C*                      Newton iteration                         *
C*****************************************************************
C
         IF (LUN .GT. 0) THEN
            I = 1
            WRITE (LUN,1000)
            WRITE (LUN,1100) IT,I,X(1)
            WRITE (LUN,1200) (I,X(I),I=2,N)
         ENDIF
C
C        determining the functional value and its euclidean norm
C        at the starting point
C
         CALL FX (N,X,F)
         FNORM2 = FENORM (N,F)
  100    CONTINUE
            IT = IT + 1
C
C           determining the functional value at the new
C           approximate solution. However, this only involves
C           relabelling if we have used damping, since the norms
C           necessary for the damped algorithm automatically
C           presuppose knowledge of the functional value
C
            IF (K .GT. 0) THEN
               DO 10 I = 1,N
                 F(I) = WORK(IF+I-1)
   10          CONTINUE
            ENDIF
C
C           the euclidean norm of the function value at the new
C           approximation, FNORM2, becomes the euclidean norm of
C           the functional value for the old approximate solution
C
            XANRM2 = FNORM2
C
C*****************************************************************
C*                  solving of DF * DELTA X = F                  *
C*****************************************************************
C
C           if necessary estimation and multiplication
C           of the Jacobi matrix
C
            IUP = IUP + 1
            IF (IUP .EQ. IUPD) THEN
               IUP = 0
               CALL FDIFQU (FX,N,X,F,WORK(IF),DF,LDDF)
C
C              1. LR factorization
C
               CALL GAUSSP (N,DF,LDDF,IWORK,MARK,WORK(IFGAUS))
            ENDIF
C
C           checking for singularity
C
            IF (MARK .NE. 0) THEN
C
C              2. solving the linear system of equations with
C                 the LR factors from GAUSSP
C
               CALL GAUSSS (N,DF,LDDF,IWORK,F,WORK(IDELTA))
C
C*****************************************************************
C*  iteration step without damping, saving of the old X, and     *
C*  determining the norm of F(X + DELTA X)                       *
C*****************************************************************
C
               DO 20 I = 1,N
                  WORK(IX   +I-1) = X(I)
                  WORK(IXOLD+I-1) = X(I)
                  X(I) = X(I) - WORK(IDELTA+I-1)
   20          CONTINUE
C
C              determination of the norm of F(X + DELTA X)
C
               CALL FX (N,X,F)
               XNNRM2 = FENORM (N,F)
               XNNRMH = XNNRM2
C
C*****************************************************************
C*                          damping                              *
C*****************************************************************
C
               K = 0
  200          IF (K .EQ. KMAX .OR. XANRM2 .GT. XNNRM2) GOTO 300
                  K = K + 1
C
C                 Newton step with damping
C
                  DO 30 I = 1,N
                     WORK(IDELTA+I-1) = 0.5D0 * WORK(IDELTA+I-1)
                     WORK(IX+I-1) = WORK(IXOLD+I-1) - WORK(IDELTA+I-1)
   30             CONTINUE
C
C                 determining the norm of F(X + 1/2**K * DELTA X)
C
                  CALL FX (N,WORK(IX),WORK(IF))
                  XNNRM2 = FENORM (N,WORK(IF))
               GOTO 200
  300          CONTINUE
C
C              if XANRM2 > XNNRM2, computations are continued
C              with damping; K=0 indicates that there was
C              no damping
C
               IF (XANRM2 .GT. XNNRM2 .AND. K .GT. 0) THEN
                  DO 40 I = 1,N
                     X(I) = WORK(IX+I-1)
   40             CONTINUE
               ELSE
                  K = 0
               ENDIF
C
C*****************************************************************
C*         test for accuracy and if warranted : stop             *
C*****************************************************************
C
C              1. checking the second break-off criterion,
C                 i.e., determining XNORM2
C
               DO 50 I = 1,N
                  WORK(IX+I-1) = X(I) - WORK(IXOLD+I-1)
   50          CONTINUE
               XNORM2 = FENORM (N,X)
               IF (XNORM2 .GT. 0) THEN
                  XNORM2 = FENORM (N,WORK(IX)) / XNORM2
               ELSE
                  XNORM2 = FENORM (N,WORK(IX))
               ENDIF
C
C              2. checking the third break-off criterion,
C                 i.e., determining FNORM2
C
               IF (K .GT. 0) THEN
                  FNORM2 = XNNRM2
               ELSE
                  FNORM2 = XNNRMH
               ENDIF
C
C              we make use of the smallest of the two error estimates
C
               IF (XNORM2 .GT. FNORM2) THEN
                  RNORM2 = FNORM2
               ELSE
                  RNORM2 = XNORM2
               ENDIF
               IF (LUN .GT. 0) THEN
                  I = 1
                  WRITE (LUN,1100) IT,I,X(1),RNORM2,K
                  WRITE (LUN,1200) (I,X(I),I=2,N)
               ENDIF
               IF (RNORM2 .LE. EPS) THEN
                  MAXIT = IT
                  IERR  = 0
               ENDIF
            ELSE
               IERR  = 2
               MAXIT = IT
            ENDIF
         IF (IT .LT. MAXIT .AND. IERR .EQ. 1) GOTO 100
      ELSE
         IERR  = 3
         MAXIT = 0
      ENDIF
C
C*****************************************************************
C*                    output formats                             *
C*****************************************************************
C
 1000 FORMAT (1X,'ITERATION STEP',10X,'APROXIMATION',14X,
     1           'ERROR ESTIMATE  K')
 1100 FORMAT (1X,I6,2X,I6,8X,D22.15,6X,D22.15,2X,I3)
 1200 FORMAT (1X, 6X,2X,I6,8X,D22.15)
C
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE FDIFQU (FX,N,X,F,FXJPH,DF,LDDF)
C[IX{FDIFQU}*)
C
C*****************************************************************
C                                                                *
C Approximates the  Jacobi matrix for the function               *
C                                                                *
C                  F1 (X(1),...,X(N))                            *
C                  F2 (X(1),...,X(N))                            *
C                  - - - - - - - - - -                           *
C                  FN (X(1),...,X(N))                            *
C                                                                *
C by forward difference quotients.                               *
C[BE*)
C                                                                *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C FX      : SUBROUTINE, that has to be provided by the user.     *
C           It finds the functional value at X. The program      *
C           computes an approximate 1st derivative of FX at X.   *
C           In the calling program FX has to be defined as       *
C           EXTERNAL of the form:                                *
C               SUBROUTINE FX (N,X,F)                            *
C               DOUBLE PRECISION X(N),F(N)                       *
C               F(1) = F1 (X(1),...,X(N))                        *
C               F(2) = F2 (X(1),...,X(N))                        *
C               - - - - - - - - - - - - -                        *
C               F(N) = FN (X(1),...,X(N))                        *
C               RETURN                                           *
C               END                                              *
C N       : number of component functions and variables of FX    *
C X       : N-vector X(1:N); location where the 1st derivative   *
C           is to be estimated                                   *
C F       : N-vector F(1:N); functional value of FX at X         *
C FXJPH   : N-vector FXJPH(1:N); provides storage space for the  *
C           function values
C               F1(X(1),...,X(J-1),X(J)+H,X(J+1),...,X(N))       *
C               F2(X(1),...,X(J-1),X(J)+H,X(J+1),...,X(N))       *
C               - - - - - - - - - - - - - - - - - - - - - -      *
C               FN(X(1),...,X(J-1),X(J)+H,X(J+1),...,X(N))       *
C LDDF    : leading dimension of DF as defined in the calling    *
C           program. LDDF >= N                                   *
C                                                                *
C                                                                *
C OUTPUT PARAMETER:                                              *
C =================                                              *
C DF      : 2-dimensional array DF(1:LDDF,1:N); the computed     *
C           Jacobi matrix at X                                   *
C                                                                *
C                                                                *
C LOCAL VARIABLES:                                               *
C ================                                               *
C H       : stepsize for the forward difference quotients, H is  *
C           taken as the square root of the machine constant     *
C HBEST   : logical variable for determining H. Initially HBEST  *
C           is set to .TRUE., so that on the initial call of     *
C           FDIFQU the constant H is determined. Then HBEST is   *
C           set to .FALSE., and H and HBEST are stored unchanged *
C           by further calls of FDIFQU                           *
C HAEPSM  : auxiliary variable for determining the machine       *
C           constant                                             *
C EPSM    : machine constant                                     *
C XJ      : auxiliary variable used for forming the difference   *
C           quotient with respect to X(J)                        *
C I       : control variable                                     *
C J       : control variable                                     *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: MACHPD                                  *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Thomas Eul                                         *
C  date     : 02.07.1985                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     declarations
C
      INTEGER N, LDDF
      DOUBLE PRECISION X(N), F(N), FXJPH(N), DF(LDDF,N)
C
C     local variables
C
      INTEGER I, J
      LOGICAL HBEST
C
      SAVE H, HBEST
C
      DATA HBEST/.TRUE./
C
      IF (HBEST) THEN
C
C*****************************************************************
C*          determining the square root of the machine constant  *
C*****************************************************************
C
C        EPSM represents the smallest positive number for which
C        (1.0+EPSM) .GT. 1.0 . EPSM is determined as a power of 1./2.
C
         HAEPSM = 1.0D0
   10    CONTINUE
            HAEPSM = 0.5D0 * HAEPSM
         IF (MACHPD(1.0D0+HAEPSM) .EQ. 1) GOTO 10
         EPSM = 2.0D0 * HAEPSM
         H = DSQRT(EPSM)
         HBEST = .FALSE.
      ENDIF
C
      DO 30 J = 1,N
         XJ = X(J)
         X(J) = XJ + H
         CALL FX (N,X,FXJPH)
         X(J) = XJ
         DO 20 I = 1,N
            DF(I,J) = (FXJPH(I)-F(I))/H
   20    CONTINUE
   30 CONTINUE
      RETURN
      END
