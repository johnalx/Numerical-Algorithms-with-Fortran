C[BA*)
C[KA{F 5}
C[  {Iterative Methods for Linear Systems}
C[  {Iterative Methods for Linear Systems}*)
C[FE{F 5.4}
C[  {The Gau"s-Seidel Iteration}
C[  {The Gau"s-Seidel Iteration}*)
C[LE*)
      SUBROUTINE ADSOR(A,N,IA,B,X,KADAPT,EPS,KMAX,IMETH,ISWITC,
C[IX{ADSOR}*)
     *                 OMEGA,WORK,RES,ITNUMB,IERR)
C
C*****************************************************************
C                                                                *
C  This program solves an inhomogeneous linear system AX = B of  *
C  equations with a nonsingular system matrix A. The method of   *
C  Jacobi is used jointly with relaxation, where the relaxation  *
C  parameter OMEGA is adjusted during the iteration (adaptive    *
C  SOR method).                                                  *
C[BE*)
C  For a suitable choice of parameters (refer to the remark      *
C  below), this program can perform the Gau·-Seidel method or    *
C  a non-adaptive SOR method.                                    *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  A      : 2-dimensional array A(1:IA,1:N), containing the      *
C           system matrix for the linear equations               *
C  N      : size of the linear system                            *
C  IA     : leading dimension of A, as specified in the calling  *
C           program                                              *
C  B      : N-vector B(1:N), the right hand side of the system   *
C  X      : N-vector X(1:N) containing the starting value for    *
C           iteration                                            *
C  KADAPT : Number of iterations, after which the relaxation     *
C           parameter is to be redefined                         *
C  EPS    : desired accuracy; the iteration is stopped when the  *
C           maximum norm of the relative error does not exceed   *
C           EPS                                                  *
C  KMAX   : Maximal number of iterations allowed                 *
C  IMETH  : parameter that determines the method used:           *
C           = 0, adaptive SOR method                             *
C           = 1, SOR method for a given relaxation parameter     *
C           = 2, Gau·-Seidel method                              *
C  ISWITC : parameter that determines the convergence criterion  *
C           to be used:                                          *
C           = 0, none                                            *
C           = 1, row sum criterion                               *
C           = 2, column sum criterion                            *
C           = 3, criterion of Schmidt and v. Mises               *
C  OMEGA  : in case IMETH=1, the optimal relaxation parameter    *
C           must be part of the input; otherwise only the name   *
C           must be declared in the callimng program.            *
C                                                                *
C                                                                *
C  REMARKS:                                                      *
C  ========                                                      *
C  For the adaptive SOR method (IMETH=0) we recommend to set     *
C  KADAPT=4 or KADAPT=5.                                         *
C  If the optimal relaxationcoefficient Wopt is known for A, then*
C  one should set IMETH=1 and OMEGA = Wopt, i.e., the SOR method *
C  with given optimal relaxation coefficient should be used.     *
C  If IMETH=2, then the program performs the Gau·-Seidel method. *
C                                                                *
C                                                                *
C  AUXILIARY PARAMETERS:                                         *
C  =====================                                         *
C  WORK   : 2-dim. array  WORK(1:N,1:3)                          *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  A      : 2-dim. array A(1:IA,1:N), the input matrix A is over-*
C           written by: A(I,J)=A(I,J)/A(I,I) for I,J=1, ..., N   *
C  B      : N-vector B(1:N), the right hand side is replaced by  *
C           B(I)=B(I)/A(I,I); I=1,N                              *
C  OMEGA  : - if IMETH = 0, the program returns the adaptively   *
C             computed relaxations parameter.                    *
C           - if IMETH = 1, the optimal relaxation parameter     *
C             is returned as put in externally.                  *
C           - if IMETH = 2, then on output OMEGA = 1.            *
C  X      : N-vector X(1:N) that contains the solution vector    *
C  RES    : N-vector RES(1:N) containing the residuum B - AX;    *
C           the residuum is available even if the desired        *
C           accuracy EPS could not be achieved with the given    *
C           maximum number of iterations.                        *
C  ITNUMB : num,bert of iterations actually performed            *
C  IERR   : error parameter:                                     *
C           = 0, the desired convergence criterium has not been  *
C                met                                             *
C           = 1, the solution X has been found                   *
C           = 2, the desired accuracy has not been achieved after*
C                KMAX iterations                                 *
C           = 3, input data incorrect                            *
C           = 4, system matrix A is numerically singular         *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Required subroutines: GAUSEI, MNORM, CONV, RESID, MACHPD      *
C                                                                *
C*****************************************************************
C                                                                *
C  Author    : Gisela Engeln-MÅllges                             *
C  Date      : 06.09.1992                                        *
C  Source    : FORTRAN 77                                        *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C  Declarations
C
      DOUBLE PRECISION A(1:IA,1:N),B(1:N),X(1:N),WORK(1:N,1:3),
     *                 RES(1:N),EPS,OMEGA,FMACHP,HELP,DIFFN,Q,
     *                 RELERR,SUM,XN
C
C  Checking the inputs EPS, KMAX, IMETH and ISWITC
C
      IF(EPS .LE. 0.0D0 .OR. KMAX .LT. 1 .OR. ISWITC .LT. 0 .OR.
     *   ISWITC .GT. 3 .OR. IMETH .LT. 0 .OR. IMETH .GT. 2) THEN
         IERR=3
         RETURN
      ENDIF
C
C  Initialize the parameters KADAPT and OMEGA depending on the method
C
      IF(IMETH .EQ. 0) THEN
         OMEGA=1.0D0
      ELSE IF(IMETH .EQ. 1) THEN
         KADAPT=KMAX
      ELSE IF(IMETH .EQ. 2) THEN
         KADAPT=KMAX
         OMEGA=1.0D0
      ENDIF
C
C  Compute the machine constant and initialize the relative error bound
C
      FMACHP=1.0D0
   10 FMACHP=0.5D0*FMACHP
      IF(MACHPD(1.0D0+FMACHP) .EQ. 1) GOTO 10
      RELERR=FMACHP*8.0D0
C
C  Initialize
C
      Q=1.0D0
      ITNUMB=0
C
C  Check whether A is singular; if so, set IERR = 4.
C
      DO 20 I=1,N
         SUM=DABS(A(I,1))
         DO 30 K=2,N
            SUM=SUM+DABS(A(I,K))
   30    CONTINUE
         IF(SUM .EQ. 0.0D0) THEN
            IERR=4
            RETURN
         ELSE IF(DABS(A(I,I))/SUM .LT. RELERR) THEN
            IERR=4
            RETURN
         ENDIF
   20 CONTINUE
C
C  Redefine the entries in A and B:   A(I,J) := A(I,J)/A(I,I)
C  and B(I) := B(I)/A(I,I) .
C
      DO 40 I=1,N
         HELP=1.0D0/A(I,I)
         DO 50  J=1,N
            A(I,J)=A(I,J)*HELP
  50     CONTINUE
         B(I)=B(I)*HELP
  40  CONTINUE
C
C  Check for convergence
C
      IF(ISWITC .NE. 0) THEN
         CALL CONV(ISWITC,A,N,IA,IERR)
         IF(IERR .EQ. 0) RETURN
      ENDIF
C
C  The vector RES serves as auxiliary storage for the previous solution
C  vektor. Initially RES contains the staring vector.
C
      DO 60 I=1,N
         RES(I)=X(I)
  60  CONTINUE
C
C  One iteration with the Gau·-Seidel method gives the first iterate X
C
      CALL GAUSEI(A,N,IA,B,OMEGA,X)
C
C  Up the iteration counter
C
      ITNUMB=ITNUMB+1
C
C  Compute the difference of the last two iterates
C
      DO 70  I=1,N
         WORK(I,1)=X(I)-RES(I)
  70  CONTINUE
C
C  Iteration loop for the chosen method
C
      DO 80  K=1,KMAX-1
C
C  Check break-off criterion
C
         CALL MNORM(WORK(1,1),N,DIFFN)
         CALL MNORM(X,N,XN)
         IF(DIFFN .LE. EPS*XN) THEN
            IERR=1
            ITNUMB=K
            CALL RESID(A,N,IA,B,X,RES)
            RETURN
         ENDIF
         IF(K .EQ. KMAX-1) THEN
            ITNUMB=KMAX
            IERR=2
            CALL RESID(A,N,IA,B,X,RES)
            RETURN
         ENDIF
C
C  RES contains the previous iterate
C
         DO 90  I=1,N
            RES(I)=X(I)
  90     CONTINUE
C
C  One iteration step using Gau·-Seidel for a fixed OMEGA
C
         CALL GAUSEI(A,N,IA,B,OMEGA,X)
C
C  Compute the difference of the last two iterates
C
         DO 100  I=1,N
            WORK(I,2)=X(I)-RES(I)
 100     CONTINUE
C
C  If the number of performed iterations K is divisible by KADAPT,
C  then we compute Q in order to adjust the relaxation parameter;
C  Q is an estimate of the spectral radius of the iteration matrix.
C
         IF(MOD(K,KADAPT) .EQ. 0) THEN
            DO 110  I=1,N
               IF(DABS(WORK(I,1)) .LT. FMACHP) THEN
                  WORK(I,3)=1.0D0
               ELSE
                  WORK(I,3)=WORK(I,2)/WORK(I,1)
               ENDIF
 110        CONTINUE
            CALL MNORM(WORK(1,3),N,Q)
C
C  If Q > 1, then the iteration counter is upped by one and
C  the next Gau·-Seidel step is executed; otherwise a new
C  relaxation parameter is calculated.
C
            IF(Q .LE. 1.0D0) THEN
               Q=MAX(Q,OMEGA-1.0D0)
               OMEGA=2.0D0/(1.0D0+DSQRT(1.0D0-((Q+OMEGA-1.0D0)
     *                                        /OMEGA)**2/Q))
            ENDIF
         ENDIF
C
C  The difference vector of the last two iterations is replaced
C  by the one of the previous two iterations for the approximate solution
C
         DO 120  I=1,N
            WORK(I,1)=WORK(I,2)
 120     CONTINUE
  80  CONTINUE
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE GAUSEI(A,N,IA,B,OMEGA,X)
C[IX{GAUSEI}*)
C
C*****************************************************************
C                                                                *
C  This subroutine performs one iteration with the Gau·-Seidel   *
C  method for a given relaxation parameter.                      *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  A      : 2-dim. array A(1:IA, 1:N), that contains the         *
C           modified system matrix A : A(I,J)=A(I,J)/A(I,I) for  *
C           I,J=1, ..., N                                        *
C  N      : order of the system                                  *
C  IA     : leading dimension of A, as specified in the calling  *
C           program                                              *
C  B      : N-vector B(1:N) with the modified right hand side:   *
C           B(I)=B(I)/A(I,I); I=1, ..., N                        *
C  OMEGA  : relaxation parameter                                 *
C  X      : N-vector X(1:N) containing the starting vector for   *
C           the iteration                                        *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  X      : N-vector X(1:N) containing the next iteration vector *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Required subroutines: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : Gisela Engeln-MÅllges                              *
C  Date     : 06.09.1992                                         *
C  Source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      DOUBLE PRECISION A(1:IA,1:N),B(1:N),X(1:N),OMEGA,S
C
      DO 10  I=1,N
         S=B(I)
         DO 20  J=1,N
            S=S-A(I,J)*X(J)
  20     CONTINUE
         X(I)=X(I)+OMEGA*S
  10  CONTINUE
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE MNORM(X,N,XNORM)
C[IX{MNORM}*)
C
C*****************************************************************
C                                                                *
C  This subroutine calculates the maximum norm XNORM of an       *
C  N-vector X.                                                   *
C                                                                *
C----------------------------------------------------------------*
C[BE*)
C                                                                *
C  Required subroutines: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : Gisela Engeln-MÅllges                              *
C  Date     : 06.09.1992                                         *
C  Source   : FORTRAN 77                                         *
C                                                                *
C*****************************************************************
C
      DOUBLE PRECISION X(1:N),XNORM
C
      XNORM=DABS(X(1))
      DO 10  I=2,N
         XNORM=DMAX1(XNORM,DABS(X(I)))
  10  CONTINUE
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE CONV(ISWITC,A,N,IA,IERR)
C[IX{CONV}*)
C
C*****************************************************************
C                                                                *
C  This subroutine helps check convergence.                      *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  ISWITC : Parameter that determines the convergence criterion  *
C           to be checked:                                       *
C           = 0, none                                            *
C           = 1, row sum criterion                               *
C           = 2, column sum criterion                            *
C           = 3, criterion of Schmidt and v. Mises               *
C  A      : 2-dim. array A(1:IA, 1:N), containing the matrix for *
C           which we want to check convergence of the iterates   *
C           from the various SOR algorithms                      *
C  N      : order of the matrix A                                *
C  IA     : leading dimension of A, as prescribed in the calling *
C           program                                              *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  IERR   : error parameter:                                     *
C           = 0, the desired convergence criterion has not been  *
C                met                                             *
C           = 1, the desired criterion is satified               *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Required subroutines: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : Gisela Engeln-MÅllges                              *
C  Date     : 06.09.1992                                         *
C  Source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      DOUBLE PRECISION A(1:IA,1:N),SUM
C
C  Row sum criterion
C
      IF(ISWITC .EQ. 1) THEN
         DO 10  I=1,N
            SUM=-1.0D0
            DO 20  J=1,N
               SUM=SUM+DABS(A(I,J))
  20        CONTINUE
            IF(SUM .LT. 1.0D0) THEN
               IERR=1
            ELSE
               IERR=0
               RETURN
            ENDIF
  10     CONTINUE
C
C  Column sum criterion
C
      ELSE IF(ISWITC .EQ. 2) THEN
         DO 30  J=1,N
            SUM=-1.0D0
            DO 40  I=1,N
               SUM=SUM+DABS(A(I,J))
  40        CONTINUE
            IF(SUM .LT. 1.0D0) THEN
               IERR=1
            ELSE
               IERR=0
               RETURN
            ENDIF
  30     CONTINUE
C
C  Criterion of Schmidt and v. Mises
C
      ELSE IF(ISWITC .EQ. 3) THEN
         SUM=-N
         DO 50  I=1,N
            DO 60  J=1,N
               SUM=SUM+A(I,J)*A(I,J)
  60        CONTINUE
  50     CONTINUE
         SUM=DSQRT(SUM)
         IF(SUM .LT. 1.0D0) THEN
            IERR=1
         ELSE
            IERR=0
            RETURN
         ENDIF
      ENDIF
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE RESID(A,N,IA,B,X,RES)
C[IX{RESID}*)
C
C*****************************************************************
C                                                                *
C  This subroutine computes the residuum  RES = B - AX, where    *
C  both A and B are given in modified form.                      *
C                                                                *
C----------------------------------------------------------------*
C[BE*)
C                                                                *
C  Required subroutines: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : Gisela Engeln-MÅllges                              *
C  Date     : 09.06.1992                                         *
C  Source   : FORTRAN 77                                         *
C                                                                *
C*****************************************************************
C
      DOUBLE PRECISION A(1:IA,1:N), B(1:N), X(1:N), RES(1:N),DSUM
C
      DO 10  I=1,N
         DSUM=B(I)
         DO 20 J=1,N
            DSUM=DSUM-A(I,J)*X(J)
  20     CONTINUE
         RES(I)=DSUM
  10  CONTINUE
      RETURN
      END
