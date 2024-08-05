C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 4.15.4}{Iterative Refinement}{Iterative Refinement}*)
C[LE*)
      SUBROUTINE POSTIT(N,A0,A,LDA,IPIVOT,Y,X,EPS,MAXIT,NUMIT,
C[IX{POSTIT}*)
     +                  IERR,Z,R,RS)
C
C*****************************************************************
C                                                                *
C  POSTIT performs iterative refinement after GAUSSP and GAUSSS  *
C  have been executed.                                           *
C[BE*)
C  The iteration is stopped if a set maximum                     *
C  MAXIT of iteration steps has been performed or if the relative*
C  improvement satisfies  MAX(ABS(Z(I)))/MAX(ABS(X(I))) < EPS    *
C  for I=1, ..., N.                                              *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  N        : order of the matrices A and A0.                    *
C  A0       : 2-dimensional array A0(1:LDA,1:N); the matrix      *
C             A(ORG).                                            *
C  A        : 2-dimensional array A(1:LDA,1:N) containing the    *
C             factors  L  and  R  with  P * A(ORG) = L * R.      *
C             P = permutation matrix. A is an output of          *
C             SUBROUTINE GAUSSP.                                 *
C  LDA      : leading dimension of A and A0 as defined in the    *
C             calling program.                                   *
C  IPIVOT   : N-vector IPIVOT(1:N); it indicates the row permuta-*
C             tions in  P * A(ORG)  compared with A(ORG). It is  *
C             an output vector of SUBROUTINE GAUSSP.             *
C  Y        : N-vector Y(1:N); the right side of the system of   *
C             equations.                                         *
C  X        : N-vector X(1:N); the solution of the system of     *
C             equations; X serves as starting vector for the     *
C             iterative refinement.                              *
C  EPS      : error bound for the relative improvement;          *
C             if EPS < 4 * machine constant, the program inter-  *
C             ally sets EPS = 4 * machine constant.              *
C  MAXIT    : maximum number of iterations.                      *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  X        : N-vector X(1:N); the solution of the linear system.*
C  EPS      : error bound actually used.                         *
C  NUMIT    : number of iteration steps executed.                *
C  IERR     : error parameter.                                   *
C             = 0 : program stopped since the relative           *
C                   improvement is < EPS.                        *
C             = 1 : the set accuracy was not achieved after      *
C                   MAXIT iterations.                            *
C             = 2 : same as IERR=1, except that the components of*
C                   the correction vector have begone to increase*
C                   indicating divergence and ill-conditioning of*
C                   A0.                                          *
C                                                                *
C                                                                *
C  AUXILIARY VECTORS:                                            *
C  ==================                                            *
C  Z        : N-vector Z(1:N).                                   *
C  R        : N-vector R(1:N) in DOUBLE PRECISION.               *
C  RS       : N-vector RS(1:N).                                  *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: GAUSSS, MACHPD                          *
C                                                                *
C*****************************************************************
C                                                                *
C  authors  : Gisela Engeln-Muellges, Guido Dubois               *
C  date     : 04.25.88                                           *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A0(1:LDA,1:N),A(1:LDA,1:N),Y(1:N),X(1:N),
     +                 Z(1:N),RS(1:N)
      DOUBLE PRECISION R(1:N)
      INTEGER IPIVOT(1:N)
C
C  local storage of the minimally acceptable error bound EPSMIN
C  in case that the subroutine is called more than once.
C
      SAVE EPSMIN,IFLAG
      DATA IFLAG /0/
C
C  calculating the machine constant and initializing the minimally
C  acceptable error bound EPSMIN for the relative improvement.
C
      IF(IFLAG .EQ. 0) THEN
         IFLAG=1
         FMACHP=1.0D0
   10    FMACHP=0.5D0*FMACHP
         IF(MACHPD(1.0D0+FMACHP) .EQ. 1) GOTO 10
         EPSMIN=8.0D0*FMACHP
      END IF
      IF(EPS .LT. EPSMIN) EPS=EPSMIN
      NUMIT=0
      ZMA=1.0D20
      IERR=0
   50 NUMIT=NUMIT+1
C
C  calculation of the residual vector.
C
      DO 20 I=1,N
         R(I)=Y(I)
         DO 30 K=1,N
            R(I)=R(I)-A0(I,K)*X(K)
   30    CONTINUE
         RS(I)=R(I)
   20 CONTINUE
C
C  calculating the correction vector and improving
C  the approximate solution.
C
      CALL GAUSSS(N,A,LDA,IPIVOT,RS,Z)
      XM=0.0D0
      ZM=0.0D0
      DO 40 I=1,N
         XM=DMAX1(XM,DABS(X(I)))
         ZM=DMAX1(ZM,DABS(Z(I)))
         X(I)=X(I)+Z(I)
   40 CONTINUE
      IF(ZM/XM .LT. EPS) RETURN
      IF(NUMIT .GE. MAXIT) THEN
         IERR=1
         IF(ZM .GT. ZMA) IERR=2
         RETURN
      END IF
      ZMA=ZM
      GOTO 50
      END
