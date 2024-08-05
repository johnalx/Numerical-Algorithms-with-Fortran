C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{}
C[  {Condition Estimate according to Forsythe/Moler}
C[  {Condition Estimate according to Forsythe/Moler}*)
C[LE*)
      SUBROUTINE CONDAE(N,AO,A,IA,IPIVOT,Y,X,Z,R,CONDA)
C[IX{CONDAE}*)
C
C*****************************************************************
C                                                                *
C  The SUBROUTINE CONDAE finds an estimate of the matrix         *
C  condition number                                              *
C          COND(A) =  Max norm of A * Max norm of A inverse,     *
C  following the condition estimate of  FORSYTHE and MOLER.      *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  N        : order of the matrices A and AO.                    *
C  AO       : DOUBLE PRECISION array AO(1:IA,1:N) containing the *
C             orginal matrix A = A(ORG).                         *
C  A        : DOUBLE PRECISION array A(1:IA,1:N), containing the *
C             output matrix of the SUBROUTINE GAUSSP, i.e., the  *
C             factors L and R with  P*A(ORG)=L*R for a           *
C             permutation matrix P.                              *
C  IA       : leading dimension of A and AO, as stipulated by    *
C             the calling program.                               *
C  IPIVOT   : INTEGER N-vector IPIVOT(1:N), output of the        *
C             SUBROUTINE GAUSS, i.e., information on P.          *
C  Y        : DOUBLE PRECISION vector Y(1:N) containing the right*
C             hand side of the linear system  A(ORG)*X=Y.        *
C  X        : DOUBLE PRECISION vector X(1:N) with the solution   *
C             of the linear system from  GAUSS.                  *
C                                                                *
C                                                                *
C  AUXILIARY PARAMETERS:                                         *
C  =====================                                         *
C  Z        : ]  DOUBLE PRECISION vectors ..(1:N).               *
C  R        : ]                                                  *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  CONDA    : Estimate for  COND(A).                             *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Required subroutines:  GAUSSS, MACHPD                         *
C                                                                *
C*****************************************************************
C                                                                *
C  Author    : Gisela Engeln-MÅllges                             *
C  Date      : 06.07.90                                          *
C  Source    : FORTRAN 77                                        *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C Declarations
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION AO(1:IA,1:N), A(1:IA,1:N), Y(1:N), X(1:N),
     +                 Z(1:N), R(1:N)
      INTEGER IPIVOT(1:N)
C
C Compute the max norm of the solution  X  of  A X = Y
C using the Gau· algorithm
C
      XMAX=DABS(X(1))
      DO 10 I=2,N
         XMAX=DMAX1(XMAX,DABS(X(I)))
10    CONTINUE
C
C Find the machine constant
C
      FMACHP=1.0D0
20    FMACHP=0.5D0*FMACHP
      IF(MACHPD(1.0D0+FMACHP) .EQ. 1) GOTO 20
      EPS=2.0D0*FMACHP
C
C Find the residuum  Y - A * X; form A * X  in double precision,
C then round the residuum to single precision
C
      DO 30 I=1,N
         R(I)=Y(I)
         DO 40 K=1,N
            R(I)=R(I)-AO(I,K)*X(K)
40       CONTINUE
30    CONTINUE
C
C Calculate the first correction Z to the solution  and its max norm
C
      CALL GAUSSS(N,A,IA,IPIVOT,R,Z)
      ZMAX=DABS(Z(1))
      DO 50 I=2,N
         ZMAX=DMAX1(ZMAX,DABS(Z(I)))
50    CONTINUE
C
C Estimate the condition number  COND(A)=NORM(A)*NORM(INV(A))
C from the sizes of Z, X and EPS
C
      CONDA=ZMAX/(XMAX*EPS)
C
      RETURN
      END
