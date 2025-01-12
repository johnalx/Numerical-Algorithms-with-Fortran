C[BA*)
C[LE*)
      SUBROUTINE CURVPG (NP,N,M,K,KV,DP,X,XP,D,E,IERR)
C[IX{CURVPG}*)

C********************************************************************
C                                                                   *
C   The subroutine CURVPG computes NP points of a closed uniform    *
C   b spline curve of order K.                                      *
C[BE*)
C                                                                   *
C   INPUT PARAMETERS                                                *
C   ================                                                *
C   NP      : desired number of points on the curve                 *
C   N       : number of DE BOOR pointsis N - 1, N >= 2              * 
C   M       : Dimension of DE BOOR points (M >= 2)                  *
C   K       : Order of the b spline, 3 <= K <= N+1                  *
C   KV      : INTEGER vector KV(1:N+2*(K-1)) with the nodes         *
C   DP      : DOUBLE PRECISION array DP(0:N+K-1,1:M); only N+1      *
C             points used on input                                  *
C                                                                   *
C   AUX ARRAYS                                                      *
C   ===========                                                     *
C   X       : DOUBLE PRECISION vector X(1:M)                        *
C   D,E     : DOUBLE PRECISION arrays ..(1:K,1:M)                   *
C                                                                   *
C   OUTPUT PARAMETERS                                               *
C   =================                                               *
C   XP      : DOUBLE PRECISION array XP(1:NP,1:M) with the computed *
C             points; each row of XP contains the M coordinates of  *
C             one point                                             *
C   IERR    : error parameter                                       *
C             IERR=0, all o k                                       *
C             IERR=1, Input incorrect                               *
C                                                                   *
C-------------------------------------------------------------------*
C                                                                   *
C   Required subroutines: KNOTVG, DEBOOR                            *
C                                                                   *
C********************************************************************
C                                                                   *
C   Author      : Reinhold Wodicka, Bj�rn Terwege                    *
C   Date        : 06.06.1995                                        *
C   Source code : FORTRAN 77                                        *
C                                                                   *
C[BA*)
C********************************************************************
C[BE*)
C
C
C  Declarations
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION KV(1:N+2*(K-1)),DP(0:N+K-1,1:M),X(1:M),XP(1:NP,1:M)
     f ,D(1:K,1:M),E(1:K,1:M)
      IERR=0
C
C  Check input
C
      IF(N.LT.2.OR.K.LT.3.OR.K.GT.N+1) THEN
          IERR=1
          RETURN
      ENDIF
C
C  Append the extra DE BOOR points to DP
C
      J=0
      DO 10 I=N+1,N+K-1
         DO 20 J=1,M
           DP(I,J)=DP(I-N-1,J)
20       CONTINUE
10    CONTINUE
C
C  Call SUBROUTINE KNOTVG to compute the nodes
C
      CALL KNOTVG(N,K,KV)
C
C  Determine step size
C
      DT=DBLE(N+1)/DBLE(NP-1)
      T=DBLE(K-1)
      IR=K-1
C
C  compute points of the curve
C
      DO 30 I=1,NP
         CALL DEBOOR(N+K-1,M,DP,K,KV,T,IR,D,E,X,IERR)
C
C  Save the coordinates of the newly computed point X in row I of XP
C
         DO 40 L=1,M
            XP(I,L)=X(L)
40       CONTINUE
         T=DMIN1(T+DT,DBLE(N+K))
25       IF(T.GT.DBLE(IR+1)) THEN
            IR=IR+1
            GOTO 25
         ENDIF
30    CONTINUE
      RETURN
      END
