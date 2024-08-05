C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 12.4.1}{B-Spline Curves}{B-Spline Curves}*)
C[LE*)
      SUBROUTINE CURVP (NP,N,M,K,KV,DP,X,XP,D,E,IERR,IX)
C[IX{CURVP}*)
C
C*****************************************************************
C                                                                *
C  The subroutine CURVP computes at most NP+1 points on an open  *
C  uniform B-Spline curve of order K.                            *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  NP     : number of desired points on the curve:               *
C           maximally NP+1 points and at least 2*(N-K+2)+1 points*
C           are computed, i.e., NP=MAX(NP,2*(N-K+2))             *
C  N      : N+1 is the number of DE BOOR points                  *
C  M      : dimension of the space for the DE BOOR points        *
C           ( M >= 2)                                            *
C  K      : order of the B-Spline curve (3 <= K <= N+1)          *
C  KV     : INTEGER vector KV(1:N+K-1) for the node vectors      *
C  DP     : DOUBLE PRECISION array DP(0:N,1:M) containing the    *
C           N+1 DE BOOR points, N >= 2                           *
C                                                                *
C                                                                *
C  AUXILIARY VARIABLES:                                          *
C  ====================                                          *
C  X      : DOUBLE PRECISION vector X(1:M)                       *
C  D,E    : DOUBLE PRECISION 2 dimensional arrays D(1:K,1:M),    *
C           E(1:K,1:M)                                           *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  XP     : DOUBLE PRECISION array XP(1:NP+1,1:M) containing the *
C           computed points on the curve (each row contains the M*
C           coordinates of a point; their maximal number is NP+1)*
C  IERR   : error parameter                                      *
C           IERR=0, everything is o.k.                           *
C           IERR=1, input conditions violated                    *
C  IX     : Number of points stored in XP                        *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  required subroutines: DEBOOR, KNOTVO                          *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : Gisela Engeln-MÅllges                              *
C  Date     : 11.30.1991                                         *
C  Source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C     Declarations
C
      DOUBLE PRECISION DP(0:N,12:M),X(1:M),XP(1:NP+1,1:M),D(1:K,1:M),
     +                 E(1:K,1:M), H, T
      INTEGER KV(1:N+K-1),IX
C
C     Checking the input data
C
      IF (N .LT. 2 .OR. K .LT. 3 .OR. K .GT. N+1) THEN
         IERR=1
         RETURN
      ENDIF
C
C     Call SUBROUTINE KNOTVO in order to compute the nodes
C
      CALL KNOTVO(N,K,KV)
C
C     Compute the number of points between two adjacent nodes
C     of the node vector and the step size H
C
      NI=INT(NP/(N-K+2))
      NI=MAX(NI,2)
      H=1.0D0/DBLE(NI)
C
C     Call SUBROUTINE DEBOOR in order to compute the coordinates
C     of the first point of the curve; store in the first row of XP
C
      CALL DEBOOR(N,M,DP,K,KV,DBLE(K-1),K-1,D,E,X,IERR)
C
      DO 10 I=1,M
         XP(1,I)=X(I)
   10 CONTINUE
C
C     Compute the subsequent points of the curve
C
      IX=1
      DO 20 IR=K-1,N
         T=DBLE(IR)
         DO 30 L=1,NI
            T=T+H
            CALL DEBOOR(N,M,DP,K,KV,T,IR,D,E,X,IERR)
            IX=IX+1
C
C           Store the coordinates of the computed point X on the curve
C           in row IX of XP, 2 <= IX <= NP+1
C
            DO 40 I=1,M
               XP(IX,I)=X(I)
   40       CONTINUE
   30    CONTINUE
   20 CONTINUE
      RETURN
      END
