C[BA*)
C[LE*)
      SUBROUTINE KNOTVG (N,K,KV)
C[IX{KNOTVG}*)
C[LE*)
C********************************************************************
C                                                                   *
C   This program computes the node vector KV of a closed uniform B  *
C   spline.                                                         *
C[BE*)
C                                                                   *
C   INPUT PARAMETERS                                                *
C   ================                                                *
C   N       : N+1 is the number of DE BOOR points, N >= 2           *
C   K       : degree of the B spline curve, 2 <= K <= N+1           *
C                                                                   *
C   OUTPUT PARAMETERS                                               *
C   =================                                               *
C   KV      : INTEGER vector KV(1:N+K-1), with the node vector of   *
C             order K                                               *
C                                                                   *
C-------------------------------------------------------------------*
C                                                                   *
C   Required subroutines: none                                      *
C                                                                   *
C********************************************************************
C                                                                   *
C   Authors     : Reinhold Wodicka, Bjoern Terwege                  *
C   Date        : 05.07.1995                                        *
C   Source code : FORTRAN 77                                        *
C                                                                   *
C[BA*)
C********************************************************************
C[BE*)
C
      INTEGER KV(1:N+2*(K-1))
c
c     Create the node vector KV of a closed uniform B spline
c
      DO 10 J=1,N+2*(K-1)
          KV(J)=j
10    CONTINUE
      RETURN
      END
