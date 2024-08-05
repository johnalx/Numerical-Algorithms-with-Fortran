C[BA*)
C[LE*)
      SUBROUTINE KNOTVO (N,K,KV)
C[IX{KNOTVO}*)
C
C*****************************************************************
C                                                                *
C                                                                *
C   This program computes the node vector KV of an open uniform  *
C   B spline.                                                    *
C[BE*)
C                                                                *
C   INPUT PARAMETERS                                             *
C   ================                                             *
C   N       : N+1 is the number of DE BOOR points, N >= 2        *
C   K       : degree of the B spline curve, 3 <= K <= N+1        *
C                                                                *
C   OUTPUT PARAMETERS                                            *
C   =================                                            *
C   KV      : INTEGER vector KV(1:N+K+1), with the node vector   *
C             of order K                                         *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C   Required subroutines: none                                   *
C                                                                *
C*****************************************************************
C                                                                *
C  Author      : Gisela Engeln-Muellge                           *
C  Date        : 11.30.91                                        *
C  Source code : FORTRAN 77                                      *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      INTEGER KV(1:N+1+K)
C
C      Create the node vector KV of an open uniform B spline
C
      DO 10 J=1,K-1
         KV(J)=K-1
  10  CONTINUE
      DO 20 J=K,N
         KV(J)=J
  20  CONTINUE
      DO 30 J=N+1,N+K-1
         KV(J)=N+1
  30  CONTINUE
      RETURN
      END
