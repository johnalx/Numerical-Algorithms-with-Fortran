C[BA*)
C[LE*)
      SUBROUTINE WGKNOT (N,WORK,IERR)
C[IX{WGKNOT}*)
C
C*****************************************************************
C                                                                *
C  This subroutine determines the weights and the nodes of the   *
C  CLENSHAW-CURTIS quadrature formula of local error order N+3   *
C  for the reference interval [-1,1].                            *
C[BE*)
C                                                                *
C  INPUT PARAMETER:                                              *
C  ================                                              *
C  N    : N+1 denotes the number of nodes and weights,           *
C         N >= 2, and N must be even.                            *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  WORK : 2-dimensional array WORK(0:N,2) containing the weights *
C         and nodes of the quadrature formula for the interval   *
C         [-1,1]; the weights appear in the first column, with   *
C         nodes in the second column of WORK.                    *
C  IERR : error parameter.                                       *
C           IERR = 0: everything o.k.                            *
C           IERR = 1: condition for N is not met                 *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required : none                                   *
C                                                                *
C*****************************************************************
C                                                                *
C  author    : Gisela Engeln-Muellges                            *
C  date      : 05.10.1989                                        *
C  source    : FORTRAN 77                                        *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C  declarations
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  WORK(0:N,2)
C
C  initializing
C
      PI = 4.0D0 * DATAN(1.0D0)
C
C  checking the validity of N
C
      IERR = 1
      IF (N .LT. 2 .OR. MOD(N,2) .NE. 0) RETURN
C
C  determine the weights and nodes
C
      IERR      = 0
      DUMMY     = N * N - 1
      WORK(0,1) = 1.0D0 / DUMMY
      WORK(N,1) = WORK(0,1)
      DUMMY1    = PI / N
      WORK(0,2) = 1.0D0
      WORK(N,2) = -1.0D0
      DO 10 K = 1, N-1
         DUMMY2 = 2.0D0 * (DUMMY- (-1) ** K) / (N * DUMMY)
         DUMMY3 = 0.0D0
         DO 20 L = 1, N/2-1
            DUMMY3 = DUMMY3 + DCOS(2.0D0*L*K*DUMMY1) / (4*L*L-1)
 20      CONTINUE
         WORK(K,1) = DUMMY2 - 4.0D0 / N * DUMMY3
         WORK(K,2) = DCOS(K * DUMMY1)
 10   CONTINUE
      RETURN
      END
