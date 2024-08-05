C[BA*)
C[LE*)
      SUBROUTINE CLENSH (FCT,N,Z,M,WORK,QCCN,IERR)
C[IX{CLENSH}*)
C
C*****************************************************************
C                                                                *
C  This subroutine computes the definite integral I(FCT;A,B) of  *
C  FCT over the interval [A,B] for the partition                 *
C     Z : A = Z(0) < Z(1) < .... < Z(M) = B                      *
C  using a summed CLENSHAW-CURTIS formula.                       *
C[BE*)
C  If M = 1, the integral is evaluated directly over the whole   *
C  interval [A,B].                                               *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  FCT    : function to be integrated. Its format is as follows  *
C                  DOUBLE PRECISION FUNCTION FCT(X).             *
C           In the calling program it has to be defined as       *
C           EXTERNAL (or as INTRINSIC, if a FTN5-standard-       *
C           function is used).                                   *
C  N      : N+1 denotes the number of nodes used in the reference*
C           interval, N+2 is the global error order, N >= 2, even*
C  Z      : vector Z(0:M) containing the nodes Z(I), I=0,...,M of*
C           the partition.                                       *
C  M      : M describes the number of sub-intervals [Z(I),Z(I+1)]*
C           I=0, 1, ..., M-1 with M >= 1.                        *
C  WORK   : 2-dimensional array  WORK(0:N,2) containing the N+1  *
C           weights of the CLENSHAW-CURTIS formula in its 1st    *
C           column and the N+1 nodes in the 2nd column.          *
C           In the calling program, the  SUBROUTINE WGKNOT must  *
C           be used before calling this subroutine in order that *
C           the weights and nodes become available.              *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  QCCN   : approximate value for the integral.                  *
C  IERR   : error parameter.                                     *
C           IERR = 0: everything o.k.                            *
C           IERR = 1: condition for N is not met                 *
C           IERR = 2: condition for M is not met                 *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Gisela Engeln-Muellges                             *
C  date     : 04.19.1988                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C  declarations
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  Z(0:M),WORK(0:N,2)
C
C  testing the input for correctness
C
      IERR = 1
      IF (N .LT. 2 .OR. MOD(N,2) .NE. 0) RETURN
      IERR = 2
      IF (M .LT. 1) RETURN
C
C  determine the approximate value for the integral
C
      IERR = 0
      QCCN  = 0.0D0
      DO 10 I = 0, M-1
         DUMMY1 = 0.5D0 * (Z(I+1) - Z(I))
         DUMMY2 = 0.5D0 * (Z(I+1) + Z(I))
         SUM   = 0.0D0
         DO 20 K = 0, N
            SUM = SUM + WORK(K,1) * FCT (DUMMY1*WORK(K,2)+DUMMY2)
 20      CONTINUE
         QCCN = QCCN + SUM*DUMMY1
C
 10   CONTINUE
      RETURN
      END
