C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 15.9}
C[  {Clenshaw-Curtis Quadrature Formulas}
C[  {Clenshaw-Curtis Quadrature Formulas}*)
C[LE*)
      SUBROUTINE CCFERR (FCT,N,Z,Z2,M,WORK,QCCNZ,QCCNZ2,ERREST,
C[IX{CCFERR}*)
     +                   QCCNST,IERR)
C
C*****************************************************************
C                                                                *
C  This subroutine computes an approximate value for the         *
C  definite integral I(FCT; Z(0),Z(M)) of the function FCT over  *
C  the interval [Z(0),Z(M)] using the CLENSHAW-CURTIS formula    *
C[BE*)
C  of global error order N+2, where N >= 2 is even, for a        *
C  partition Z : Z(0) < Z(1) < .... < Z(M), and for the          *
C  partition Z/2  obtained by halving the local stepsizes.       *
C  Both approximate values help give an estimate for the global  *
C  procedural error of the method and are used to find an im-    *
C  proved approximate value for the integral.                    *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  FCT    : function to be integrated. It is formatted as        *
C                  DOUBLE PRECISION FUNCTION FCT(X)              *
C           and has to be defined as EXTERNAL in the calling     *
C           program (or as INTRINSIC, if a FTN5 standard         *
C           function is used).                                   *
C  N      : N+1 is the number of nodes used in the reference     *
C           interval, N+2 is the global error order. N >= 2, even*
C  Z      : vector Z(0:M) containing the nodes Z(I), I=0, ..., M *
C           of the partition                                     *
C  M      : M denotes the number of sub-intervals [Z(I),Z(I+1)], *
C           I=0, 1, ..., M-1, with M >= 1.                       *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  WORK   : 2-dimensional array WORK(0:N,2) containing the       *
C           weights and nodes of the quadrature formula          *
C           for the reference interval [-1,1]; the weights are   *
C           in the first column, with the nodes in the second    *
C           column.                                              *
C  QCCNZ  : approximate value for I(FCT; Z(0),Z(M)) for the      *
C           partition  Z.                                        *
C  QCCNZ2 : approximate value for I(FCT; Z(0),Z(M)) for the      *
C           partition  Z/2.                                      *
C  ERREST : estimate for the global procedural error of the      *
C           approximate value QCCNZ2 of the integral using the   *
C           partition Z/2.                                       *
C  QCCNST : improved approximate value for I(FCT; Z(0),Z(M)).    *
C  IERR   : error parameter.                                     *
C           IERR = 0: everything o.k.                            *
C           IERR = 1: condition for N not met                    *
C           IERR = 2: condition for M not met                    *
C                                                                *
C  AUXILIARY PARAMETERS:                                         *
C  =====================                                         *
C  Z2     : vector Z2(0:2*M), containing the partition Z/2.      *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: WGKNOT, CLENSH                          *
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
      DIMENSION Z(0:M),Z2(0:2*M),WORK(0:N,2)
      EXTERNAL FCT
C
C  test that M >= 1
C
      IERR = 2
      IF (M .LT. 1) RETURN
C
C  determine the N+1 nodes and N+1 weights of the CLENSHAW-
C  CURTIS formula for the reference interval [-1,1]
C
      CALL WGKNOT (N,WORK,IERR)
      IF (IERR .NE. 0) RETURN
C
C  determine the partition Z/2
C
      DO 10 I=0,M-1
         Z2(2*I) = Z(I)
         Z2(2*I+1) = 0.5D0*(Z(I)+Z(I+1))
   10 CONTINUE
      Z2(2*M) = Z(M)
C
C  determine the approximate value QCCNZ for the
C  integral of FCT over [Z(0),Z(M)] using the CLENSHAW-CURTIS
C  formula for the partition Z
C
      CALL CLENSH (FCT,N,Z,M,WORK,QCCNZ,IERR)
      IF (IERR .NE. 0) RETURN
C
C  determine the approximate value QCCNZ2 for the
C  integral of FCT over [Z(0),Z(M)] using the CLENSHAW-CURTIS
C  formula for the partition Z/2
C
      CALL CLENSH (FCT,N,Z2,2*M,WORK,QCCNZ2,IERR)
      IF (IERR .NE. 0) RETURN
C
C  determine an estimate for the global
C  procedural error
C
      ERREST = (QCCNZ2-QCCNZ)/(2.0D0**(N+2)-1.0D0)
C
C  determine an improved approximatation QCCNST for
C  the integral
C
      QCCNST = QCCNZ2 + ERREST
      RETURN
      END
