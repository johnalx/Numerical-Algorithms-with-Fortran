![          {Clenshaw--Curtis Quadrature Formulas}*)                    
      SUBROUTINE CCFERR (FCT, N, Z, Z2, M, WORK, QCCNZ, QCCNZ2, ERREST, &
      QCCNST, IERR)                                                     
!                                                                       
!*****************************************************************      
!                                                                *      
!  This subroutine computes an approximate value for the         *      
!  definite integral I(FCT; Z(0),Z(M)) of the function FCT over  *      
!  the interval [Z(0),Z(M)] using the CLENSHAW-CURTIS formula of *      
!  global error order N+2, where N >= 2 is even, for a partition *      
!  Z : Z(0) < Z(1) < .... < Z(M), and for the partition Z/2  ob- *      
!  tained by halving the local stepsizes.                        *      
!  Both approximate values help give an estimate for the global  *      
!  procedural error of the method and are used to find an im-    *      
!  proved approximate value for the integral.                    *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  FCT    : function to be integrated. It is formatted as        *      
!                  DOUBLE PRECISION FUNCTION FCT(X)              *      
!           and has to be defined as EXTERNAL in the calling     *      
!           program (or as INTRINSIC, if a FTN5 standard         *      
!           function is used).                                   *      
!  N      : N+1 is the number of nodes used in the reference     *      
!           interval, N+2 is the global error order. N >= 2, even*      
!  Z      : vector Z(0:M) containing the nodes Z(I), I=0, ..., M *      
!           of the partition                                     *      
!  M      : M denotes the number of sub-intervals [Z(I),Z(I+1)], *      
!           I=0, 1, ..., M-1, with M >= 1.                       *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  WORK   : 2-dimensional array WORK(0:N,2) containing the       *      
!           weights and nodes of the quadrature formula          *      
!           for the reference interval [-1,1]; the weights are   *      
!           in the first column, with the nodes in the second    *      
!           column.                                              *      
!  QCCNZ  : approximate value for I(FCT; Z(0),Z(M)) for the      *      
!           partition  Z.                                        *      
!  QCCNZ2 : approximate value for I(FCT; Z(0),Z(M)) for the      *      
!           partition  Z/2.                                      *      
!  ERREST : estimate for the global procedural error of the      *      
!           approximate value QCCNZ2 of the integral using the   *      
!           partition Z/2.                                       *      
!  QCCNST : improved approximate value for I(FCT; Z(0),Z(M)).    *      
!  IERR   : error parameter.                                     *      
!           IERR = 0: everything o.k.                            *      
!           IERR = 1: condition for N not met                    *      
!           IERR = 2: condition for M not met                    *      
!                                                                *      
!  AUXILIARY PARAMETERS:                                         *      
!  =====================                                         *      
!  Z2     : vector Z2(0:2*M), containing the partition Z/2.      *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: WGKNOT, CLENSH                          *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Gisela Engeln-Muellges                             *      
!  date     : 04.19.1988                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
!  declarations                                                         
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION Z (0:M), Z2 (0:2 * M), WORK (0:N, 2) 
      EXTERNAL FCT 
!                                                                       
!  test that M >= 1                                                     
!                                                                       
      IERR = 2 
      IF (M.LT.1) RETURN 
!                                                                       
!  determine the N+1 nodes and N+1 weights of the CLENSHAW-             
!  CURTIS formula for the reference interval [-1,1]                     
!                                                                       
      CALL WGKNOT (N, WORK, IERR) 
      IF (IERR.NE.0) RETURN 
!                                                                       
!  determine the partition Z/2                                          
!                                                                       
      DO 10 I = 0, M - 1 
         Z2 (2 * I) = Z (I) 
         Z2 (2 * I + 1) = 0.5D0 * (Z (I) + Z (I + 1) ) 
   10 END DO 
      Z2 (2 * M) = Z (M) 
!                                                                       
!  determine the approximate value QCCNZ for the                        
!  integral of FCT over [Z(0),Z(M)] using the CLENSHAW-CURTIS           
!  formula for the partition Z                                          
!                                                                       
      CALL CLENSH (FCT, N, Z, M, WORK, QCCNZ, IERR) 
      IF (IERR.NE.0) RETURN 
!                                                                       
!  determine the approximate value QCCNZ2 for the                       
!  integral of FCT over [Z(0),Z(M)] using the CLENSHAW-CURTIS           
!  formula for the partition Z/2                                        
!                                                                       
      CALL CLENSH (FCT, N, Z2, 2 * M, WORK, QCCNZ2, IERR) 
      IF (IERR.NE.0) RETURN 
!                                                                       
!  determine an estimate for the global                                 
!  procedural error                                                     
!                                                                       
      ERREST = (QCCNZ2 - QCCNZ) / (2.0D0** (N + 2) - 1.0D0) 
!                                                                       
!  determine an improved approximatation QCCNST for                     
!  the integral                                                         
!                                                                       
      QCCNST = QCCNZ2 + ERREST 
      RETURN 
      END SUBROUTINE CCFERR                         
