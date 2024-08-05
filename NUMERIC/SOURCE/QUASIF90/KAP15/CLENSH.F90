      SUBROUTINE CLENSH (FCT, N, Z, M, WORK, QCCN, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This subroutine computes the definite integral I(FCT;A,B) of  *      
!  FCT over the interval [A,B] for the partition                 *      
!     Z : A = Z(0) < Z(1) < .... < Z(M) = B                      *      
!  using a summed CLENSHAW-CURTIS formula.                       *      
!  If M = 1, the integral is evaluated directly over the whole   *      
!  interval [A,B].                                               *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  FCT    : function to be integrated. Its format is as follows  *      
!                  DOUBLE PRECISION FUNCTION FCT(X).             *      
!           In the calling program it has to be defined as       *      
!           EXTERNAL (or as INTRINSIC, if a FTN5-standard-       *      
!           function is used).                                   *      
!  N      : N+1 denotes the number of nodes used in the reference*      
!           interval, N+2 is the global error order, N >= 2, even*      
!  Z      : vector Z(0:M) containing the nodes Z(I), I=0,...,M of*      
!           the partition.                                       *      
!  M      : M describes the number of sub-intervals [Z(I),Z(I+1)]*      
!           I=0, 1, ..., M-1 with M >= 1.                        *      
!  WORK   : 2-dimensional array  WORK(0:N,2) containing the N+1  *      
!           weights of the CLENSHAW-CURTIS formula in its 1st    *      
!           column and the N+1 nodes in the 2nd column.          *      
!           In the calling program, the  SUBROUTINE WGKNOT must  *      
!           be used before calling this subroutine in order that *      
!           the weights and nodes become available.              *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  QCCN   : approximate value for the integral.                  *      
!  IERR   : error parameter.                                     *      
!           IERR = 0: everything o.k.                            *      
!           IERR = 1: condition for N is not met                 *      
!           IERR = 2: condition for M is not met                 *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
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
      DIMENSION Z (0:M), WORK (0:N, 2) 
!                                                                       
!  testing the input for correctness                                    
!                                                                       
      IERR = 1 
      IF (N.LT.2.OR.MOD (N, 2) .NE.0) RETURN 
      IERR = 2 
      IF (M.LT.1) RETURN 
!                                                                       
!  determine the approximate value for the integral                     
!                                                                       
      IERR = 0 
      QCCN = 0.0D0 
      DO 10 I = 0, M - 1 
         DUMMY1 = 0.5D0 * (Z (I + 1) - Z (I) ) 
         DUMMY2 = 0.5D0 * (Z (I + 1) + Z (I) ) 
         SUM = 0.0D0 
         DO 20 K = 0, N 
            SUM = SUM + WORK (K, 1) * FCT (DUMMY1 * WORK (K, 2) +       &
            DUMMY2)                                                     
   20    END DO 
         QCCN = QCCN + SUM * DUMMY1 
!                                                                       
   10 END DO 
      RETURN 
      END SUBROUTINE CLENSH                         
