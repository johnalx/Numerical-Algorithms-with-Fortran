![             Condition Number}*)                                      
      SUBROUTINE HACOND (N, A0, A, LDA, MARK, HCOND) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  SUBROUTINE HACOND calculates the Hadamard condition number    *      
!  of the matrix A0 with A0 = A(ORG).  The determinant of A0 is  *      
!  calculated via the product of the diagonal elements of the    *      
!  upper triangular factor R from SUBROUTINE GAUSSP, where       *      
!    P * A(ORG) = L * R, for  P a permutation matrix.            *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N        : order of the square matrices A and A0.             *      
!  A0       : 2-dimensional array A0(1:LDA,1:N); the matrix      *      
!             A(ORG).                                            *      
!  A        : 2-dimensional array A(1:LDA,1:N)  containing the   *      
!             factors  L  and  R  with  P * A(ORG) = L * R.      *      
!             P = permutation array; A is overwritten by L and R.*      
!             It is the output array of SUBROUTINE GAUSSP.       *      
!  LDA      : leading dimension of A and A0 as defined in the    *      
!             calling program.                                   *      
!  MARK     : = 1, even number of row permutations.              *      
!             =-1, odd number of row permutations.               *      
!             = 0, matrix A is singular. A is the output of      *      
!             SUBROUTINE GAUSSP. The determinant is given as :   *      
!                DET(A(ORG)) = MARK * A(1,1) * ... * A(N,N).     *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  HCOND    : Hadamard condition number of A0.                   *      
!             empirical advice: If                               *      
!             HCOND < 0.01         : badly conditioned matrix A0 *      
!                                    (The smaller HCOND, the     *      
!                                    worse A0 is conditiond).    *      
!             0.01 <= HCOND <= 0.1 : no precise conclusion.      *      
!             HCOND > 0.1          : well conditioned matrix A0. *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  authors   : Gisela Engeln-Muellges, Guido Dubois              *      
!  date      : 04.25.88                                          *      
!  source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION A0 (1:LDA, 1:N), A (1:LDA, 1:N) 
      HCOND = 0.0D0 
      IF (MARK.EQ.0) RETURN 
      HCOND = 1.0D0 
      DO 10 I = 1, N 
         ZSNORM = 0.0D0 
         DO 20 K = 1, N 
            ZSNORM = ZSNORM + A0 (I, K) * A0 (I, K) 
   20    END DO 
         HCOND = HCOND * A (I, I) / DSQRT (ZSNORM) 
   10 END DO 
      HCOND = DABS (HCOND) 
      RETURN 
      END SUBROUTINE HACOND                         
