      SUBROUTINE GAUSRS (N, A, LDA, M, RS, XL, MARK, D, IPIVOT) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Solving a linear systems of equations  A * XL = RS  for M     *      
!  right hand sides using the Gauss-elimination method with      *      
!  scaling and column pivot search .                             *      
!  If the system has the form                                    *      
!         A * XL = I  , where I = identity matrix and A, XL, I   *      
!  are all (NxN) matrices, then the solution XL is the matrix    *      
!  inverse of A.                                                 *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETER:                                              *      
!  ================                                              *      
!  N        : order of the system of equations.                  *      
!  A        : 2-dimensional array A(1:LDA,1:N), containing the   *      
!             LDAxN matrix A common to all M systems of equations*      
!             (A = A(ORG)).                                      *      
!  LDA      : leading dimension of A, RS and XL, as defined in   *      
!             the calling program.                               *      
!  M        : number of right hand sides and hence the number of *      
!             solution vectors.                                  *      
!  RS       : 2-dimensional array RS(1:LDA,1:M), that is formed  *      
!             with the M right hand sides as columns.            *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  A        : 2-dimensional array A(1:LDA,1:N), containing the   *      
!             factors L and R with  P * A(ORG) = L * R. Here     *      
!             P = permutation matrix, L = unit lower triangular  *      
!             matrix and R = upper triangular matrix.            *      
!  XL       : 2-dimensional array XL(1:LDA,1:M) that contains    *      
!             the M solution vectors as columns for each of the  *      
!             M systems of equations.                            *      
!  MARK     : = 1, even number of row permutations.              *      
!             =-1, odd number of row permutations.               *      
!             = 0, input matrix A is numerically singular.       *      
!             The determinant of A is :                          *      
!                DET(A(ORG)) = MARK * A(1,1) * ... * A(N,N).     *      
!  D        : N-vector D(1:N); the reciprocals of the row sum    *      
!             norms of A(ORG), used for scaling:                 *      
!             D(I) = 1./(ABS(A(I,1)) + ... + ABS(A(I,N)))  for   *      
!             I = 1, ..., N.                                     *      
!  IPIVOT   : N-vector IPIVOT(1:N); it indicates the row per-    *      
!             mutations and thus defines the permutation matrix  *      
!             P. If e.g. IPIVOT(2) = 7, then the 7th row in      *      
!             of A(ORG) will become the 2nd row of P * A(ORG).   *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: GAUSSP, GAUSSS                          *      
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
      DOUBLEPRECISION A (1:LDA, 1:N), RS (1:LDA, 1:M), XL (1:LDA, 1:M), &
      D (1:N)                                                           
      INTEGER IPIVOT (1:N) 
!                                                                       
!  Factoring the matrix A by applying SUBROUTINE GAUSSP.                
!                                                                       
      CALL GAUSSP (N, A, LDA, IPIVOT, MARK, D) 
!                                                                       
!  Updating and bachsubstitution using SUBROUTINE GAUSSS in order to    
!  calculate the solution vectors for the M systems of equations.       
!                                                                       
      IF (MARK.NE.0) THEN 
         DO 10 K = 1, M 
            CALL GAUSSS (N, A, LDA, IPIVOT, RS (1, K), XL (1, K) ) 
   10    END DO 
      ENDIF 
      RETURN 
      END SUBROUTINE GAUSRS                         
