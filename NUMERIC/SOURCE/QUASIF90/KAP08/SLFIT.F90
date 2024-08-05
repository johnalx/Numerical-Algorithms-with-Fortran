![  {Householder Transformation for Linear Least Squares}               
![  {Solving Linear Least Squares Problems using Householder            
![   Transformations}*)                                                 
      SUBROUTINE SLFIT (X, Y, W, IWFL, FCT, PSI, LDA, M, N, INUM, A, D, &
      C, FV, SQERR, IERR)                                               
!                                                                       
!*****************************************************************      
!                                                                *      
!  SUBROUTINE SLFIT solves the linear discrete least squares     *      
!  problem for M+1 nodes (X(I), Y(I)) by computing a linear      *      
!  approximation function                                        *      
!         F(X) = C(0) * F1(X) + ... + C(N) * FN(X).              *      
!  The model functions F1(X), ..., FN(X) must be provided by the *      
!  user in form of a SUBROUTINE.                                 *      
!  SLFIT determines the optimal coefficients C(I), I=0, ..., N   *      
!  in the least squares sense.                                   *      
!  The problem is reduced to a linear minimization problem which *      
!  is solved by applying Householder transformations to the      *      
!  system matrix of the normal equations.                        *      
!  For INUM user-given values the approximating function is then *      
!  evaluated.                                                    *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!                                                                *      
!  X      (N+1)-vector X(0:M) containing the X-values of the     *      
!         nodes                                                  *      
!  Y      (N+1)-vector Y(0:M) containing the Y-values at the     *      
!         nodes                                                  *      
!  W      (N+1)-vector W(0:M) containing the positive weights    *      
!  IWFL   if IWFL = 0, the nodes are weighed according to W      *      
!         Otherwise the nodes are all equally weighed by 1, and  *      
!         W does not need to be defined explicitly.              *      
!  FCT    A SUBROUTINE that has to be provided by the user. It   *      
!         defines the model functions and has the form:          *      
!                                                                *      
!         SUBROUTINE FCT (X, N, F)                               *      
!         IMPLICIT DOUBLE PRECISION (A-H,O-Z)                    *      
!         INTEGER N                                              *      
!         DIMENSION F(0:N), X                                    *      
!         -----------------                                      *      
!         F(0) = value for the 0-th model function at X          *      
!         .                                                      *      
!         .                                                      *      
!         .                                                      *      
!         F(N) = value of the N-th model function at X           *      
!         -----------------                                      *      
!         RETURN                                                 *      
!         END                                                    *      
!                                                                *      
!         In the calling program FCT has to be declared as       *      
!         EXTERNAL.                                              *      
!  PSI    INUM-vector PSI(1:INUM) containing the values where    *      
!         the approximating function is to be evaluated.         *      
!  LDA    leading dimension of array A as defined in the calling *      
!         program; LDA has to be >= M.                           *      
!  M      M+1 = number of nodes                                  *      
!  N      N+1 = number of model functions                        *      
!  INUM   for INUM values the approximation function is to be    *      
!         evaluated.                                             *      
!                                                                *      
!                                                                *      
!  AUXILIARY PARAMETERS:                                         *      
!  =====================                                         *      
!                                                                *      
!  A      2-dim. array A(0:LDA,0:N+1) for SUBROUTINE SHOUSE      *      
!  D      (N+1)-vector D(0:N) for SUBROUTINE SHOUSE              *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!                                                                *      
!  C      (N+1)-vector C(0:N) containing the optimal coefficients*      
!         with respect to the least square error                 *      
!  FV     INUM-vector FV(1:INUM) containing the values  of the   *      
!         approximating function at the INUM given locations PSI.*      
!  IERR  error parameter                                         *      
!         = 0  everything  o.k.                                  *      
!         = 1  error in the input parameters                     *      
!         = 2  approximating function cannot be determined,      *      
!              since numerically A does not have maximal rank.   *      
!              I.e., the model functions are numerically linearly*      
!              dependent.                                        *      
!  SQERR least square error                                      *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: SLPRE, SHOUSE, SENORM                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Ilona Westermann                                   *      
!  date     : 09.01.1987                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION X (0:M), Y (0:M), A (0:LDA, 0:N + 1), D (0:N), C (0:N), &
      PSI (INUM), FV (INUM), W (0:M)                                    
      EXTERNAL FCT 
      IERR = 0 
!                                                                       
!     checking the input parameters                                     
!                                                                       
      IF (M.GE.N.AND.LDA.GE.M.AND.N.GE.0.AND.M.GE.0) THEN 
!                                                                       
!        minimize (A * X - B) via the following steps:                  
!                                                                       
!        1. determine the matrix A in SUBROUTINE SLPRE                  
!                                                                       
         CALL SLPRE (X, W, IWFL, FCT, LDA, M, N, D, A) 
!                                                                       
!        2. store the vector B (if necessary modified by the            
!           weights) in the last column of A                            
!                                                                       
         IF (IWFL.EQ.0) THEN 
            DO 20 I = 0, M 
               A (I, N + 1) = Y (I) * DSQRT (W (I) ) 
   20       END DO 
         ELSE 
            DO 30 I = 0, M 
               A (I, N + 1) = Y (I) 
   30       END DO 
         ENDIF 
!                                                                       
!        3. solve the minimization problem using the SUBROUTINE SHOUSE  
!                                                                       
         CALL SHOUSE (A, LDA, M, N, D, C, MARK) 
!                                                                       
!        test for singularity                                           
!                                                                       
         IF (MARK.EQ.0) THEN 
                                                                        
!                                                                       
!           determine the value of the approximating function           
!           at the locations desired                                    
!                                                                       
            DO 40 I = 1, INUM 
               CALL FCT (PSI (I), N, D) 
               FV (I) = 0.0D0 
               DO 40 J = 0, N 
                  FV (I) = FV (I) + C (J) * D (J) 
   40       CONTINUE 
!                                                                       
!           determine least square error                                
!                                                                       
            SQERR = DSQRT (SENORM (A (N + 1, N + 1), M - (N + 1) ) ) 
         ELSE 
            IERR = 2 
         ENDIF 
      ELSE 
         IERR = 1 
      ENDIF 
      RETURN 
      END SUBROUTINE SLFIT                          
