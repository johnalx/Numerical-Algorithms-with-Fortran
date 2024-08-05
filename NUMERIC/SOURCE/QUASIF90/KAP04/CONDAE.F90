      SUBROUTINE CONDAE (N, AO, A, IA, IPIVOT, Y, X, Z, R, CONDA) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  The SUBROUTINE CONDAE finds an estimate of the matrix         *      
!  condition number                                              *      
!          COND(A) =  Max norm of A * Max norm of A inverse,     *      
!  following the condition estimate of  FORSYTHE and MOLER.      *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N        : order of the matrices A and AO.                    *      
!  AO       : DOUBLE PRECISION array AO(1:IA,1:N) containing the *      
!             orginal matrix A = A(ORG).                         *      
!  A        : DOUBLE PRECISION array A(1:IA,1:N), containing the *      
!             output matrix of the SUBROUTINE GAUSSP, i.e., the  *      
!             factors L and R with  P*A(ORG)=L*R for a           *      
!             permutation matrix P.                              *      
!  IA       : leading dimension of A and AO, as stipulated by    *      
!             the calling program.                               *      
!  IPIVOT   : INTEGER N-vector IPIVOT(1:N), output of the        *      
!             SUBROUTINE GAUSS, i.e., information on P.          *      
!  Y        : DOUBLE PRECISION vector Y(1:N) containing the right*      
!             hand side of the linear system  A(ORG)*X=Y.        *      
!  X        : DOUBLE PRECISION vector X(1:N) with the solution   *      
!             of the linear system from  GAUSS.                  *      
!                                                                *      
!                                                                *      
!  AUXILIARY PARAMETERS:                                         *      
!  =====================                                         *      
!  Z        : ]  DOUBLE PRECISION vectors ..(1:N).               *      
!  R        : ]                                                  *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  CONDA    : Estimate for  COND(A).                             *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines:  GAUSSS, MACHPD                         *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author    : Gisela Engeln-MÅllges                             *      
!  Date      : 06.07.90                                          *      
!  Source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
! Declarations                                                          
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION AO (1:IA, 1:N), A (1:IA, 1:N), Y (1:N), X (1:N),  &
      Z (1:N), R (1:N)                                                  
      INTEGER IPIVOT (1:N) 
!                                                                       
! Compute the max norm of the solution  X  of  A X = Y                  
! using the Gau· algorithm                                              
!                                                                       
      XMAX = DABS (X (1) ) 
      DO 10 I = 2, N 
         XMAX = DMAX1 (XMAX, DABS (X (I) ) ) 
   10 END DO 
!                                                                       
! Find the machine constant                                             
!                                                                       
      FMACHP = 1.0D0 
   20 FMACHP = 0.5D0 * FMACHP 
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 20 
      EPS = 2.0D0 * FMACHP 
!                                                                       
! Find the residuum  Y - A * X; form A * X  in double precision,        
! then round the residuum to single precision                           
!                                                                       
      DO 30 I = 1, N 
         R (I) = Y (I) 
         DO 40 K = 1, N 
            R (I) = R (I) - AO (I, K) * X (K) 
   40    END DO 
   30 END DO 
!                                                                       
! Calculate the first correction Z to the solution  and its max norm    
!                                                                       
      CALL GAUSSS (N, A, IA, IPIVOT, R, Z) 
      ZMAX = DABS (Z (1) ) 
      DO 50 I = 2, N 
         ZMAX = DMAX1 (ZMAX, DABS (Z (I) ) ) 
   50 END DO 
!                                                                       
! Estimate the condition number  COND(A)=NORM(A)*NORM(INV(A))           
! from the sizes of Z, X and EPS                                        
!                                                                       
      CONDA = ZMAX / (XMAX * EPS) 
!                                                                       
      RETURN 
      END SUBROUTINE CONDAE                         
