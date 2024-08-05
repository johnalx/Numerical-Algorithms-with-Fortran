      SUBROUTINE CFSP3D (N, XN, FN, W, ALPHA, BETA, A, B, C, D, DL1,    &
      DL2, DU1, DU2, RS, H1, H, IERR)                                   
!                                                                       
!*****************************************************************      
!                                                                *      
!  CFSP3D computes the coefficients A(I), B(I), C(I), D(I),      *      
!  I=0, 1. ..., N-1, of a cubic fitting spline with prescribed   *      
!  third end point derivative.                                   *      
!  The spline is represented in the form:                        *      
!                                                                *      
!  S(X) = A(I) + B(I)(X-XN(I)) + C(I)(X-XN(I))**2 +              *      
!                              + D(I)(X-XN(I))**3                *      
!                                                                *      
!  for X in the interval [XN(I),XN(I+1)], I=0, 1, ..., N-1.      *      
!                                                                *      
!                                                                *      
!  ASSUMPTIONS:    1.         N > 4                              *      
!  ============    2.     XN(I) < XN(I+1), I=0, 1, ..., N-1      *      
!                  3.      W(I) > 0.0    , I=0, 1, ..., N        *      
!                                                                *      
!                                                                *      
!  REMARK:  CFSP3D should not be called directly, but rather via *      
!  =======  the subroutine CFSPNP. The subroutine CFSPNP also    *      
!           checks the above assumptions.                        *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N  :  Index of the last node                                  *      
!  XN :  vector XN(0:N); XN(I) is the Ith node, I = 0, ..., N    *      
!  FN :  vector FN(0:N); FN(I) is the data at the node XN(I)     *      
!  W  :  vector W(0:N);  W(I) is the weight of FN(I)             *      
!                                                                *      
!  ALPHA :  third end point derivative at XN(0)                  *      
!  BETA  :  third end point derivative at XN(N)                  *      
!                                                                *      
!                                                                *      
!  AUXILIARY VARIABLES:                                          *      
!  ====================                                          *      
!  H   :]                                                        *      
!  H1  :]  N-vectors H(0:N-1), H1(0:N-1)                         *      
!                                                                *      
!  DL1 :]                                                        *      
!  DL2 :]                                                        *      
!  DU1 :]  (N-1)-vectors dimensioned as  ..(1:N-1)               *      
!  DU2 :]                                                        *      
!  RS  :]                                                        *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  A    :  Vector A(0:N) ]  The entries in positions 0 to N-1    *      
!  B    :  Vector B(0:N) ]  contain the spline coefficients for  *      
!  C    :  Vector C(0:N) ]  S. The entries in A(N), B(N), C(N)   *      
!  D    :  Vector D(0:N) ]  and D(N) are auxiliary variables.    *      
!                                                                *      
!  IERR :  error parameter                                       *      
!          =  0 :  All is o.k.                                   *      
!          = -1 :  N < 5                                         *      
!          =  1 :  FDIAG did not run correctly (matrix singular) *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  required subroutines: FDIAG                                   *      
!                                                                *      
!                                                                *      
!  Reference: Engeln-Mllges, G.; Reutter, F., [ENGE87].         *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : Gnter Palm                                        *      
!  Date     : 04.18.1988                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION XN (0:N), FN (0:N), W (0:N), A (0:N), B (0:N),    &
      C (0:N), D (0:N), DL1 (1:N - 1), DL2 (1:N - 1), DU1 (1:N - 1),    &
      DU2 (1:N - 1), RS (1:N - 1), H1 (0:N - 1), H (0:N - 1)            
!                                                                       
!-----Computing the auxiliary variables                                 
!                                                                       
      DO 10 I = 0, N - 1, 1 
         H (I) = XN (I + 1) - XN (I) 
         H1 (I) = 1.0D0 / H (I) 
         C (I) = H1 (I) * H1 (I) 
         B (I) = 6.0D0 / W (I) 
   10 END DO 
      B (N) = 6.0D0 / W (N) 
!                                                                       
      DO 20 I = 0, N - 2, 1 
         D (I) = H1 (I) + H1 (I + 1) 
   20 END DO 
!                                                                       
!-----Compute the system matrix elements (main and two co-diagonals)    
!     and the right-hand side for  A*C=RS with a five-diagonal matrix A.
!                                                                       
!     second co-diagonal                                                
!                                                                       
      DO 30 I = 3, N - 1, 1 
         DL2 (I) = B (I - 1) * H1 (I - 2) * H1 (I - 1) 
         DU2 (I - 2) = DL2 (I) 
   30 END DO 
!                                                                       
!     first co-diagonal                                                 
!                                                                       
      DUMMY1 = H (1) - B (2) * H1 (1) * D (1) 
      DL1 (2) = DUMMY1 - B (1) * C (1) 
      DU1 (1) = DUMMY1 - B (1) * H1 (1) * D (0) 
      DO 40 I = 3, N - 2, 1 
         K = I - 1 
         DL1 (I) = H (K) - B (K) * H1 (K) * D (K - 1) - B (I) * H1 (K)  &
         * D (K)                                                        
         DU1 (K) = DL1 (I) 
   40 END DO 
      DUMMY1 = H (N - 2) - B (N - 2) * H1 (N - 2) * D (N - 3) 
      DL1 (N - 1) = DUMMY1 - B (N - 1) * H1 (N - 2) * D (N - 2) 
      DU1 (N - 2) = DUMMY1 - B (N - 1) * C (N - 2) 
!                                                                       
!     main diagonal                                                     
!                                                                       
      A (1) = 3.0D0 * H (0) + 2.0D0 * H (1) + B (1) * H1 (1) * D (0)    &
      + B (2) * C (1)                                                   
      DO 50 I = 2, N - 2, 1 
         K = I - 1 
         A (I) = 2.0D0 * (H (K) + H (I) ) + B (K) * C (K) + B (I)       &
         * D (K) * D (K) + B (I + 1) * C (I)                            
   50 END DO 
      A (N - 1) = 2.0D0 * H (N - 2) + 3.0D0 * H (N - 1) + B (N - 2)     &
      * C (N - 2) + B (N - 1) * H1 (N - 2) * D (N - 2)                  
!                                                                       
!     right-hand side                                                   
!                                                                       
      C (0) = 0.5D0 * ALPHA 
      C (N) = 0.5D0 * BETA 
!                                                                       
      DUMMY2 = (FN (2) - FN (1) ) * H1 (1) 
      DUMMY1 = (FN (3) - FN (2) ) * H1 (2) 
      RS (1) = 3.0D0 * (DUMMY2 - (FN (1) - FN (0) ) * H1 (0) ) + C (0)  &
      * (H (0) * H (0) - B (0) * H1 (0) - B (1) * D (0) )               
      RS (2) = 3.0D0 * (DUMMY1 - DUMMY2) + C (0) * B (1) * H1 (1) 
      DO 60 I = 3, N - 3, 1 
         DUMMY2 = (FN (I + 1) - FN (I) ) * H1 (I) 
         RS (I) = 3.0D0 * (DUMMY2 - DUMMY1) 
         DUMMY1 = DUMMY2 
   60 END DO 
      DUMMY2 = (FN (N - 1) - FN (N - 2) ) * H1 (N - 2) 
      RS (N - 2) = 3.0D0 * (DUMMY2 - DUMMY1) - C (N) * B (N - 1)        &
      * H1 (N - 2)                                                      
      RS (N - 1) = 3.0D0 * ( (FN (N) - FN (N - 1) ) * H1 (N - 1)        &
      - DUMMY2) - C (N) * (H (N - 1) * H (N - 1) - B (N - 1) * D (N - 2)&
      - B (N) * H1 (N - 1) )                                            
!                                                                       
!-----Compute the coefficients C(1) to C(N-1) by                        
!     solving the linear system                                         
!                                                                       
      CALL FDIAG (N - 1, DL2, DL1, A (1), DU1, DU2, RS, C (1), IFLAG) 
      IF (IFLAG.NE.1) THEN 
         IF (IFLAG.EQ.0) THEN 
            IERR = 1 
         ELSE 
            IERR = - 1 
         ENDIF 
         RETURN 
      ENDIF 
      IERR = 0 
!                                                                       
!-----Computing the remaining spline coefficients                       
!                                                                       
      C (0) = C (1) - C (0) * H (0) 
      C (N) = C (N - 1) + C (N) * H (N - 1) 
!                                                                       
      A (0) = FN (0) + B (0) / 3.0D0 * H1 (0) * (C (0) - C (1) ) 
      DO 70 I = 1, N - 1, 1 
         A (I) = FN (I) - B (I) / 3.0D0 * (C (I - 1) * H1 (I - 1)       &
         - D (I - 1) * C (I) + C (I + 1) * H1 (I) )                     
   70 END DO 
      A (N) = FN (N) - B (N) / 3.0D0 * H1 (N - 1) * (C (N - 1) - C (N) ) 
!                                                                       
      DO 80 I = 0, N - 1, 1 
         B (I) = H1 (I) * (A (I + 1) - A (I) ) - H (I) / 3.0D0 *        &
         (C (I + 1) + 2.0D0 * C (I) )                                   
         D (I) = H1 (I) / 3.0D0 * (C (I + 1) - C (I) ) 
   80 END DO 
      RETURN 
      END SUBROUTINE CFSP3D                         
