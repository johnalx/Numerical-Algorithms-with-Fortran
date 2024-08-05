      SUBROUTINE ISPLNK (N, XN, FN, B, C, D, H, DM, RS, IERR) 
!                                                                       
!*******************************************************************    
!                                                                  *    
!  ISPLNK computes the coefficients B(I), C(I), D(I) for I=0,1,.., *    
!  N-1 of a cubic interpolating spline with 'not a node' end       *    
!  point derivatives.                                              *    
!                                                                  *    
!  The spline has the form:                                        *    
!                                                                  *    
!  S(X) = FN(I) + B(I)(X-XN(I)) + C(I)(X-XN(I))**2 +               *    
!                               + D(I)(X-XN(I))**3                 *    
!                                                                  *    
!  for X in the interval [XN(I),XN(I+1)] for I=0,1,...,N-1.        *    
!                                                                  *    
!                                                                  *    
!  ASSUMPTIONS:   1.         N > 2                                 *    
!  ============   2.     XN(I) < XN(I+1), I=0,1,...,N-1            *    
!                                                                  *    
!                                                                  *    
!  NOTE:  ISPLNK should not be used by itself, but rather via the  *    
!  =====  SUBROUTINE ISPLNP, which also checks the input data.     *    
!                                                                  *    
!                                                                  *    
!  INPUT PARAMETERS:                                               *    
!  =================                                               *    
!  N  :  index of the final node                                   *    
!  XN :  vector XN(0:N); the nodes XN(I), I = 0,1,...,N            *    
!  FN :  vector FN(0:N); the functional values FN(I) = FN( XN(I) ) *    
!                                                                  *    
!                                                                  *    
!  AUXILIARY VARIABLES:                                            *    
!  ====================                                            *    
!  H  :   N-vector H(0:N-1)                                        *    
!  DM :   N-1-vector DM(1:N-1)                                     *    
!  RS :   N-1-vector RS(1:N-1)                                     *    
!                                                                  *    
!                                                                  *    
!  OUTPUT PARAMETERS:                                              *    
!  ==================                                              *    
!  FN :  ]  N+1-vectors ..(0:N);                                   *    
!  B  :  ]  The first N entries of B, C and D are the spline       *    
!  C  :  ]  coefficients for S. B(N), C(N), D(N) are auxiliary     *    
!  D  :  ]  variables.                                             *    
!  IERR :  error parameter                                         *    
!          =  0 :  All is ok                                       *    
!          = -1 :  N < 3                                           *    
!          =  1 :  crash in SUBROUTINE TRDIG, system matrix        *    
!                  numerically singular                            *    
!                                                                  *    
!------------------------------------------------------------------*    
!                                                                  *    
!  Subroutines required: TRDIG                                     *    
!                                                                  *    
!                                                                  *    
!                                                                  *    
!*******************************************************************    
!                                                                  *    
!  Author   : Gnter Palm                                          *    
!  Date     : 04.15.1988                                           *    
!  Source   : FORTRAN 77                                           *    
!                                                                  *    
!*******************************************************************    
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION XN (0:N), FN (0:N), B (0:N), C (0:N), D (0:N),    &
      H (0:N - 1), DM (1:N - 1), RS (1:N - 1)                           
!                                                                       
!-----Computing auxiliary variables                                     
!                                                                       
      DO 10 I = 0, N - 1, 1 
         H (I) = XN (I + 1) - XN (I) 
   10 END DO 
!                                                                       
!-----Compute the main and both co diagonals of the system              
!     matrix and the right hand side of A*C=RS with a                   
!     tridiagonal matrix A                                              
!                                                                       
!     co diagonals                                                      
!                                                                       
      D (1) = H (1) - H (0) 
      B (2) = H (0) 
      DO 20 I = 2, N - 3, 1 
         D (I) = H (I) 
         B (I + 1) = H (I) 
   20 END DO 
      D (N - 2) = H (N - 2) 
      B (N - 1) = H (N - 2) - H (N - 1) 
!                                                                       
!     main diagonal                                                     
!                                                                       
      DM (1) = H (0) + 2.0D0 * H (1) 
      DO 30 I = 2, N - 2, 1 
         DM (I) = 2.0D0 * (H (I - 1) + H (I) ) 
   30 END DO 
      DM (N - 1) = 2.0D0 * H (N - 2) + H (N - 1) 
!                                                                       
!     right hand side                                                   
!                                                                       
      DUMMY1 = (FN (2) - FN (1) ) / H (1) 
      RS (1) = 3.0D0 * H (1) / (H (1) + H (0) ) * (DUMMY1 - (FN (1)     &
      - FN (0) ) / H (0) )                                              
      DO 40 I = 2, N - 2, 1 
         DUMMY2 = (FN (I + 1) - FN (I) ) / H (I) 
         RS (I) = 3.0D0 * (DUMMY2 - DUMMY1) 
         DUMMY1 = DUMMY2 
   40 END DO 
      RS (N - 1) = 3.0D0 * H (N - 2) / (H (N - 2) + H (N - 1) ) *       &
      ( (FN (N) - FN (N - 1) ) / H (N - 1) - DUMMY1)                    
!                                                                       
!-----Solve the linear system to find C(1), ..., C(N-1)                 
!                                                                       
      CALL TRDIG (N - 1, B (1), DM, D (1), RS, C (1), IFLAG) 
      IF (IFLAG.NE.1) THEN 
         IF (IFLAG.EQ.0) THEN 
            IERR = 1 
         ELSE 
            IERR = - 1 
         ENDIF 
      ENDIF 
      IERR = 0 
!                                                                       
!-----Compute the remaining spline coefficients                         
!                                                                       
      C (0) = C (1) + H (0) / H (1) * (C (1) - C (2) ) 
      C (N) = C (N - 1) + H (N - 1) / H (N - 2) * (C (N - 1) - C (N - 2)&
      )                                                                 
!                                                                       
      DO 50 I = 0, N - 1, 1 
         B (I) = (FN (I + 1) - FN (I) ) / H (I) - H (I) / 3.0D0 *       &
         (C (I + 1) + 2.0D0 * C (I) )                                   
         D (I) = (C (I + 1) - C (I) ) / (3.0D0 * H (I) ) 
   50 END DO 
      RETURN 
      END SUBROUTINE ISPLNK                         
