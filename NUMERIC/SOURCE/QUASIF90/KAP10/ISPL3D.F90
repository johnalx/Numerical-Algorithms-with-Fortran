      SUBROUTINE ISPL3D (N, XN, FN, ALPHA, BETA, B, C, D, H, RS, IERR) 
!                                                                       
!*******************************************************************    
!                                                                  *    
!  ISPL3D computes the coefficients B(I), C(I), D(I) for I=0,1,.., *    
!  N-1 of a cubic interpolating spline with prescribed third end   *    
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
!  NOTE:  ISPL3D should not be used by itself, but rather via the  *    
!  =====  SUBROUTINE ISPLNP, which also checks the input data.     *    
!                                                                  *    
!                                                                  *    
!  INPUT PARAMETERS:                                               *    
!  =================                                               *    
!  N  :  index of the final node                                   *    
!  XN :  vector XN(0:N); the nodes XN(I), I = 0,1,...,N            *    
!  FN :  vector FN(0:N); the functional values FN(I) = FN( XN(I) ) *    
!                                                                  *    
!  ALPHA :  first derivative at XN(0)                              *    
!  BETA  :  first derivative at XN(N)                              *    
!                                                                  *    
!                                                                  *    
!  AUXILIARY VARIABLES:                                            *    
!  ====================                                            *    
!  H  :   N-vector H(0:N-1)                                        *    
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
!          =  1 :  crash in SUBROUTINE TRDSY, system matrix        *    
!                  numerically singular                            *    
!                                                                  *    
!------------------------------------------------------------------*    
!                                                                  *    
!  Subroutines required: TRDSY                                     *    
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
      H (0:N - 1), RS (1:N - 1)                                         
!                                                                       
!-----Compute auxiliary variables                                       
!                                                                       
      DO 10 I = 0, N - 1, 1 
         H (I) = XN (I + 1) - XN (I) 
   10 END DO 
!                                                                       
!-----Compute the main and co diagonal elements of the system           
!     matrix and the right hand side of A*C=RS with a                   
!     symmetric and tridiagonal matrix A                                
!                                                                       
!     co diagonal                                                       
!                                                                       
      DO 20 I = 1, N - 2, 1 
         D (I) = H (I) 
   20 END DO 
!                                                                       
!     main diagonal                                                     
!                                                                       
      B (1) = 3.0D0 * H (0) + 2.0D0 * H (1) 
      DO 30 I = 2, N - 2, 1 
         B (I) = 2.0D0 * (H (I - 1) + H (I) ) 
   30 END DO 
      B (N - 1) = 2.0D0 * H (N - 2) + 3.0D0 * H (N - 1) 
!                                                                       
!     right hand side                                                   
!                                                                       
      C (0) = 0.5D0 * ALPHA * H (0) 
      C (N) = 0.5D0 * BETA * H (N - 1) 
!                                                                       
      DUMMY1 = (FN (2) - FN (1) ) / H (1) 
      RS (1) = 3.0D0 * (DUMMY1 - (FN (1) - FN (0) ) / H (0) ) + C (0)   &
      * H (0)                                                           
      DO 40 I = 2, N - 2, 1 
         DUMMY2 = (FN (I + 1) - FN (I) ) / H (I) 
         RS (I) = 3.0D0 * (DUMMY2 - DUMMY1) 
         DUMMY1 = DUMMY2 
   40 END DO 
      RS (N - 1) = 3.0D0 * ( (FN (N) - FN (N - 1) ) / H (N - 1) -       &
      DUMMY1) - C (N) * H (N - 1)                                       
!                                                                       
!-----Solve the linear system to obtain C(1), ... ,C(N-1)               
!                                                                       
      CALL TRDSY (N - 1, B (1), D (1), RS, C (1), IFLAG) 
      IF (IFLAG.NE.1) THEN 
         IF (IFLAG.EQ. - 2) THEN 
            IERR = - 1 
         ELSE 
            IERR = 1 
         ENDIF 
         RETURN 
      ENDIF 
      IERR = 0 
!                                                                       
!-----Compute the remaining spline coefficients                         
!                                                                       
      C (0) = C (1) - C (0) 
      C (N) = C (N - 1) + C (N) 
!                                                                       
      DO 50 I = 0, N - 1, 1 
         B (I) = (FN (I + 1) - FN (I) ) / H (I) - H (I) / 3.0D0 *       &
         (C (I + 1) + 2.0D0 * C (I) )                                   
         D (I) = (C (I + 1) - C (I) ) / (3.0D0 * H (I) ) 
   50 END DO 
      RETURN 
      END SUBROUTINE ISPL3D                         
