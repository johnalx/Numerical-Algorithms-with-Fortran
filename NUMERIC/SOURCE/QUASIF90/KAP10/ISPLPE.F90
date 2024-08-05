      SUBROUTINE ISPLPE (N, XN, FN, MREP, B, C, D, H, DU, DM, RC, RS,   &
      IERR)                                                             
!                                                                       
!*******************************************************************    
!                                                                  *    
!  ISPLPE computes the coefficients B(I), C(I), D(I) for I=0,1,.., *    
!  N-1 of a cubic periodic interpolating spline.                   *    
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
!                 3.     FN(0) = FN(N)                             *    
!                                                                  *    
!                                                                  *    
!  NOTE:  ISPLPE should not be used by itself, but rather via the  *    
!  =====  SUBROUTINE ISPLNP, or ISPLPA for parametric or ISPLTR for*    
!         transformed parametric splines.                          *    
!         These subroutines also check the input data.             *    
!                                                                  *    
!                                                                  *    
!  INPUT PARAMETERS:                                               *    
!  =================                                               *    
!  N  :  index of the final node                                   *    
!  XN :  vector XN(0:N); the nodes XN(I), I = 0,1,...,N            *    
!  FN :  vector FN(0:N); the functional values FN(I) = FN( XN(I) ) *    
!                                                                  *    
!                                                                  *    
!  MREP  : index for repeated call of the SUBROUTINE:              *    
!          MREP = 1: In order to compute C(I) the system matrix and*    
!                    its factorization must be evaluated in SUB-   *    
!                    ROUTINE CYTSY anew.                           *    
!          MREP = 2: Only the right hand side of the system needs  *    
!                    to be recomputed. The vectors RC, DU and DM as*    
!                    computed during the first call in SUBROUTINE  *    
!                    CYTSYS help us find the solution.             *    
!                    (This prevents a duplicate decomposition in   *    
!                    case of parametric splines).                  *    
!                    The entries of H, RC, DU and DM must not be   *    
!                    altered after the first call.                 *    
!                                                                  *    
!                                                                  *    
!  AUXILIARY VARIABLES:                                            *    
!  ====================                                            *    
!  H  :   N-vector H(0:N-1)                                        *    
!  DU : ]                                                          *    
!  DM : ] N-1-vectors ..(1:N-1)                                    *    
!  RS : ]                                                          *    
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
!          = -4 :  FN(0) is not equal to FN(N)                     *    
!          = -5 :  erroneous value for MREP                        *    
!          =  1 :  crash in SUBROUTINE CYTSY, system matrix        *    
!                  numerically singular                            *    
!                                                                  *    
!------------------------------------------------------------------*    
!                                                                  *    
!  Subroutines required: CYTSY, CYTSYS                             *    
!                                                                  *    
!                                                                  *    
!                                                                  *    
!*******************************************************************    
!                                                                  *    
!  Author   : Gnter Palm                                          *    
!  Date     : 15.04.1988                                           *    
!  Source   : FORTRAN 77                                           *    
!                                                                  *    
!*******************************************************************    
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION XN (0:N), FN (0:N), B (0:N), C (0:N), D (0:N),    &
      H (0:N), DU (1:N), DM (1:N), RC (1:N), RS (1:N)                   
!                                                                       
!-----check periodicity of FN                                           
!                                                                       
      IERR = - 4 
      IF (FN (N) .NE.FN (0) ) RETURN 
!                                                                       
!-----check subsequent call                                             
!                                                                       
      IERR = - 5 
      IF (MREP.NE.1.AND.MREP.NE.2) RETURN 
!                                                                       
!-----Compute auxiliarty variables and the matrix elements              
!     (co and main diagonal) for the system matrix on first call        
!                                                                       
      IF (MREP.EQ.1) THEN 
!                                                                       
!       Auxiliary variables                                             
!                                                                       
         DO 10 I = 0, N - 1, 1 
            H (I) = XN (I + 1) - XN (I) 
   10    END DO 
         H (N) = H (0) 
!                                                                       
!       co diagonal                                                     
!                                                                       
         DO 20 I = 1, N - 1, 1 
            DU (I) = H (I) 
   20    END DO 
         DU (N) = H (0) 
!                                                                       
!       main diagonal                                                   
!                                                                       
         DO 30 I = 1, N, 1 
            DM (I) = 2.0D0 * (H (I - 1) + H (I) ) 
   30    END DO 
      ENDIF 
!                                                                       
!-----Compute the right hand side                                       
!                                                                       
      DUMMY1 = (FN (1) - FN (0) ) / H (0) 
      DO 40 I = 1, N - 1, 1 
         DUMMY2 = (FN (I + 1) - FN (I) ) / H (I) 
         RS (I) = 3.0D0 * (DUMMY2 - DUMMY1) 
         DUMMY1 = DUMMY2 
   40 END DO 
      RS (N) = 3.0D0 * ( (FN (1) - FN (0) ) / H (0) - DUMMY1) 
!                                                                       
!-----Solve the linear system to find C(1) , ..., C(N-1)                
!                                                                       
      IF (MREP.EQ.1) THEN 
!                                                                       
!       ... factorize on first call                                     
!                                                                       
         CALL CYTSY (N, DM, DU, RC, RS, C (1), IFLAG) 
         IF (IFLAG.NE.1) THEN 
            IF (IFLAG.EQ. - 2) THEN 
               IERR = - 1 
            ELSE 
               IERR = 1 
            ENDIF 
            RETURN 
         ENDIF 
      ELSE 
!                                                                       
!       ... without factorization on subsequent calls                   
!                                                                       
         CALL CYTSYS (N, DM, DU, RC, RS, C (1) ) 
      ENDIF 
      IERR = 0 
!                                                                       
!-----Compute remaining spline coefficients                             
!                                                                       
      C (0) = C (N) 
!                                                                       
      DO 50 I = 0, N - 1, 1 
         B (I) = (FN (I + 1) - FN (I) ) / H (I) - H (I) / 3.0D0 *       &
         (C (I + 1) + 2.0D0 * C (I) )                                   
         D (I) = (C (I + 1) - C (I) ) / (3.0D0 * H (I) ) 
   50 END DO 
      RETURN 
      END SUBROUTINE ISPLPE                         
