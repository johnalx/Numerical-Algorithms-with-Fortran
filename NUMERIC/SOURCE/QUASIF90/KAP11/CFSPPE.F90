      SUBROUTINE CFSPPE (N, XN, FN, W, MREP, A, B, C, D, H, H1, H2, H3, &
      RS, WORK, IERR)                                                   
!                                                                       
!*****************************************************************      
!                                                                *      
!  'CFSPPE' computes the coefficients A(I), B(I), C(I), D(I) for *      
!  I = 0, 1, ..., N-1 of a periodic cubic least squares spline.  *      
!  The spline function is represented in the form:               *      
!                                                                *      
!  S(X) = A(I) + B(I)(X-XN(I)) + C(I)(X-XN(I))**2 +              *      
!                              + D(I)(X-XN(I))**3                *      
!                                                                *      
!  for X in the interval [XN(I),XN(I+1)] for I=0,1,...,N-1.      *      
!                                                                *      
!                                                                *      
!  ASSUMPTIONS:        1.         N > 5                          *      
!  ================    2.     XN(I) < XN(I+1), I=0,1,...,N-1     *      
!                      3.      W(I) > 0.0    , I=0,1,...,N       *      
!                      4.      W(0) = W(N)                       *      
!                      5.     FN(0) = FN(N)                      *      
!                                                                *      
!                                                                *      
!  BEMERKUNG:  'CFSPPE' should not be called by itself, but      *      
!  ==========  rather using the subroutine 'CFSPNP' - or in      *      
!              case of parametric or transformed parametric      *      
!              splines via the subroutines 'CFSPPA' or 'CFSPTR'. *      
!              These subroutines check the assumptions 1 to 3.   *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N  :  Index of final node                                     *      
!  XN :  (N+1)-vector XN(0:N); XN(I) = nodes, I = 0,1,...,N      *      
!  FN :  (N+1)-vector FN(0:N); FN(I) = function values at XN(I)  *      
!  W  :  (N+1)-vector W(0:N);   W(I) = weights for FN(I)         *      
!                                                                *      
!  MREP :  Marker for repeated call of this subroutine:          *      
!          MREP = 1: We must form the system matrix and its      *      
!                    complete LR decomposition in order to       *      
!                    compute C(I) in the subroutine NCYFSY.      *      
!          MREP = 2: We only have to compute the right hand side.*      
!                    We can use the vector WORK from our first   *      
!                    call of NCYFSS to find the solution.        *      
!                    (This avoids a duplicate LR decomposition   *      
!                    for parametric splines). If MREP = 2, the   *      
!                    entries of the vectors H, H1, H2, H3 and    *      
!                    WORK must not be altered after the first    *      
!                    call of NCYFSS !                            *      
!  REMARK:    For parametric splines with differing weights      *      
!             (WX(I) different from WY(I) for at least one index *      
!             I) put MREP = 1 always!                            *      
!                                                                *      
!                                                                *      
!  AUX VARIABLES:                                                *      
!  ==============                                                *      
!  H    :]                                                       *      
!  H1   :]  (N+1)-vectors ..(0:N)                                *      
!  H2   :]                                                       *      
!  H3   :]                                                       *      
!                                                                *      
!  RS   :  N-vector RS(1:N)                                      *      
!  WORK :  vector WORK(1:9*N-14)                                 *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  A    :  ] (N+1)-vectors    The elements in positions 0 to N-1 *      
!  B    :  ] ..(0:N)          are the coefficients of the spline *      
!  C    :  ]                  function S.                        *      
!  D    :  ]                  A(N), B(N), C(N), D(N) are aux     *      
!                             variables.                         *      
!  IERR :  error parameter                                       *      
!          =  0 :  All ok                                        *      
!          = -1 :  N < 6                                         *      
!          = -5 :  FN(0) not equal to FN(N),  or                 *      
!                   W(0) not equal to W(N)                       *      
!          = -6 :  wrong value assigned to MREP                  *      
!          =  1 :  fatal error in NCYFSY                         *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Subroutines required : NCYFSY, NCYFSS                         *      
!                                                                *      
!                                                                *      
!  Reference : Engeln-Mllges, G.; Reutter, F., see [ENGE87].    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author      : Gnter Palm                                     *      
!  Date        : 18.04.1988                                      *      
!  Source code : FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION XN (0:N), FN (0:N), W (0:N), A (0:N), B (0:N),    &
      C (0:N), D (0:N), H (0:N), H1 (0:N), H2 (0:N), H3 (0:N), RS (1:N),&
      WORK (1:9 * N - 14)                                               
!                                                                       
!-----Check input for periodicity---------------                        
!                                                                       
      IERR = - 5 
      IF (FN (N) .NE.FN (0) ) RETURN 
      IF (W (N) .NE.W (0) ) RETURN 
!                                                                       
!-----Check marker for a repeated call---------------                   
!                                                                       
      IERR = - 6 
      IF (MREP.NE.1.AND.MREP.NE.2) RETURN 
!                                                                       
!-----Compute aux variables and matrix elements-----------              
!     Main diagonal, lower, and upper co-diagonals of system matrix     
!     for linear system  (this matrix is symmetric, almost cyclic       
!     and five-diagonal),  in case of a first call.                     
!                                                                       
      IF (MREP.EQ.1) THEN 
!                                                                       
!       Aux variables                                                   
!                                                                       
         DO 10 I = 0, N - 1, 1 
            H (I) = XN (I + 1) - XN (I) 
            H1 (I) = 1.0D0 / H (I) 
            C (I) = H1 (I) * H1 (I) 
            H2 (I) = 6.0D0 / W (I) 
   10    END DO 
         H (N) = H (0) 
         H1 (N) = H1 (0) 
         C (N) = C (0) 
         H2 (N) = H2 (0) 
!                                                                       
         DO 20 I = 0, N - 1, 1 
            H3 (I) = H1 (I) + H1 (I + 1) 
   20    END DO 
!                                                                       
!       Upper co-diagonal                                               
!                                                                       
         DO 30 I = 1, N - 1, 1 
            D (I) = H2 (I + 1) * H1 (I) * H1 (I + 1) 
   30    END DO 
         D (N) = H2 (1) * H1 (0) * H1 (1) 
!                                                                       
!       Lower co-diagonal                                               
!                                                                       
         DO 40 I = 1, N - 1, 1 
            B (I) = H (I) - H2 (I) * H1 (I) * H3 (I - 1) - H2 (I + 1)   &
            * H1 (I) * H3 (I)                                           
   40    END DO 
         B (N) = H (0) - H2 (0) * H1 (0) * H3 (N - 1) - H2 (1) * H1 (0) &
         * H3 (0)                                                       
!                                                                       
!       Main diagonal                                                   
!                                                                       
         DO 50 I = 1, N - 1, 1 
            K = I - 1 
            A (I) = 2.0D0 * (H (K) + H (I) ) + H2 (K) * C (K) + H2 (I)  &
            * H3 (K) * H3 (K) + H2 (I + 1) * C (I)                      
   50    END DO 
         A (N) = 2.0D0 * (H (N - 1) + H (N) ) + H2 (N - 1) * C (N - 1)  &
         + H2 (N) * H3 (N - 1) * H3 (N - 1) + H2 (1) * C (0)            
      ENDIF 
!                                                                       
!-----Compute right hand side--------------------------------           
!                                                                       
      DUMMY1 = (FN (1) - FN (0) ) * H1 (0) 
      DO 60 I = 1, N - 1, 1 
         DUMMY2 = (FN (I + 1) - FN (I) ) * H1 (I) 
         RS (I) = 3.0D0 * (DUMMY2 - DUMMY1) 
         DUMMY1 = DUMMY2 
   60 END DO 
      RS (N) = 3.0D0 * ( (FN (1) - FN (0) ) * H1 (0) - DUMMY1) 
!                                                                       
!-----Find coefficients C(1) to C(N) by------------                     
!     solving the linear system ...                                     
!                                                                       
      IF (MREP.EQ.1) THEN 
!                                                                       
!       ... form LR decomposition                                       
!                                                                       
         CALL NCYFSY (N, A (1), B (1), D (1), RS, C (1), WORK (1),      &
         WORK (N + 1), WORK (2 * N + 1), WORK (3 * N + 1), WORK (4 * N -&
         3), WORK (5 * N - 6), WORK (6 * N - 6), WORK (7 * N - 6),      &
         WORK (8 * N - 10), IERR)                                       
         IF (IERR.NE.0) RETURN 
      ELSE 
!                                                                       
!       ... without decomposing anew                                    
!                                                                       
         IERR = 0 
         CALL NCYFSS (N, RS, C (1), WORK (1), WORK (N + 1), WORK (2 * N &
         + 1), WORK (3 * N + 1), WORK (4 * N - 3), WORK (5 * N - 6),    &
         WORK (6 * N - 6), WORK (7 * N - 6), WORK (8 * N - 10) )        
      ENDIF 
!                                                                       
!-----Compute remaining spline coefficients---------------              
!                                                                       
      C (0) = C (N) 
!                                                                       
      A (0) = FN (0) - H2 (0) / 3.0D0 * H1 (0) * (C (1) - C (0) )       &
      + H2 (N) / 3.0D0 * H1 (N - 1) * (C (N) - C (N - 1) )              
      DO 70 I = 1, N - 1, 1 
         A (I) = FN (I) - H2 (I) / 3.0D0 * (C (I - 1) * H1 (I - 1)      &
         - H3 (I - 1) * C (I) + C (I + 1) * H1 (I) )                    
   70 END DO 
      A (N) = A (0) 
!                                                                       
      DO 80 I = 0, N - 1, 1 
         B (I) = H1 (I) * (A (I + 1) - A (I) ) - H (I) / 3.0D0 *        &
         (C (I + 1) + 2.0D0 * C (I) )                                   
         D (I) = H1 (I) / 3.0D0 * (C (I + 1) - C (I) ) 
   80 END DO 
      RETURN 
      END SUBROUTINE CFSPPE                         
