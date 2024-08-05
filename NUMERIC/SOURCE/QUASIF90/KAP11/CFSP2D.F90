      SUBROUTINE CFSP2D (N, XN, FN, W, ALPHA, BETA, MREP, A, B, C, D, H,&
      H1, H2, DM, DU1, DU2, RS, IERR)                                   
!                                                                       
!*****************************************************************      
!                                                                *      
!  CFSP2D computes the coefficients A(I), B(I), C(I), D(I),      *      
!  I=0, 1. ..., N-1, of a cubic fitting spline with prescribed   *      
!  second end point derivative.                                  *      
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
!  REMARK:  CFSP2D should not be called directly, but rather via *      
!  =======  the subroutine CFSPNP, or in case of parametric      *      
!           splines via the subroutine CFSPPA. The subroutines   *      
!           CFSPNP and CFSPPA also check the above assumptions.  *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N  :  Index of the last node                                  *      
!  XN :  vector XN(0:N); XN(I) is the Ith node, I = 0, ..., N    *      
!  FN :  vector FN(0:N); FN(I) is the data at the node XN(I)     *      
!  W  :  vector W(0:N);  W(I) is the weight of FN(I)             *      
!                                                                *      
!  ALPHA :  second end point derivative at XN(0)                 *      
!  BETA  :  second end point derivative at XN(N)                 *      
!                                                                *      
!           (For ALPHA = BETA = 0.0 one will obtain a natural    *      
!            fitting spline.)                                    *      
!                                                                *      
!  MREP  :  indicator used for repeated call of the subroutine:  *      
!           MREP = 1: The system matrix elements must be computed*      
!                     This matrix must be factored via subroutine*      
!                     FDISY in order to find C(I).               *      
!           MREP = 2: The right-hand side only needs to be com-  *      
!                     puted. We can use the vectors DM, DU1 and  *      
!                     DU2 computed during the first pass of sub- *      
!                     routine FDISYS to find the solution.       *      
!                     This avoids a repeat factorization in case *      
!                     of parametric splines.                     *      
!                     The elements in H, H1, H2, DM, DU1 and DU2 *      
!                     must not be altered after the first call.  *      
!  REMARK: For parametric splines with differing weights WX(I)   *      
!          not equal to WY(I) for at least one index I, one      *      
!          must work with MREP = 1.                              *      
!                                                                *      
!                                                                *      
!  AUXILIARY VARIABLES:                                          *      
!  ====================                                          *      
!  H   :]                                                        *      
!  H1  :]  N-vectors H(0:N-1), H1(0:N-1), H2(0:N-1)              *      
!  H2  :]                                                        *      
!                                                                *      
!  DM  :]                                                        *      
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
!          =  1 :  FDISY did not run correctly (matrix singular) *      
!          = -6 :  Wrong input for MREP                          *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  required subroutines: FDISY, FDISYS                           *      
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
      C (0:N), D (0:N), H (0:N - 1), H1 (0:N - 1), H2 (0:N - 1),        &
      DM (1:N - 1), DU1 (1:N - 1), DU2 (1:N - 1), RS (1:N - 1)          
!                                                                       
!-----Check MREP for repeated calls                                     
!                                                                       
      IERR = - 6 
      IF (MREP.NE.1.AND.MREP.NE.2) RETURN 
!                                                                       
!-----Compute auxiliary values and system matrix elements               
!     i.e., its main and two co-diagonals, if in the first pass         
!                                                                       
      IF (MREP.EQ.1) THEN 
!                                                                       
!       Auxiliary variables                                             
!                                                                       
         DO 10 I = 0, N - 1, 1 
            H (I) = XN (I + 1) - XN (I) 
            H1 (I) = 1.0D0 / H (I) 
            C (I) = H1 (I) * H1 (I) 
            B (I) = 6.0D0 / W (I) 
   10    END DO 
         B (N) = 6.0D0 / W (N) 
!                                                                       
         DO 20 I = 0, N - 2, 1 
            H2 (I) = H1 (I) + H1 (I + 1) 
   20    END DO 
!                                                                       
!       second co-diagonal                                              
!                                                                       
         DO 30 I = 1, N - 3, 1 
            DU2 (I) = B (I + 1) * H1 (I) * H1 (I + 1) 
   30    END DO 
!                                                                       
!       first co-diagonal                                               
!                                                                       
         DO 40 I = 1, N - 2, 1 
            DU1 (I) = H (I) - B (I) * H1 (I) * H2 (I - 1) - B (I + 1)   &
            * H1 (I) * H2 (I)                                           
   40    END DO 
!                                                                       
!       main diagonal                                                   
!                                                                       
         DO 50 I = 1, N - 1, 1 
            K = I - 1 
            DM (I) = 2.0D0 * (H (K) + H (I) ) + B (K) * C (K) + B (I)   &
            * H2 (K) * H2 (K) + B (I + 1) * C (I)                       
   50    END DO 
      ENDIF 
!                                                                       
!-----Compute the right-hand side                                       
!                                                                       
      C (0) = 0.5D0 * ALPHA 
      C (N) = 0.5D0 * BETA 
!                                                                       
      DUMMY2 = (FN (2) - FN (1) ) * H1 (1) 
      DUMMY1 = (FN (3) - FN (2) ) * H1 (2) 
      RS (1) = 3.0D0 * (DUMMY2 - (FN (1) - FN (0) ) * H1 (0) ) - C (0)  &
      * (H (0) - 6.0D0 / W (0) * H1 (0) * H1 (0) - 6.0D0 / W (1)        &
      * H1 (0) * H2 (0) )                                               
      RS (2) = 3.0D0 * (DUMMY1 - DUMMY2) - C (0) * (6.0D0 / W (1) )     &
      * H1 (0) * H1 (1)                                                 
      DO 60 I = 3, N - 3, 1 
         DUMMY2 = (FN (I + 1) - FN (I) ) * H1 (I) 
         RS (I) = 3.0D0 * (DUMMY2 - DUMMY1) 
         DUMMY1 = DUMMY2 
   60 END DO 
      DUMMY2 = (FN (N - 1) - FN (N - 2) ) * H1 (N - 2) 
      RS (N - 2) = 3.0D0 * (DUMMY2 - DUMMY1) - C (N) * (6.0D0 / W (N -  &
      1) ) * H1 (N - 2) * H1 (N - 1)                                    
      RS (N - 1) = 3.0D0 * ( (FN (N) - FN (N - 1) ) * H1 (N - 1)        &
      - DUMMY2) - C (N) * (H (N - 1) - 6.0D0 / W (N - 1) * H1 (N - 1)   &
      * H2 (N - 2) - 6.0D0 / W (N) * C (N) )                            
!                                                                       
!-----Compute the coefficients C(1) to C(N-1) from                      
!     the system of equations                                           
!                                                                       
      IF (MREP.EQ.1) THEN 
!                                                                       
!       In case we must decompose the system matrix                     
!                                                                       
         CALL FDISY (N - 1, DM, DU1, DU2, RS, C (1), IFLAG) 
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
!       When no factorization is necessary                              
!                                                                       
         CALL FDISYS (N - 1, DM, DU1, DU2, RS, C (1) ) 
      ENDIF 
      IERR = 0 
!                                                                       
!-----Compute the remaining spline coefficients                         
!                                                                       
      A (0) = FN (0) + 2.0D0 / W (0) * H1 (0) * (C (0) - C (1) ) 
      DO 70 I = 1, N - 1, 1 
         A (I) = FN (I) - 2.0D0 / W (I) * (C (I - 1) * H1 (I - 1)       &
         - H2 (I - 1) * C (I) + C (I + 1) * H1 (I) )                    
   70 END DO 
      A (N) = FN (N) - 2.0D0 / W (N) * H1 (N - 1) * (C (N - 1) - C (N) ) 
!                                                                       
      DO 80 I = 0, N - 1, 1 
         B (I) = H1 (I) * (A (I + 1) - A (I) ) - H (I) / 3.0D0 *        &
         (C (I + 1) + 2.0D0 * C (I) )                                   
         D (I) = H1 (I) / 3.0D0 * (C (I + 1) - C (I) ) 
   80 END DO 
      RETURN 
      END SUBROUTINE CFSP2D                         
