      SUBROUTINE RALGIP (X, A, N, L, L1, IERR, NLEFT, Z) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Rational Lagrange interpolation. The interpolation         *      
!     function defined by the nodes X1, X2, X3, ...  and the     *      
!     functional values is the ratio of two polynomials. The     *      
!     degree of the nominator polynomial can be chosen by the    *      
!     user.                                                      *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     X      : (N+1)-vector X(0:N); the nodes X(I) for I=0,...,N *      
!     A      : (N+1)-vector A(0:N); the functional values at     *      
!              the nodes  X(I)                                   *      
!     N      : largest index used for the nodes                  *      
!     L      : degree of the denominator polynomial; L < N       *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     X      : as above, but possibly reordered                  *      
!     A      : (N+1)-vector A(0:N); the coefficients of the      *      
!              rational interpolation function                   *      
!     L1     : (N+1)-vector L1(0:N); contains information about  *      
!              the coefficients. This is needed for use in the   *      
!              FUNCTION RLIFCT.                                  *      
!     IERR   : = 0  everything o.k.                              *      
!              = 1  problem not solvable (check the nodes)       *      
!     NLEFT  : final index of the nodes for which we still have  *      
!              to interpolate                                    *      
!                                                                *      
!                                                                *      
!     AUXILIARY VARIABLES:                                       *      
!     ====================                                       *      
!     Z      : N-vector Z(0:N-1)                                 *      
!                                                                *      
!     The interpolating function has the continued fraction form:*      
!                                                                *      
!     F(X) =  A(N)  +                                            *      
!                                                                *      
!             + (X-X(N)) * A(N-1) +                              *      
!                                                                *      
!                          (X-X(N)) * (X-X(N-1))                 *      
!             + ------------------------------------------------ *      
!                                          (X-X(N-2))*(X-X(N-3)) *      
!               A(N-2) + (X-X(N-2))*A(N-3)+--------------------- *      
!                                          A(N-4)  +   . . .     *      
!                                                                *      
!     F(X) can be evaluated using the FUNCTION RLIFCT.           *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: RLIFCT, RLIORD                          *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Helmut Werner                                      *      
!  editor   : Christiane Beer                                    *      
!  date     : 1983                                               *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION X (0:N), A (0:N), Z (0:N - 1) 
      INTEGER L1 (0:N) 
      IERR = 0 
      EPS = 1.0D-11 
!                                                                       
!     determine the degree of the denominator polynomial                
!                                                                       
      L2 = L 
      M1 = N - L 
      NLEFT = N 
      IF (M1.LE.0) THEN 
         IERR = 1 
         RETURN 
      ENDIF 
!                                                                       
!     the degree of the numerator polynomial exceeds that of the denomin
!     a number of nodes equal to the difference between the denominator 
!     the numerator degree is interpolated by forming divided difference
!                                                                       
  100 DO 10 K = 1, (L2 - M1) 
         CALL RLIORD (X, A, NLEFT, XJ, AJ) 
         DO 20 I = 0, NLEFT - 1 
!                                                                       
!   stop if two nodes are identical                                     
!                                                                       
            IF ( (X (I) - XJ) .EQ.0.0D0) THEN 
               IERR = 1 
               RETURN 
            ENDIF 
            A (I) = (A (I) - AJ) / (X (I) - XJ) 
   20    END DO 
         L1 (NLEFT - 1) = 1 
         NLEFT = NLEFT - 1 
   10 END DO 
      IF (NLEFT.LT.0) GOTO 200 
!                                                                       
!   the subsequent node is interpolated by forming                      
!   inverse divided differences                                         
!                                                                       
      CALL RLIORD (X, A, NLEFT, XJ, AJ) 
      I1 = 0 
      DO 30 I = 0, NLEFT - 1 
         A2 = A (I) - AJ 
         X2 = X (I) - XJ 
!                                                                       
!   test for automatically interpolated nodes                           
!                                                                       
         IF (DABS (A2) .LE.DABS (X2) * EPS) THEN 
            Z (I1) = X (I) 
            I1 = I1 + 1 
         ELSE 
            A (I - I1) = X2 / A2 
            X (I - I1) = X (I) 
         ENDIF 
   30 END DO 
!                                                                       
!   automatically interpolated nodes are included in the                
!   divided differences                                                 
!                                                                       
      DO 40 K = 0, I1 - 1 
         X (NLEFT - 1) = Z (K) 
         A (NLEFT - 1) = 0.0D0 
         DO 50 I = 0, NLEFT - 1 
            A (I) = A (I) * (X (I) - X (NLEFT) ) 
   50    END DO 
         L1 (NLEFT - 1) = 1 
         NLEFT = NLEFT - 1 
   40 END DO 
!                                                                       
!   determine the degrees of the factors of the interpolating           
!   rational function that still must be determined                     
!                                                                       
      IF (NLEFT.GT.0) THEN 
         NLEFT = NLEFT - 1 
         L1 (NLEFT) = - 1 
         L2 = M1 
         M1 = NLEFT - L2 
      ENDIF 
  200 IF (M1.LT.0) THEN 
         IERR = 1 
         RETURN 
      ENDIF 
!                                                                       
!   check whether interpolation is complete                             
!                                                                       
      IF (NLEFT.GT.0) GOTO 100 
!                                                                       
!   finish interpolation                                                
!   The following tests whether all nodes have been                     
!   used in the interpolation                                           
!                                                                       
      FM = DABS (A (N) ) 
      DO 60 I = 0, N - 1 
         FM = FM + DABS (A (I) ) 
   60 END DO 
      K = 0 
      DO 70 I = 0, N - 1 
         IF (L1 (I) .LT.0) K = I + 1 
         IF (K.NE.0) THEN 
            F = RLIFCT (X (I + 1), X, A, L1, I) 
            IF (DABS (F) .LE.FM * EPS) THEN 
               IERR = 1 
               RETURN 
            ENDIF 
         ENDIF 
   70 END DO 
      RETURN 
      END SUBROUTINE RALGIP                         
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION RLIFCT (XW, X, A, L1, N) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Evaluate the rational interpolation function               *      
!     (see SUBROUTINE RALGIP)                                    *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     XW     : point where the rational interpolation function   *      
!              is to be evaluated                                *      
!     X      : (N+1)-vector X(0:N); the nodes, I=0,...,N         *      
!     A      : (N+1)-vector A(0:N); the coefficients of the      *      
!              rational interpolating function                   *      
!     L1     : (N+1)-vector L1(0:N); labels indicating whether   *      
!              we need to divide or mulptiply after evaluation   *      
!              of a Horner scheme                                *      
!     N      : final index of the nodes                          *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETER:                                          *      
!     =================                                          *      
!     RLIFCT : functional value at XW                            *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Helmut Werner                                      *      
!  editor   : Christiane Beer                                    *      
!  date     : 1983                                               *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION X (0:N), A (0:N) 
      INTEGER L1 (0:N) 
      F = A (0) 
      DO 10 I = 1, N 
         IF (L1 (I - 1) .GE.0) THEN 
            F = A (I) + (XW - X (I) ) * F 
         ELSE 
            F = A (I) + (XW - X (I) ) / F 
         ENDIF 
   10 END DO 
      RLIFCT = F 
      RETURN 
      END FUNCTION RLIFCT                           
!                                                                       
!                                                                       
      SUBROUTINE RLIORD (X, A, N, XJ, AJ) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Determines the functional value of least magnitude         *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     X      : (N+1)-vector X(0:N); the nodes X(I), I=0,...,N    *      
!     A      : (N+1)-vector A(0:N); the functional values at X(I)*      
!     N      : final index of the nodes                          *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     X      : as above, permutation of the node to be           *      
!              interpolated next with the last one used          *      
!     A      : see above, swap A(I) analogous to the X(I)        *      
!     XJ     : node with the functional value of smallest        *      
!              magnitude                                         *      
!     AJ     : largest absolute functional value                 *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Helmut Werner                                      *      
!  editor   : Christiane Beer                                    *      
!  date     : 1983                                               *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION X (0:N), A (0:N) 
      AJ = A (N) 
      J1 = N 
      DO 10 I = 0, N - 1 
         IF (DABS (AJ) .GT.DABS (A (I) ) ) THEN 
            J1 = I 
            AJ = A (I) 
         ENDIF 
   10 END DO 
      XJ = X (J1) 
      X (J1) = X (N) 
      A (J1) = A (N) 
      X (N) = XJ 
      A (N) = AJ 
      RETURN 
      END SUBROUTINE RLIORD                         
