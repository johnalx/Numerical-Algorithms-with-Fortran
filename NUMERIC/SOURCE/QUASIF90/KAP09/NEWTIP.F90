![KA{P 9}{Polynomial and Rational Interpolation}                        
![  {Polynomial and Rational Interpolation}*)                           
![  {Newton Formula for Arbitrary Nodes}                                
![  {Newton Formula for Arbitrary Nodes}*)                              
      SUBROUTINE NEWTIP (N, X, Y, B, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     NEWTIP determines the coefficients of the interpolation    *      
!     polynomial in Newton's form.                               *      
!     Functional values of this polynomial can then be determined*      
!     by FUNCTION NIPFCT.                                        *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     N       : degree of the interpolation polynomial           *      
!               (= index of the last node, if counted from 0 on) *      
!     X       : (N+1)-vector X(0:N) ] value pairs (X(I),Y(I)) to *      
!     Y       : (N+1)-vector Y(0:N) ] be interpolated by a       *      
!                                   ] polynomial                 *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     B       : (N+1)-vector B(0:N), the coefficients of the     *      
!               polynomial in the following form:                *      
!                  NP(X) = B(0) + B(1)*(X-X(0)) +                *      
!                               + B(2)*(X-X(0))*(X-X(1)) + ...   *      
!                           ... + B(N)*(X-X(0))*...*(X-X(N-1))   *      
!     IERR    : error code                                       *      
!               = 0 : no error                                   *      
!               = 1 : N < 0                                      *      
!               = 2 : there are two identical nodes              *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Elmar Pohl                                         *      
!  date     : 09.28.1985                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION X (0:N), Y (0:N), B (0:N) 
      IERR = 0 
      IF (N.LT.0) THEN 
         IERR = 1 
         RETURN 
      ENDIF 
      DO 10 I = 0, N 
         B (I) = Y (I) 
   10 END DO 
      DO 20 I = 1, N 
         DO 30 K = N, I, - 1 
            H = X (K) - X (K - I) 
            IF (H.EQ.0.0D0) THEN 
               IERR = 2 
               RETURN 
            ENDIF 
            B (K) = (B (K) - B (K - 1) ) / H 
   30    END DO 
   20 END DO 
      RETURN 
      END SUBROUTINE NEWTIP                         
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION NIPFCT (X0, X, B, N) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     NIPFCT determines the functional value of the Newton       *      
!     interpolation polynomial with the coefficients B(I) for    *      
!     I=0,1, ...,N at X0 by a generalized Horner scheme.         *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     X0      : X-value where the polynomial is to be evaluated  *      
!     X       : (N+1)-vector X(0:N) of nodes X(I), I=0,1,...,N,  *      
!               for the interpolating polynomial                 *      
!     B       : (N+1)-vector B(0:N), the coefficients of the     *      
!               Newton interpolation polynomial, in the following*      
!               form:                                                   
!                  NP(X) = B(0) + B(1)*(X-X(0)) +                *      
!                               + B(2)*(X-X(0))*(X-X(1)) + ...   *      
!                           ... + B(N)*(X-X(0))*...*(X-X(N-1))   *      
!     N       : degree of the interpolating polynomial           *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETER:                                          *      
!     =================                                          *      
!     NIPFCT  : NP(X0), the value of the interpolating polynomial*      
!               at X0                                            *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Elmar Pohl                                         *      
!  date     : 09.28.1985                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION B (0:N), X (0:N) 
      NIPFCT = B (N) 
      DO 10 I = N - 1, 0, - 1 
         NIPFCT = NIPFCT * (X0 - X (I) ) + B (I) 
   10 END DO 
      RETURN 
      END FUNCTION NIPFCT                           
