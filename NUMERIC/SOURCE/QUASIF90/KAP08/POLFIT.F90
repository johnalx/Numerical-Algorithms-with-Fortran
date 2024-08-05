![  {Discrete Least Squares using Orthogonal Polynomials}               
![  {Discrete Least Squares via Algebraic Polynomials                   
![   using Orthogonal Polynomials}*)                                    
      SUBROUTINE POLFIT (N, X, F, W, NDEG, ALPHA, B, C, SUMERS, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     This program performs a discrete polynomial approximation  *      
!     using orthogonal polynomials:                              *      
!     For N+1 given tuples ( X(J), F(J) ), J = 0, 1, ..., N,     *      
!     the polynomial coefficients ALPHA(0), ALPHA(1), .. ,       *      
!     ALPHA(NDEG) are computed so that the weighted least squares*      
!     error                                                      *      
!                                                                *      
!      SUMERS =                                                  *      
!     ( SUM( J=0 to N ) W(J)*( F(J) - POLY(X(J)) )**2 )** (1/2)  *      
!                                                                *      
!     becomes minimal for the positive weights W(0), W(1), ..,   *      
!     W(N). Here                                                 *      
!                                                                *      
!       POLY(X) = SUM( K=0 to NDEG ) ALPHA(K) * Q(K,X)           *      
!                                                                *      
!     denotes the approximating polynomial of degree NDEG.       *      
!     The discrete orthogonal polynomials Q(K,X) of degree       *      
!     K satisfy the recursion                                    *      
!                                                                *      
!       Q(K,X) = ( X - B(K) ) * Q(K-1,X) - C(K) * Q(K-2,X),      *      
!       for K >= 2 with  Q(0,X) := 1  and  Q(1,X) := X - B(1).   *      
!                                                                *      
!     By using the values for B(1), B(2), .. , B(NDEG) and       *      
!     C(2), .., C(NDEG), the approximating polynomial can be     *      
!     evaluated efficiently using the function POEVAL.           *      
!     ( the standard FORTRAN conventions apply to the variable   *      
!       types used for the transfer parameters).                 *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!        N       : the number of tuples used is N + 1            *      
!                  (compare with the dimension of the vectors    *      
!                   X, F, W)                                     *      
!      X(0:N)    : nodes used                                    *      
!      F(0:N)    : functional values at the nodes                *      
!      W(0:N)    : weights with W(J) > 0 for all J; for equal    *      
!                  weighting set W(J) = 1 for all J.             *      
!      NDEG      : degree of the approximating polynomial,       *      
!                  0 <= NDEG <= N.                               *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  ALPHA(0:NDEG) : coefficients of the approximating polynomial  *      
!                  with respect to the orthogonal polynomials    *      
!                  Q(K,X)                                        *      
!  B(1:NDEG)   )   coefficients for the recursive calculation    *      
!  C(2:NDEG)   )   of the polynomials Q(K,X)                     *      
!                  ( if NDEG < 2 or NDEG < 1 the coefficients    *      
!                    in C or in B and C are not needed).         *      
!  SUMERS        : weighed least squares error of the computed   *      
!                  approximating polynomial at the given nodes   *      
!  IERR          : error parameter                               *      
!                  = 0 : everything o.k.                         *      
!                  = 1 : the inequalities 0 <= NDEG <= N are not *      
!                        true                                    *      
!                  = 2 : not all weights are positive            *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: SCALP, Q, POEVAL                        *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Klaus Niederdrenk                                  *      
!  date     : 05.27.87                                           *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION X (0:N), F (0:N), W (0:N) 
      DIMENSION ALPHA (0:NDEG), B (1:NDEG), C (2:NDEG) 
!                                                                       
! ***  checking the input parameters                                    
!                                                                       
      IERR = 0 
      IF (N.LT.NDEG.OR.NDEG.LT.0) THEN 
         IERR = 1 
      ELSE 
         DO 10 J = 0, N 
            IF (W (J) .LE.0.0D0) IERR = 2 
   10    END DO 
      ENDIF 
      IF (IERR.NE.0) RETURN 
!                                                                       
! *** calculating the coefficients ALPHA(K), B(K) and C(K)              
! *** by applying the function SCALP                                    
!                                                                       
      QKMIN1 = SCALP (0, N, X, F, W, B, C, '(QK,QK)') 
      B (1) = SCALP (0, N, X, F, W, B, C, '(XQK,QK)') / QKMIN1 
      ALPHA (0) = SCALP (0, N, X, F, W, B, C, '(F,QK)') / QKMIN1 
      DO 20 K = 2, NDEG 
         QKMIN2 = QKMIN1 
         QKMIN1 = SCALP (K - 1, N, X, F, W, B, C, '(QK,QK)') 
         B (K) = SCALP (K - 1, N, X, F, W, B, C, '(XQK,QK)') / QKMIN1 
         C (K) = QKMIN1 / QKMIN2 
         ALPHA (K - 1) = SCALP (K - 1, N, X, F, W, B, C, '(F,QK)')      &
         / QKMIN1                                                       
   20 END DO 
      DUMMY = SCALP (NDEG, N, X, F, W, B, C, '(F,QK)') 
      ALPHA (NDEG) = DUMMY / SCALP (NDEG, N, X, F, W, B, C, '(QK,QK)') 
!                                                                       
! *** determine the weighted least squares error                        
!                                                                       
      SUMERS = 0.0D0 
      DO 30 J = 0, N 
         DUMMY = F (J) - POEVAL (X (J), NDEG, ALPHA, B, C) 
         SUMERS = SUMERS + W (J) * DUMMY * DUMMY 
   30 END DO 
      SUMERS = DSQRT (SUMERS) 
!                                                                       
      RETURN 
      END SUBROUTINE POLFIT                         
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION SCALP (K, N, X, F, W, B, C, CHOICE) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     This program uses the formula                              *      
!                                                                *      
!       ( G1 , G2 ) = SUMME( J=0 to N ) W(J)*G1(X(J))*G2(X(J))   *      
!                                                                *      
!     to compute the weighted scalar products required for the   *      
!     discrete polynomial approximation.                         *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: Q                                       *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Klaus Niederdrenk                                  *      
!  date     : 05.27.87                                           *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION X (0:N), F (0:N), W (0:N) 
      DIMENSION B (1: * ), C (2: * ) 
      CHARACTER ( * ) CHOICE 
!                                                                       
      SCALP = 0.0D0 
      IF (CHOICE.EQ.'(QK,QK)') THEN 
         DO 10 J = 0, N 
            DUMMY = Q (K, X (J), B, C) 
            SCALP = SCALP + W (J) * DUMMY * DUMMY 
   10    END DO 
      ELSEIF (CHOICE.EQ.'(XQK,QK)') THEN 
         DO 20 J = 0, N 
            DUMMY = Q (K, X (J), B, C) 
            SCALP = SCALP + W (J) * X (J) * DUMMY * DUMMY 
   20    END DO 
      ELSEIF (CHOICE.EQ.'(F,QK)') THEN 
         DO 30 J = 0, N 
            SCALP = SCALP + W (J) * F (J) * Q (K, X (J), B, C) 
   30    END DO 
      ENDIF 
!                                                                       
      RETURN 
      END FUNCTION SCALP                            
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION Q (K, X, B, C) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     This program evaluates the polynomial Q(K,X) of degree K   *      
!     at the given node X by applying the two-stage recursion    *      
!                                                                *      
!       Q(K,X) = ( X - B(K) ) * Q(K-1,X) - C(K) * Q(K-2,X),      *      
!       for K >= 2 with Q(0,X) = 1 and Q(1,X) = X - B(1).        *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Klaus Niederdrenk                                  *      
!  date     : 05.27.87                                           *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION B (1: * ), C (2: * ) 
!                                                                       
      IF (K.EQ.0) THEN 
         Q = 1.0D0 
      ELSEIF (K.EQ.1) THEN 
         Q = X - B (1) 
      ELSE 
         QMIN2 = 1.0D0 
         QMIN1 = X - B (1) 
         DO 10 I = 2, K 
            Q = (X - B (I) ) * QMIN1 - C (I) * QMIN2 
            QMIN2 = QMIN1 
            QMIN1 = Q 
   10    END DO 
      ENDIF 
!                                                                       
      RETURN 
      END FUNCTION Q                                
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION POEVAL (X, NDEG, ALPHA, B, C) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     In a Horner like way, this program evaluates the           *      
!     approximating polynomial                                   *      
!                                                                *      
!       POLY(X) = SUMME( K=0 to NDEG ) ALPHA(K) * Q(K,X)         *      
!                                                                *      
!     of degree NDEG, that was determined by SUBROUTINE POLFIT,  *      
!     at the given node X.                                       *      
!     ( the standard FORTRAN conventions apply to the variable   *      
!       types for the transfer parameters).                      *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!        X       : value at which the approximating polynomial   *      
!                  is to be evaluated                            *      
!      NDEG      : degree of the approximating polynomial        *      
!   ALPHA(0:NDEG): the coefficients of the approximating         *      
!                  polynomial determined by SUBROUTINE POLFIT    *      
!                  with respect to the orthogonal polynomials    *      
!                  Q(K,X)                                        *      
!   B(1:NDEG)    : ) coefficients for the recursion              *      
!   C(2:NDEG)    : )                                             *      
!                                                                *      
!  OUTPUT PARAMETER:                                             *      
!  =================                                             *      
!     POEVAL     : value of the approximating polynomial at X    *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Klaus Niederdrenk                                  *      
!  date     : 05.27.87                                           *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION ALPHA (0:NDEG), B (1:NDEG), C (2:NDEG) 
      IF (NDEG.EQ.0) THEN 
         SK = ALPHA (0) 
      ELSEIF (NDEG.EQ.1) THEN 
         SK = ALPHA (0) + ALPHA (1) * (X - B (1) ) 
      ELSE 
         SKP2 = ALPHA (NDEG) 
         SKP1 = ALPHA (NDEG - 1) + SKP2 * (X - B (NDEG) ) 
         DO 10 K = NDEG - 2, 0, - 1 
            SK = ALPHA (K) + SKP1 * (X - B (K + 1) ) - SKP2 * C (K + 2) 
            SKP2 = SKP1 
            SKP1 = SK 
   10    END DO 
      ENDIF 
      POEVAL = SK 
!                                                                       
      RETURN 
      END FUNCTION POEVAL                           
