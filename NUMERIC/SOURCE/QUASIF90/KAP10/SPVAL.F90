      DOUBLEPRECISION FUNCTION SPVAL (X, N, XN, A, B, C, D) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  The FUNCTION SPVAL determines the functional value of a cubic *      
!  spline S of the form:                                         *      
!                                                                *      
!  S(X) = P(I)(X) = A(I) + B(I)*(X-XN(I)) + C(I)*(X-XN(I))**2 +  *      
!                                         + D(I)*(X-XN(I))**3    *      
!                                                                *      
!  for X in the interval [XN(I),XN(I+1)], where I=0,1,...,N-1.   *      
!                                                                *      
!  If X < XN(0) we evaluate the first polynom P(0), if X > XN(N) *      
!  we evaluate the polynom P(N-1).                               *      
!  We do not check the input X in any way.                       *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  X  :  value for which we want to find the value of the spline *      
!  N  :  index of the final node                                 *      
!  XN :  vector of nodes XN(0:N)                                 *      
!  A  :  ] N+1-vectors ..(0:N);                                  *      
!  B  :  ] the spline coefficients A(I), B(I), C(I), D(I) for    *      
!  C  :  ] I = 0, ..., N-1; used for storage for I = N.          *      
!  D  :  ]                                                       *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : GÅnter Palm                                        *      
!  Date     : 01.06.1991                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
!  Declarations                                                         
!                                                                       
      DOUBLEPRECISION XN (0:N), A (0:N), B (0:N), C (0:N), D (0:N) 
      SAVE I 
!                                                                       
!  Initializing                                                         
!                                                                       
      DATA I / 0 / 
      IF (I.GE.N) I = 0 
!                                                                       
!  In a repeat call we determine the interval [XN(I), XN(I+1)]          
!  that contains X only if X does not lie in the same interval          
!  as last used.                                                        
!                                                                       
      IF (X.LT.XN (I) .OR.X.GE.XN (I + 1) ) THEN 
         I = 0 
         K = N 
   10    M = (I + K) / 2 
         IF (X.LT.XN (M) ) THEN 
            K = M 
         ELSE 
            I = M 
         ENDIF 
         IF (K.GT.I + 1) GOTO 10 
      ENDIF 
!                                                                       
!  Compute X-XN(I) in order to evaluate the polynomial                  
!                                                                       
      XL = X - XN (I) 
!                                                                       
!  Compute the value of the spline via a Horner scheme                  
!                                                                       
      SPVAL = ( (D (I) * XL + C (I) ) * XL + B (I) ) * XL + A (I) 
      RETURN 
      END FUNCTION SPVAL                            
