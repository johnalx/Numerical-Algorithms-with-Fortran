![KA{P 14}{Numerical Differentiation}                                   
![        {Numerical Differentiation}*)                                 
![          {Differentiation by Using Interpolating                     
![           Polynomials}*)                                             
      SUBROUTINE SPLFVD (X, N, XN, A, B, C, D, S, S1, S2, S3) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Subroutine SPLFVD determines the functional value as well as  *      
!  the 1st, 2nd and 3rd derivatives of a cubic spline function   *      
!  S given in the following form                                 *      
!                                                                *      
!  S(X) = P(I)(X) = A(I) + B(I)*(X-XN(I)) + C(I)*(X-XN(I))**2 +  *      
!                        + D(I)*(X-XN(I))**3,                    *      
!                                                                *      
!  where X lies in the interval [XN(I), XN(I+1)], I=0, ..., N-1. *      
!                                                                *      
!  For X < XN(0) the boundary polynomial P(0) is evaluated,      *      
!  for X > XN(N) we use the boundary polynomial P(N-1).          *      
!  No plausibility check of the input values is performed.       *      
!  This program may be used in connection with other spline      *      
!  routines.                                                     *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETER:                                              *      
!  ================                                              *      
!  X  : value where S(X), S'(X), S''(X), S'''(X) are to be       *      
!       determined.                                              *      
!  N  : index of the last node.                                  *      
!  XN : (N+1)-vector XN(0:N) containing the nodes XN(I),         *      
!       I=0, ..., N.                                             *      
!  A  : ] (N+1)-vectors A(0:N), ..., D(0:N), containing the      *      
!  B  : ] spline coefficients A(I), B(I), C(I), D(I) for         *      
!  C  : ] I=0, ..., N-1 and an auxiliary variable in position    *      
!  D  : ] I=N.                                                   *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  S  : functional value S(X) of the spline function S at X.     *      
!  S1 : 1st derivative S'(X) of the spline function S at X.      *      
!  S2 : 2nd derivative S''(X) of the spline function S at X.     *      
!  S3 : 3rd derivative S'''(X) of the spline function S at X.    *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Gisela Engeln-Muellges                             *      
!  date     : 04.12.1988                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
!  declarations                                                         
!                                                                       
      DOUBLEPRECISION XN (0:N), A (0:N), B (0:N), C (0:N), D (0:N) 
      SAVE I 
!                                                                       
!  initializing                                                         
!                                                                       
      DATA I / 0 / 
      IF (I.GE.N) I = 0 
!                                                                       
!  If this is a repeated call of the subroutine, the                    
!  loop for determining the interval [XN(I), XN(I+1)],                  
!  that contains X, is only executed if X is not                        
!  inside the same interval as during the previous call.                
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
!  determine the linear factor X-XN(I) needed for the polynomial evaluat
!                                                                       
      XL = X - XN (I) 
!                                                                       
!  determine the functional and derivative values via a HORNER          
!  scheme using certain auxiliary variables                             
!                                                                       
      DUMMY1 = 3.0D0 * D (I) 
      DUMMY2 = 2.0D0 * C (I) 
      DUMMY3 = 2.0D0 * DUMMY1 
      S = ( (D (I) * XL + C (I) ) * XL + B (I) ) * XL + A (I) 
      S1 = (DUMMY1 * XL + DUMMY2) * XL + B (I) 
      S2 = DUMMY3 * XL + DUMMY2 
      S3 = DUMMY3 
      RETURN 
      END SUBROUTINE SPLFVD                         
