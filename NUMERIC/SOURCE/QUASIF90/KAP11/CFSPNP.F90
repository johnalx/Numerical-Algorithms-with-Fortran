![KA{P 11}{Polynomial Fitting Splines}                                  
![        {Polynomial Cubic Fitting Splines for Constructing            
![         Smooth Curves}*)                                             
![  {Non--Parametric Cubic Fitting Splines}                             
![  {Non--Parametric Cubic Fitting Splines}*)                           
      SUBROUTINE CFSPNP (N, XN, FN, W, IB, ALPHA, BETA, A, B, C, D,     &
      AUXF, IERR)                                                       
!                                                                       
!*****************************************************************      
!                                                                *      
!  CFSPNP computes the coefficients A(I), B(I), C(I), D(I) for   *      
!  I=0, 1, ..., N-1 of a nonparametric cubic fitting spline.     *      
!  The end point condition is to be prescribed via the parameter *      
!  IB.                                                           *      
!  The splinefunction is represented in the form:                *      
!                                                                *      
!  S(X) = A(I) + B(I)(X-XN(I)) + C(I)(X-XN(I))**2 +              *      
!                              + D(I)(X-XN(I))**3                *      
!                                                                *      
!  for X in the interval [XN(I),XN(I+1)], I=0, 1, ..., N-1.      *      
!                                                                *      
!                                                                *      
!  ASSUMPTIONS:    1.         N > 4      , for IB = 1, 2 or 3    *      
!  ============               N > 5      , for IB =  4           *      
!                  2.     XN(I) < XN(I+1), I=0, 1, ..., N-1      *      
!                  3.      W(I) > 0.0    , I=0, 1, ..., N        *      
!                  4.      W(0) = W(N)   , for IB = 4            *      
!                  5.     FN(0) = FN(N)  , for IB = 4            *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N  :  Index of the last node                                  *      
!  XN :  vector XN(0:N); XN(I) is the Ith node, I = 0, ..., N    *      
!  FN :  vector FN(0:N); FN(I) is the data at the node XN(I)     *      
!  W  :  vector W(0:N);  W(I) is the weight of FN(I)             *      
!                                                                *      
!  IB :  determines end point condition:                         *      
!        IB = 1:  first end point derivative prescribed          *      
!        IB = 2:  second end point derivative prescribed         *      
!        IB = 3:  third end point derivative prescribed          *      
!        IB = 4:  periodic spline                                *      
!                                                                *      
!  ALPHA :  IB end point derivative at XN(0) ] for IB=1,2,3;     *      
!  BETA  :  IB end point derivative at XN(N) ] meaningless for   *      
!                                              IB=4              *      
!                                                                *      
!           (A natural fitting spline will be achieved for IB=2  *      
!            and ALPHA = BETA =0.0)                              *      
!                                                                *      
!                                                                *      
!  AUXILIARY VARIABLES:                                          *      
!  ====================                                          *      
!  AUXF :  vector AUXF(1:14*N-10)                                *      
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
!          = -1 :  N < 5  if IB = 1, 2 or 3                      *      
!                  N < 6  if IB = 4                              *      
!          = -2 :  IB < 1  or  IB > 4                            *      
!          = -3 :  Inadmissable weight W                         *      
!          = -4 :  nodes XN(I) not ordered monotonically:        *      
!                  XN(I) >= XN(I+1) for some  I=0, 1, ..., N-1   *      
!          = -5 :  IB = 4 and FN(0) not equal to FN(N) or        *      
!                  W(0) not equal to W(N)                        *      
!          =  1 :  Error in FDISY, FDIAG or NCYFSY (numerically  *      
!                  singular system matrix)                       *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  required subroutines: CFSP1D, CFSP2D, CFSP3D, CFSPPE,         *      
!                        FDISY, FDISYS, NCYFSY, NCYFSP,          *      
!                        NCYFSS, FDIAG                           *      
!                                                                *      
!                                                                *      
!  Reference: Engeln-MÅllges, G.; Reutter, F., [ENGE87].         *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : GÅnter Palm                                        *      
!  Date     : 04.18.1988                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION XN (0:N), FN (0:N), W (0:N), A (0:N), B (0:N),    &
      C (0:N), D (0:N), AUXF (1:14 * N - 10)                            
!                                                                       
!-----Checking the assumptions                                          
!                                                                       
      IERR = - 1 
      IF (N.LT.5) RETURN 
      DO 10 I = 0, N - 1, 1 
         IF (XN (I) .GE.XN (I + 1) ) THEN 
            IERR = - 4 
            RETURN 
         ENDIF 
   10 END DO 
      DO 20 I = 0, N, 1 
         IF (W (I) .LE.0.0D0) THEN 
            IERR = - 3 
            RETURN 
         ENDIF 
   20 END DO 
!                                                                       
!-----Compute the spline coefficients                                   
!                                                                       
      IF (IB.EQ.1) THEN 
         CALL CFSP1D (N, XN, FN, W, ALPHA, BETA, 1, A, B, C, D, AUXF (1)&
         , AUXF (N + 1), AUXF (2 * N + 1), AUXF (3 * N + 1), AUXF (4 *  &
         N), AUXF (5 * N - 1), AUXF (6 * N - 2), IERR)                  
      ELSEIF (IB.EQ.2) THEN 
         CALL CFSP2D (N, XN, FN, W, ALPHA, BETA, 1, A, B, C, D, AUXF (1)&
         , AUXF (N + 1), AUXF (2 * N + 1), AUXF (3 * N + 1), AUXF (4 *  &
         N), AUXF (5 * N - 1), AUXF (6 * N - 2), IERR)                  
      ELSEIF (IB.EQ.3) THEN 
         CALL CFSP3D (N, XN, FN, W, ALPHA, BETA, A, B, C, D, AUXF (1),  &
         AUXF (N), AUXF (2 * N - 1), AUXF (3 * N - 2), AUXF (4 * N - 3),&
         AUXF (5 * N - 4), AUXF (6 * N - 4), IERR)                      
      ELSEIF (IB.EQ.4) THEN 
         IF (N.LT.6) RETURN 
         CALL CFSPPE (N, XN, FN, W, 1, A, B, C, D, AUXF (1), AUXF (N +  &
         2), AUXF (2 * N + 3), AUXF (3 * N + 4), AUXF (4 * N + 5),      &
         AUXF (5 * N + 5), IERR)                                        
      ELSE 
         IERR = - 2 
      ENDIF 
      RETURN 
      END SUBROUTINE CFSPNP                         
