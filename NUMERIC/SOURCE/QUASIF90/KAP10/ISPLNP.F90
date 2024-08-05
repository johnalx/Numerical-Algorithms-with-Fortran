![KA{P 10}{Interpolating Polynomial Splines}                            
![        {Interpolating Polynomial Splines for Constructing            
![         Smooth Curves}*)                                             
![  {Non--Parametric Cubic Splines}                                     
![  {Non--Parametric Cubic Splines}*)                                   
      SUBROUTINE ISPLNP (N, XN, FN, IB, ALPHA, BETA, B, C, D, DUMMY,    &
      IERR)                                                             
!                                                                       
!*********************************************************************  
!                                                                    *  
!  ISPLNP computes the coefficients B(I), C(I), D(I) for I=0,1,...,  *  
!  N-1 of a nonparametric cubic interpolating spline for various end *  
!  point conditions.                                                 *  
!  The end point condition can be specified using the parameter IB   *  
!  The spline function S has the form:                               *  
!                                                                    *  
!  S(X) = FN(I) + B(I)(X-XN(I)) + C(I)(X-XN(I))**2 +                 *  
!                               + D(I)(X-XN(I))**3                   *  
!                                                                    *  
!  for X in the interval [XN(I),XN(I+1)], I=0,1,...,N-1.             *  
!                                                                    *  
!                                                                    *  
!  ASSUMPTIONS:    1.         N > 2                                  *  
!  ============    2.     XN(I) < XN(I+1), I=0,1,...,N-1             *  
!                  3.     FN(0) = FN(N)  , if IB = 4                 *  
!                                                                    *  
!                                                                    *  
!  INPUT PARAMETERS:                                                 *  
!  =================                                                 *  
!  N  :  index of the final node                                     *  
!  XN :  vector XN(0:N); the nodes XN(I), I = 0,1,...,N              *  
!  FN :  vector FN(0:N); the functional values FN(I) = FN( XN(I) )   *  
!                                                                    *  
!  IB :  specifies the end point condition:                          *  
!        IB = 1:  first end point derivative                         *  
!        IB = 2:  second end point derivative                        *  
!        IB = 3:  third end point derivative                         *  
!        IB = 4:  periodic spline                                    *  
!        IB = 5:  'not-a-node' - condition                           *  
!                                                                    *  
!  ALPHA :  IBth derivative at XN(0) ]  used only for IB = 1,2,3;    *  
!  BETA  :  IBth derivative at XN(N) ]  not used for IB = 4,5.       *  
!                                                                    *  
!  (A natural spline will result by setting ALPHA = BETA = 0.0 and   *  
!   IB = 2.)                                                         *  
!                                                                    *  
!                                                                    *  
!  AUXILIARY VARIABLES:                                              *  
!  ====================                                              *  
!  DUMMY :  vector DUMMY(1:5*N+1)                                    *  
!                                                                    *  
!                                                                    *  
!  OUTPUT PARAMETERS:                                                *  
!  ==================                                                *  
!  FN :  ]  N+1-vectors ..(0:N);                                     *  
!  B  :  ]  The first N entries of B, C and D are the spline         *  
!  C  :  ]  coefficients for S. B(N), C(N), D(N) are auxiliary       *  
!  D  :  ]  variables.                                               *  
!  IERR :  error parameter                                           *  
!          =  0 :  All is ok                                         *  
!          = -1 :  N < 4                                             *  
!          = -2 :  IB < 1 or IB > 5                                  *  
!          = -3 :  nodes not monotone, XN(I) => XN(I+1) for some     *  
!                  I = 0, 1, ..., N-1                                *  
!          = -4 :  IB = 4 and FN(0) not equal to FN(N)               *  
!          =  1 :  crash in SUBROUTINE TRDSY, TRDIG or CYTSY,        *  
!                  system matrix is numerically singular             *  
!                                                                    *  
!--------------------------------------------------------------------*  
!                                                                    *  
!  Subroutines required: ISPL1D, ISPL2D, ISPL3D, ISPLPE, ISPLNK,     *  
!                        TRDSY, TRDSYS, CYTSY, CYTSYS, TRDIG         *  
!                                                                    *  
!*********************************************************************  
!                                                                    *  
!  Author   : GÅnter Palm                                            *  
!  Date     : 15.04.1988                                             *  
!  Source   : FORTRAN 77                                             *  
!                                                                    *  
!*********************************************************************  
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION XN (0:N), FN (0:N), B (0:N), C (0:N), D (0:N),    &
      DUMMY (1:5 * N + 1)                                               
!                                                                       
!-----Check input data                                                  
!                                                                       
      IERR = - 1 
      IF (N.LT.3) RETURN 
      DO 10 I = 0, N - 1, 1 
         IF (XN (I) .GE.XN (I + 1) ) THEN 
            IERR = - 3 
            RETURN 
         ENDIF 
   10 END DO 
!                                                                       
!-----Compute spline coefficients for the indicated                     
!     end point condition                                               
!                                                                       
      IF (IB.EQ.1) THEN 
         CALL ISPL1D (N, XN, FN, ALPHA, BETA, 1, B, C, D, DUMMY (1),    &
         DUMMY (N + 1), DUMMY (2 * N), DUMMY (3 * N - 1), IERR)         
      ELSEIF (IB.EQ.2) THEN 
         CALL ISPL2D (N, XN, FN, ALPHA, BETA, 1, B, C, D, DUMMY (1),    &
         DUMMY (N + 1), DUMMY (2 * N), DUMMY (3 * N - 1), IERR)         
      ELSEIF (IB.EQ.3) THEN 
         CALL ISPL3D (N, XN, FN, ALPHA, BETA, B, C, D, DUMMY (1),       &
         DUMMY (N + 1), IERR)                                           
      ELSEIF (IB.EQ.4) THEN 
         CALL ISPLPE (N, XN, FN, 1, B, C, D, DUMMY (1), DUMMY (N + 2),  &
         DUMMY (2 * N + 2), DUMMY (3 * N + 2), DUMMY (4 * N + 2),       &
         IERR)                                                          
      ELSEIF (IB.EQ.5) THEN 
         CALL ISPLNK (N, XN, FN, B, C, D, DUMMY (1), DUMMY (N + 1),     &
         DUMMY (2 * N), IERR)                                           
      ELSE 
         IERR = - 2 
      ENDIF 
      RETURN 
      END SUBROUTINE ISPLNP                         
