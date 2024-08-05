![KA{P 13}{Akima and Renner Subsplines}                                 
![        {Akima and Renner Subsplines}*)                               
      SUBROUTINE AKIMA1 (N, XN, FN, NK, BETA, B, C, D, ISWTCH, HELP,    &
      IERR)                                                             
!                                                                       
!*****************************************************************      
!                                                                *      
!  The program AKIMA1 computes the coefficients  B(I), C(I) and  *      
!  D(I) for I=0, ... , N-1 of an interpolating cubic AKIMA       *      
!  spline, which can be either periodic or nonperiodic.          *      
!  The subspline has the representation                          *      
!                                                                *      
!  S(X)=FN(I)+B(I)(X-XN(I))+C(I)(X-XN(I))**2+D(I)(X-XN(I))**3    *      
!                                                                *      
!  for any point X in the subinterval [XN(I), XN(I+1)] for       *      
!  I=0,..., N-1.                                                 *      
!                                                                *      
!                                                                *      
!  ASSUMPTIONS:                                                  *      
!  ============                                                  *      
!                1. N >= 4 or NK >= 4                            *      
!                2. If 0.0 < BETA < 1.0, we must have            *      
!                   NK >= N + INT((N+1)/2), otherwise NK = N     *      
!                3. The nodes XN(I), I=0, ..., N, must be        *      
!                   ordered monotonically, i.e., XN(I) < XN(I+1) *      
!                   for I=0, ... , N-1                           *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N       : Index of the last node                              *      
!  XN      : DOUBLE PRECISION (NK+1)-vector XN(0:NK), containing *      
!            the nodes XN(I) for I=0, ... , N                    *      
!  FN      : DOUBLE PRECISION (NK+1)-vector FN(0:NK), containing *      
!            the functional values FN(I) at XN(I) for I=0,...,N  *      
!  NK      : NK = N+INT((N+1)/2) maximal number of nodes allowed *      
!            when using rounded corners. If corners are not      *      
!            rounded, we use: NK = N                             *      
!  BETA    : If 0.0 < BETA < 1.0 the corners are rounded, other- *      
!            wise corners are kept                               *      
!            In the periodic case, we do not round a corner that *      
!            may exist at XN(0) = XN(N), even if 0.0 < BETA < 1.0*      
!  ISWTCH  : = 0, nonperiodic spline                             *      
!            = 1, periodic spline                                *      
!            In the periodic case, the interval [XN(0), XN(N)]   *      
!            must be an interval of periodicity with FN(0) =     *      
!            FN(N).                                              *      
!                                                                *      
!                                                                *      
!  AUXILIARY VARIABLES:                                          *      
!  ====================                                          *      
!  HELP    : DOUBLE PRECISION array HELP(-2:NK+1,1:4)            *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  =================                                             *      
!  N       : Number of the last node. If corners are rounded,    *      
!            i.e., if 0.0 < BETA < 1.0 , then the output value   *      
!            of N can differ from its input. When corners are    *      
!            rounded, the set of nodes can maximally be enlarged *      
!            by INT((N+1)/2) new nodes.                          *      
!  XN      : DOUBLE PRECISION vector XN(0:NK) containing the     *      
!            nodes XN(I), I=0,...,N. If 0.0 < BETA < 1.0, then   *      
!            the output nodes can differ from the input nodes.   *      
!  FN      : DOUBLE PRECISION vector FN(0:NK), containing the    *      
!            functional values FN(I) at XN(I) for I = 0, ..., N. *      
!            If 0.0 < BETA < 1.0, the node and functional values *      
!            XN(I) and FN(I) can differ from their input values. *      
!  B       : DOUBLE PRECISION vector B(0:NK-1) ] B, C and D con- *      
!  C       : DOUBLE PRECISION vector C(0:NK-1) ] tain the coeffi-*      
!  D       : DOUBLE PRECISION vector D(0:NK-1) ] cients of the   *      
!                                     subspline for I=0 to NK-1. *      
!                                                                *      
!  IERR    : error parameter                                     *      
!            = 0,  all is ok                                     *      
!            =-1, N < 4 or NK < 4                                *      
!            =-2, the XN(I) are not monotonically ordered        *      
!            =-3, NK < N+INT((N+1)/2) while 0.0 < BETA < 1.0     *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutine: MACHPD                                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : Gisela Engeln-Mllges                              *      
!  Date     : 04.09.1993                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
!     Declarations                                                      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION XN (0:NK), FN (0:NK), B (0:NK - 1), C (0:NK - 1), &
      D (0:NK - 1), HELP ( - 2:NK + 1, 1:4)                             
!                                                                       
!     Check input parameters                                            
!                                                                       
      IERR = 0 
      IF (N.LT.4.OR.NK.LT.4) THEN 
         IERR = - 1 
         RETURN 
      ENDIF 
      IF (BETA.GT.0.0D0.AND.BETA.LT.1.0D0.AND.NK.LT.N + INT (0.5D0 *    &
      (N + 1) ) ) THEN                                                  
         IERR = - 3 
         RETURN 
      ENDIF 
      IF (ISWTCH.NE.0) THEN 
         ISWTCH = 1 
      ENDIF 
!                                                                       
!     Compute the machine constant                                      
!                                                                       
      FMACHP = 1.0D0 
    5 FMACHP = 0.5D0 * FMACHP 
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 5 
      FMACHP = 2.0D0 * FMACHP 
      EPS = 4.0D0 * FMACHP 
!                                                                       
!     Compute the lengths of the subintervals and store in the first    
!     column of the auxiliary array HELP                                
!                                                                       
      DO 10, I = 0, N - 1 
         HELP (I, 1) = XN (I + 1) - XN (I) 
         IF (HELP (I, 1) .LE.EPS) THEN 
            IERR = - 2 
            RETURN 
         ENDIF 
   10 END DO 
!                                                                       
!     Compute the secant slopes and store in the second column          
!     of the auxiliary array HELP                                       
!                                                                       
      DO 20, I = 0, N - 1 
         HELP (I, 2) = (FN (I + 1) - FN (I) ) / HELP (I, 1) 
   20 END DO 
!                                                                       
!     Compute the magnitude of the slope differences and store          
!     in the third column of HELP                                       
!                                                                       
      DO 30, I = 0, N - 2 
         HELP (I, 3) = DABS (HELP (I + 1, 2) - HELP (I, 2) ) 
   30 END DO 
                                                                        
      IF (DABS (FN (0) - FN (N) ) .GT.EPS) THEN 
         ISWTCH = 0 
      ENDIF 
                                                                        
      IF (BETA.GT.0.0D0.AND.BETA.LT.1.0D0) THEN 
         IF (ISWTCH.EQ.1) THEN 
!                                                                       
!     Prepare loop                                                      
!                                                                       
            I = 1 
            IMAX = N - 1 
            HELP ( - 1, 3) = DABS (HELP (0, 2) - HELP (N - 1, 2) ) 
            HELP (N - 1, 3) = HELP ( - 1, 3) 
         ELSE 
            I = 2 
            IMAX = N - 2 
         ENDIF 
                                                                        
   35    XL = HELP (I - 2, 3) + HELP (I, 3) 
         XR = HELP (I - 1, 3) 
!                                                                       
!     Eliminate existing corners                                        
!                                                                       
         IF (XL.LE.EPS.AND.XR.GT.EPS) THEN 
!                                                                       
!     Relabel points I to N                                             
!                                                                       
            DO 40, J = N, I, - 1 
               XN (J + 1) = XN (J) 
               FN (J + 1) = FN (J) 
   40       END DO 
!                                                                       
!     Reassign interval lengths and slopes from I to N-1                
!                                                                       
            DO 50, J = N - 1, I, - 1 
               HELP (J + 1, 1) = HELP (J, 1) 
               HELP (J + 1, 2) = HELP (J, 2) 
   50       END DO 
!                                                                       
!     Reassign slope differences from I to IMAX                         
!                                                                       
            DO 60, J = IMAX, I, - 1 
               HELP (J + 1, 3) = HELP (J, 3) 
   60       END DO 
!                                                                       
!     Generate new points labelled I and I+1                            
!                                                                       
            XL = HELP (I - 1, 1) 
            XR = HELP (I + 1, 1) 
            XB = BETA * DMIN1 (XL, XR) 
            XLAMDA = XB / XL 
            XMUE = XB / XR 
            XN (I) = XN (I) - XLAMDA * HELP (I - 1, 1) 
            FN (I) = FN (I) - XLAMDA * (FN (I) - FN (I - 1) ) 
            XN (I + 1) = XN (I + 1) + XMUE * HELP (I + 1, 1) 
            FN (I + 1) = FN (I + 1) + XMUE * (FN (I + 2) - FN (I + 1) ) 
!                                                                       
!     Compute new interval lengths                                      
!                                                                       
            DO 70, J = I - 1, I + 1 
               HELP (J, 1) = XN (J + 1) - XN (J) 
   70       END DO 
!                                                                       
!     Compute new slopes and slope differences                          
!                                                                       
            HELP (I, 2) = (FN (I + 1) - FN (I) ) / HELP (I, 1) 
            DO 80, J = I - 1, I 
               HELP (J, 3) = DABS (HELP (J + 1, 2) - HELP (J, 2) ) 
   80       END DO 
!                                                                       
!     Increase number of nodes                                          
!                                                                       
            N = N + 1 
            IMAX = IMAX + 1 
         ENDIF 
!                                                                       
!     Set index for next point                                          
!                                                                       
         I = I + 1 
         IF (I.LE.IMAX) GOTO 35 
      ENDIF 
                                                                        
      IF (ISWTCH.EQ.1) THEN 
!                                                                       
!     For the periodic case, form additional slope data                 
!                                                                       
         HELP ( - 2, 2) = HELP (N - 2, 2) 
         HELP ( - 1, 2) = HELP (N - 1, 2) 
         HELP (N, 2) = HELP (0, 2) 
         HELP (N + 1, 2) = HELP (1, 2) 
      ELSE 
!                                                                       
!     For the nonperiodic case, provide additional slopes as well       
!                                                                       
         HELP ( - 2, 2) = 3.0D0 * HELP (0, 2) - 2.0D0 * HELP (1, 2) 
         HELP ( - 1, 2) = 2.0D0 * HELP (0, 2) - HELP (1, 2) 
         HELP (N, 2) = 2.0D0 * HELP (N - 1, 2) - HELP (N - 2, 2) 
         HELP (N + 1, 2) = 3.0D0 * HELP (N - 1, 2) - 2.0D0 * HELP (N -  &
         2, 2)                                                          
      ENDIF 
!                                                                       
!     Compute additional slope differences                              
!                                                                       
      HELP ( - 2, 3) = DABS (HELP ( - 1, 2) - HELP ( - 2, 2) ) 
      HELP ( - 1, 3) = DABS (HELP (0, 2) - HELP ( - 1, 2) ) 
      HELP (N - 1, 3) = DABS (HELP (N, 2) - HELP (N - 1, 2) ) 
      HELP (N, 3) = DABS (HELP (N + 1, 2) - HELP (N, 2) ) 
!                                                                       
!     Compute the left and right handed slopes at the points 0 and N-1, 
!     and store in column 4 of HELP and in B                            
!                                                                       
      DO 90, I = 0, N - 1 
         XL = HELP (I - 2, 3) 
         XR = HELP (I, 3) 
         IF (XL + XR.GT.EPS) THEN 
            ALPHA = XL / (XL + XR) 
            HELP (I, 4) = HELP (I - 1, 2) + ALPHA * (HELP (I, 2)        &
            - HELP (I - 1, 2) )                                         
            B (I) = HELP (I, 4) 
         ELSE 
            HELP (I, 4) = HELP (I - 1, 2) 
            B (I) = HELP (I, 2) 
         ENDIF 
   90 END DO 
!                                                                       
!     Compute the left handed slope at the point N                      
!                                                                       
      XL = HELP (N - 2, 3) 
      XR = HELP (N, 3) 
      IF (XL + XR.GT.EPS) THEN 
         ALPHA = XL / (XL + XR) 
         HELP (N, 4) = HELP (N - 1, 2) + ALPHA * (HELP (N, 2) - HELP (N &
         - 1, 2) )                                                      
      ELSE 
         HELP (N, 4) = HELP (N - 1, 2) 
      ENDIF 
!                                                                       
!     Compute the coefficients C(I) and D(I)                            
!                                                                       
      DO 100, I = 0, N - 1 
         H = 1.0D0 / HELP (I, 1) 
         C (I) = (3.0D0 * HELP (I, 2) - 2.0D0 * B (I) - HELP (I + 1, 4) &
         ) * H                                                          
         D (I) = (B (I) + HELP (I + 1, 4) - 2.0D0 * HELP (I, 2) )       &
         * H * H                                                        
  100 END DO 
                                                                        
      RETURN 
      END SUBROUTINE AKIMA1                         
