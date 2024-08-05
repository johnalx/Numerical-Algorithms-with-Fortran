      SUBROUTINE ISPLTR (N, XN, FN, MS, PX, PY, R, B, C, D, PHIN, PHIR, &
      DUMMY, IERR)                                                      
!                                                                       
!*****************************************************************      
!                                                                *      
!  ISPLTR computes the coefficients R(I), B(I), C(I), D(I) for   *      
!  I=0,1,...,N-1 of a transformed parametric cubic interpolating *      
!  spline for a closed everywhere smooth curve K.                *      
!  This program determines the transformed coordinates PHIN(I)   *      
!  and R(I) from the given nodes XN(I), FN(I) for I=0,1,...,N.   *      
!  The new points PHIN(I) and R(I) define a nonparametric        *      
!  periodic interpolating spline S(PHI) of the form:             *      
!                                                                *      
!  S(PHI) = R(I) + B(I)(PHI-PHIN(I)) + C(I)(PHI-PHIN(I))**2 +    *      
!                                    + D(I)(PHI-PHIN(I))**3      *      
!                                                                *      
!  for PHI in the interval [PHIN(I),PHIN(I+1)] for I=0,1,...,N-1.*      
!                                                                *      
!  Since the nodes PHIN(I) must be ordered monotonically, we must*      
!  translate and rotate the origin to P=(PX,PY) by an angle PHIR.*      
!  This is clearly dependent on the original nodes XN(I), FN(I). *      
!  In order to be able to achieve monotonicity for the PHIN(I),  *      
!  the following must be true for the original nodes:            *      
!   - the new origin P must lie in the region F defines by XN(I),*      
!     FN(I) so that every ray from P intersects the boundary     *      
!     curve of F precisely once.                                 *      
!     (The coordinates PX and PY of P should be specified by the *      
!     user. See input parameter MS).                             *      
!   - the points XN(I), FN(I) must be ordered so that the        *      
!     boundary curve of F is traversed counterclockwise from     *      
!     XN(0), FN(0) to XN(N), FN(N). For periodicity of S(PHI),   *      
!     the endpoint  XN(N), FN(N) must be equal to XN(0), FN(0).  *      
!                                                                *      
!  The coordinates of PHIN(I), R(I) are computed as:             *      
!         PHIN(0) = 0,                                           *      
!         PHIN(I) = ATAN(Y'(I)/X'(I)) - PHIR , I=1,...,N-1,      *      
!         PHIN(N) = 2*PI,                                        *      
!         R(I) = SQRT(X'(I)**2 + Y'(I)**2), I=0,1,...,N,         *      
!         MIT:  PHIR = ATAN(FN(0)/XN(0)),                        *      
!               X'(I) = XN(I) - PX,  Y'(I) = FN(I) - PY.         *      
!                                                                *      
!                                                                *      
!  REMARK:  The curve K can be evaluated using SUBROUTINE        *      
!  =======  TSPANA. A table of values for K can be generated     *      
!           with SUBROUTINE TSPTAB.                              *      
!           Both subroutines need the parameters PX, PY and      *      
!           PHIR from this subroutine in order to transform      *      
!           the coordinates.                                     *      
!                                                                *      
!                                                                *      
!  ASSUMPTIOS:    1.         N > 2                               *      
!  ===========    2.     XN(0) = XN(N)                           *      
!                 3.     FN(0) = FN(N)                           *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N  :  index of the final node                                 *      
!  XN :  vector XN(0:N); the nodes for I = 0,1,...,N             *      
!  FN :  vector FN(0:N); the values at XN(I)                     *      
!                                                                *      
!  MS :  Index for the translation of the origin:                *      
!        MS > 0 : The user specifies the coordinates PX, PY      *      
!        MS = 0 : no shift, i.e., PX=PY=0                        *      
!        MS < 0 : the shift coordinates PX, PY are computed in   *      
!                 SUBROUTINE ISPLTR as:                          *      
!                    PX = (XMAX+XMIN)/2                          *      
!                    PY = (YMAX+YMIN)/2                          *      
!                    where : XMAX = DMAX1(XN(I)) ]               *      
!                            XMIN = DMIN1(XN(I)) ] I=0,1,...,N   *      
!                            YMAX = DMAX1(FN(I)) ]               *      
!                            YMIN = DMIN1(FN(I)) ]               *      
!                 Please note: This does not insure that the     *      
!                 point P will be properly centerd as required.  *      
!                 If this is not the case IERR = -3 will result. *      
!                                                                *      
!  PX : ]  the coordinates of P                                  *      
!  PY : ]  (for MS > 0)                                          *      
!                                                                *      
!                                                                *      
!  AUXILIARY VARIABLES:                                          *      
!  ====================                                          *      
!  DUMMY :  vector DUMMY(1:5*N+1)                                *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  R  :  ] vectors ..(0:N);                                      *      
!  B  :  ] the elements in positions 0 to N-1 are the            *      
!  C  :  ] coefficients of the spline function S(PHI).           *      
!  D  :  ] B(N), C(N), D(N) are used for storage.                *      
!          The elements R(I), I=0,1,...,N contain the length of  *      
!          the vectors with angle PHIN(I).                       *      
!                                                                *      
!  PHIN :  vector PHIN(0:N); PHIN(I) denotes the angle of        *      
!          (XN(I),FN(I)) as seen from P=(PX,PY).                 *      
!                                                                *      
!  PX   : ] coordinates of P                                     *      
!  PY   : ]                                                      *      
!  PHIR :  rotation angle of the coordinate system in radians    *      
!  IERR :  error parameter                                       *      
!          =  0 :  All is ok                                     *      
!          = -1 :  N < 3                                         *      
!          = -3 :  non-monotonous nodes PHIN(I);                 *      
!                  PHIN(I) >= PHIN(I+1) for some I=0,1,...,N-1   *      
!          = -4 :  XN(0) not equal to XN(N)  or                  *      
!                  FN(0) not equal to FN(N)                      *      
!          =  1 :  crash in CYTSY, system matrix numerically     *      
!                  singular.                                     *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Subroutines required: ISPLPE                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : GÅnter Palm                                        *      
!  Date     : 04.15.1988                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION XN (0:N), FN (0:N), R (0:N), B (0:N), C (0:N),    &
      D (0:N), PHIN (0:N), DUMMY (1:5 * N + 1)                          
!                                                                       
!-----Check input data                                                  
!                                                                       
      IERR = - 1 
      IF (N.LT.3) RETURN 
      IERR = - 4 
      IF (XN (0) .NE.XN (N) ) RETURN 
      IF (FN (0) .NE.FN (N) ) RETURN 
!                                                                       
!-----Initialize                                                        
!                                                                       
      PI = 4.0D0 * DATAN (1.0D0) 
      PI2 = 2.0D0 * PI 
!                                                                       
!-----Transform the coordinates                                         
!                                                                       
!     If MS differs from 0, translate the origin to P=(PX,PY),          
!     by changing the coordinates of X(I), Y(I) by PX and PY.           
!                                                                       
      IF (MS.EQ.0) THEN 
!                                                                       
!       no shift                                                        
!                                                                       
         PX = 0.0D0 
         PY = 0.0D0 
         DO 10 I = 0, N, 1 
            B (I) = XN (I) 
            C (I) = FN (I) 
   10    END DO 
      ELSE 
!                                                                       
!       shift by  (PX,PY)                                               
!                                                                       
         IF (MS.LT.0) THEN 
!                                                                       
!         Compute PX and PY                                             
!                                                                       
            XMAX = XN (0) 
            XMIN = XN (0) 
            YMAX = FN (0) 
            YMIN = FN (0) 
            DO 20 I = 1, N, 1 
               XMAX = DMAX1 (XMAX, XN (I) ) 
               XMIN = DMIN1 (XMIN, XN (I) ) 
               YMAX = DMAX1 (YMAX, FN (I) ) 
               YMIN = DMIN1 (YMIN, FN (I) ) 
   20       END DO 
            PX = (XMAX + XMIN) / 2.0D0 
            PY = (YMAX + YMIN) / 2.0D0 
         ENDIF 
         DO 30 I = 0, N, 1 
            B (I) = XN (I) - PX 
            C (I) = FN (I) - PY 
   30    END DO 
      ENDIF 
!                                                                       
!-----Compute the transformed nodes                                     
!                                                                       
!     1. Compute R(I); Return if R(I)=0, i.e., when                     
!        (PX,PY) is one of the nodes                                    
!                                                                       
      IERR = - 3 
      DO 40 I = 0, N, 1 
         R (I) = DSQRT (B (I) * B (I) + C (I) * C (I) ) 
         IF (R (I) .EQ.0.0D0) RETURN 
   40 END DO 
!                                                                       
!     2. Compute the coordinates X',Y' rotated by ALPHA from:           
!                                                                       
!           [ X']   [ COS(ALPHA)   -SIN(ALPHA) ] [ X ]                  
!           [   ] = [                          ] [   ]                  
!           [ Y']   [ SIN(ALPHA)    COS(ALPHA) ] [ Y ]                  
!                                                                       
!        with  ALPHA = -PHIR                                            
!                                                                       
      PHIR = DACOS (B (0) / R (0) ) 
      IF (C (0) .LT.0.0D0) PHIR = PI2 - PHIR 
      CA = B (0) / R (0) 
      SA = - C (0) / R (0) 
      DO 50 I = 0, N 
         D (I) = B (I) * CA - C (I) * SA 
         C (I) = B (I) * SA + C (I) * CA 
   50 END DO 
!                                                                       
!     3. Compute the angular coordinate PHIN(I); Return                 
!        if the angles do not increase strictly monotonically           
!                                                                       
      PHIN (0) = 0.0D0 
      DO 60 I = 1, N - 1 
         PHIN (I) = DACOS (D (I) / R (I) ) 
         IF (C (I) .LT.0.0D0) PHIN (I) = PI2 - PHIN (I) 
         IF (PHIN (I) .LE.PHIN (I - 1) ) RETURN 
   60 END DO 
      PHIN (N) = PI2 
!                                                                       
!-----Compute the spline coefficients                                   
!                                                                       
      CALL ISPLPE (N, PHIN, R, 1, B, C, D, DUMMY (1), DUMMY (N + 2),    &
      DUMMY (2 * N + 2), DUMMY (3 * N + 2), DUMMY (4 * N + 2), IERR)    
      RETURN 
      END SUBROUTINE ISPLTR                         
