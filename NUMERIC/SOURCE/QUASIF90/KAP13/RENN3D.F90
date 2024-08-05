      SUBROUTINE RENN3D (N, XN, YN, ZN, NK, BETA, B, C, D, T, HELP,     &
      IMARK, IERR)                                                      
!                                                                       
!*****************************************************************      
!                                                                *      
!  The programm RENN3D computes the coefficient vectors B, C and *      
!  D, as well as the lengths T(I) of the parameter intervals of  *      
!  a closed or open parametric 3-dimensional Renner spline.      *      
!  The subsplines are represented in vector form:                *      
!                                                                *      
!  [X(T)]=[PIX(T)]=[XN(I)+B(I,1)*T+C(I,1)*T**2+D(I,1)*T**3]      *      
!  [Y(T)]=[PIY(T)]=[YN(I)+B(I,2)*T+C(I,2)*T**2+D(I,2)*T**3]      *      
!  [Z(T)]=[PIZ(T)]=[ZN(I)+B(I,3)*T+C(I,3)*T**2+D(I,3)*T**3]      *      
!                                                                *      
!  for I=0, ..., N-1 and T a point in the interval [0, T(I)].    *      
!                                                                *      
!                                                                *      
!  ASSUMPTIONS:                                                  *      
!  ============    1. N >= 4 or NK >= 4                          *      
!                  2. If 0.0 < BETA < 1.0: we must have          *      
!                     NK >= N + INT((N+1)/2); otherwise NK = N   *      
!                  2. The node (XN(I), YN(I), ZN(I)) cannot be   *      
!                     equal to  (XN(I+1),YN(I+1), ZN(I+1))       *      
!                     for any I = 0, ..., N-1                    *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N       : Index of the last node                              *      
!  XN      : DOUBLE PRECISION (NK+1)-vector XN(0:NK), containing *      
!            the x-coordinates XN(I) of the nodes for I=0,..., N *      
!  YN      : DOUBLE PRECISION (NK+1)-vector YN(0:NK), containing *      
!            the y-coordinates YN(I) of the nodes for I=0,..., N *      
!  ZN      : DOUBLE PRECISION (NK+1)-vector ZN(0:NK), containing *      
!            the y-coordinates ZN(I) of the nodes for I=0,..., N *      
!  NK      : NK = N + INT((N+1)/2), the maximal number of nodes  *      
!            when using rounded corners,                         *      
!            i.e. for 0.0 < BETA < 1.0.                          *      
!            Without rounding of corners we must have: NK = N.   *      
!  BETA    : for 0.0 < BETA < 1.0 corners are rounded; otherwise *      
!            corners are kept as corners                         *      
!                                                                *      
!                                                                *      
!  AUXILIARY PARAMETERS:                                         *      
!  =====================                                         *      
!  HELP    : DOUBLE PRECISION array HELP(-2:NK+1,1:14)           *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  N       : Index of the last node. If 0.0 < BETA < 1.0 , then  *      
!            the output value for N can differ from its input.   *      
!            When corners are rounded, the node list can be en-  *      
!            larged by at most INT((N+1)/2) points.              *      
!  T       : DOUBLE PRECISION NK-vector T(0:NK-1),] nodes        *      
!            containing the parameter values      ] and          *      
!  XN      : DOUBLE PRECISION vector XN(0:NK)     ] coefficients *      
!  YN      : DOUBLE PRECISION vector YN(0:NK)     ] of           *      
!  ZN      : DOUBLE PRECISION vector ZN(0:NK)     ] the          *      
!  B       : DOUBLE PRECISION array B(0:NK-1,1:3) ] subsplines   *      
!  C       : DOUBLE PRECISION array C(0:NK-1,1:3) ] for          *      
!  D       : DOUBLE PRECISION array D(0:NK-1,1:3) ] I=0, ...,N-1.*      
!            If 0.0 < BETA < 1.0, then the nodes XN(I) and FN(I) *      
!            can be different from their input values for I=0,...*      
!            ..., N                                              *      
!                                                                *      
!  IMARK   : Pointer                                             *      
!            = 0, Subspline is not closed                        *      
!            = 1, Subspline is closed                            *      
!                                                                *      
!  IERR    : Error parameter                                     *      
!            =  0, all is ok                                     *      
!            = -1, N < 4  or  NK < 4                             *      
!            = -2, the assumption                                *      
!                  (XN(I), YN(I), ZN(I)) not equal to            *      
!                  (XN(I+1), YN(I+1), ZN(I+1)) is violated for   *      
!                  some I=0, ..., N-1                            *      
!            = -3, NK < N+INT((N+1)/2) while 0.0 < BETA < 1.0    *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: MACHPD                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author     : Gisela Engeln-Muellges                           *      
!  Date       : 4.7.1993                                         *      
!  Source code: FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
!     Declarations                                                      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION XN (0:NK), YN (0:NK), ZN (0:NK), B (0:NK - 1, 1:3)&
      , C (0:NK - 1, 1:3), D (0:NK - 1, 1:3), T (0:NK - 1), HELP ( - 2: &
      NK + 1, 1:14)                                                     
      DOUBLEPRECISION sklprd 
!                                                                       
!     Check input                                                       
!                                                                       
      IERR = 0 
      IF (N.LT.4.OR.NK.LT.4) THEN 
         IERR = - 1 
         RETURN 
      ENDIF 
!                                                                       
!     If  0.0 < BETA < 1.0  corners are rounded                         
!                                                                       
      IF (BETA.GT.0.0D0.AND.BETA.LT.1.0D0.AND.NK.LT.N + INT (0.5D0 *    &
      (N + 1) ) ) THEN                                                  
         IERR = - 3 
         RETURN 
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
!     Calculate the chordal vectors; store their three coordinates      
!     in the first three columns of the auxiliary array HELP,           
!     the fourth column of HELP contains the length of these chordal vec
!                                                                       
      DO 10, I = 0, N - 1 
         HELP (I, 1) = XN (I + 1) - XN (I) 
         HELP (I, 2) = YN (I + 1) - YN (I) 
         HELP (I, 3) = ZN (I + 1) - ZN (I) 
         HELP (I, 4) = DSQRT (HELP (I, 1) * HELP (I, 1) + HELP (I, 2)   &
         * HELP (I, 2) + HELP (I, 3) * HELP (I, 3) )                    
!                                                                       
!     we check that all chordal vectors are nonzero, thus               
!     verifying that consecutive nodes are distinct.                    
!                                                                       
         IF (HELP (I, 4) .LE.EPS) THEN 
            IERR = - 2 
            RETURN 
         ELSE 
!                                                                       
!    Columns 5 to 7 of HELP contain the three components                
!    of the the chordal unit vectors                                    
!                                                                       
            HELP (I, 5) = HELP (I, 1) / HELP (I, 4) 
            HELP (I, 6) = HELP (I, 2) / HELP (I, 4) 
            HELP (I, 7) = HELP (I, 3) / HELP (I, 4) 
         ENDIF 
   10 END DO 
!                                                                       
!     We compute the area of the parallelograms generated by            
!     two consecutive chordal vectors and store the results in          
!     column 8 of the auxiliary array HELP                              
!                                                                       
      DO 20, I = 0, N - 2 
         sklprd = HELP (I, 5) * HELP (I + 1, 5) + HELP (I, 6) * HELP (I &
         + 1, 6) + HELP (I, 7) * HELP (I + 1, 7)                        
         HELP (I, 8) = DSQRT (1.0D0 - sklprd * sklprd) 
   20 END DO 
!                                                                       
!     Determine whether the curve is closed or not                      
!                                                                       
      IMARK = 0 
      IF (DABS (XN (0) - XN (N) ) .LE.EPS.AND.DABS (YN (0) - YN (N) )   &
      .LE.EPS.AND.DABS (ZN (0) - ZN (N) ) .LE.EPS) THEN                 
         IMARK = 1 
      ENDIF 
!                                                                       
      IF (BETA.GT.0.0D0.AND.BETA.LT.1.0D0) THEN 
!                                                                       
!     If 0.0 < BETA < 1.0, round existing corners                       
!                                                                       
         IF (IMARK.EQ.1) THEN 
!                                                                       
!     Curve is closed                                                   
!                                                                       
            XL = HELP (N - 2, 8) + HELP (0, 8) 
            sklprd = HELP (N - 1, 5) * HELP (0, 5) + HELP (N - 1, 6)    &
            * HELP (0, 6) + HELP (N - 1, 7) * HELP (0, 7)               
            XR = DSQRT (1.0D0 - sklprd * sklprd) 
!                                                                       
!     Condition for a corner                                            
!                                                                       
            IF (XL.LE.EPS.AND.XR.GT.EPS) THEN 
!                                                                       
!     Relabel points                                                    
!                                                                       
               DO 30, I = 0, N - 1 
                  XN (I) = XN (I + 1) 
                  YN (I) = YN (I + 1) 
                  ZN (I) = ZN (I + 1) 
   30          END DO 
               XN (N) = XN (0) 
               YN (N) = YN (0) 
               ZN (N) = ZN (0) 
!                                                                       
!     Relabel chordal data                                              
!                                                                       
               DO 40 K = 1, 7 
                  HELP (N, K) = HELP (0, K) 
   40          END DO 
               DO 50, I = 0, N - 1 
                  DO 60, K = 1, 7 
                     HELP (I, K) = HELP (I + 1, K) 
   60             END DO 
   50          END DO 
!                                                                       
!     Relabel areas                                                     
!                                                                       
               DO 70, I = 0, N - 3 
                  HELP (I, 8) = HELP (I + 1, 8) 
   70          END DO 
               sklprd = HELP (N - 2, 5) * HELP (N - 1, 5) + HELP (N - 2,&
               6) * HELP (N - 1, 6) + HELP (N - 2, 7) * HELP (N - 1, 7) 
               HELP (N - 2, 8) = DSQRT (1.0D0 - sklprd * sklprd) 
            ENDIF 
!                                                                       
!     Prepare loop                                                      
!                                                                       
            I = 1 
            IMAX = N - 1 
            sklprd = HELP (N - 1, 5) * HELP (0, 5) + HELP (N - 1, 6)    &
            * HELP (0, 6) + HELP (N - 1, 7) * HELP (0, 7)               
            HELP ( - 1, 8) = DSQRT (1.0D0 - sklprd * sklprd) 
            HELP (N - 1, 8) = HELP ( - 1, 8) 
         ELSE 
!                                                                       
!     Curve not closed                                                  
!                                                                       
!     Prepare loop                                                      
!                                                                       
            I = 2 
            IMAX = N - 2 
         ENDIF 
                                                                        
   75    XL = HELP (I - 2, 8) + HELP (I, 8) 
         XR = HELP (I - 1, 8) 
!                                                                       
!     Round existing corners                                            
!                                                                       
         IF (XL.LE.EPS.AND.XR.GT.EPS) THEN 
!                                                                       
!     Relabel points  I to N                                            
!                                                                       
            DO 80, J = N, I, - 1 
               XN (J + 1) = XN (J) 
               YN (J + 1) = YN (J) 
               ZN (J + 1) = ZN (J) 
   80       END DO 
!                                                                       
!    Relabel chordal data for I to N-1                                  
!                                                                       
            DO 90, J = N - 1, I, - 1 
               DO 100, K = 1, 7 
                  HELP (J + 1, K) = HELP (J, K) 
  100          END DO 
   90       END DO 
!                                                                       
!     Relabel areas I to IMAX                                           
!                                                                       
            DO 110, J = IMAX, I, - 1 
               HELP (J + 1, 8) = HELP (J, 8) 
  110       END DO 
!                                                                       
!     Create two new nodes  I and I+1                                   
!                                                                       
            XL = HELP (I - 1, 4) 
            XR = HELP (I + 1, 4) 
            XB = BETA * DMIN1 (XL, XR) 
            XLAMDA = XB / XL 
            XMUE = XB / XR 
                                                                        
            XN (I) = XN (I) - XLAMDA * HELP (I - 1, 1) 
            YN (I) = YN (I) - XLAMDA * HELP (I - 1, 2) 
            ZN (I) = ZN (I) - XLAMDA * HELP (I - 1, 3) 
            XN (I + 1) = XN (I + 1) + XMUE * HELP (I + 1, 1) 
            YN (I + 1) = YN (I + 1) + XMUE * HELP (I + 1, 2) 
            ZN (I + 1) = ZN (I + 1) + XMUE * HELP (I + 1, 3) 
!                                                                       
!     Compute new chordal data at I-1, I und I+1 and store              
!                                                                       
            DO 120, J = I - 1, I + 1 
               HELP (J, 1) = XN (J + 1) - XN (J) 
               HELP (J, 2) = YN (J + 1) - YN (J) 
               HELP (J, 3) = ZN (J + 1) - ZN (J) 
               HELP (J, 4) = DSQRT (HELP (J, 1) * HELP (J, 1) + HELP (J,&
               2) * HELP (J, 2) + HELP (J, 3) * HELP (J, 3) )           
               HELP (J, 5) = HELP (J, 1) / HELP (J, 4) 
               HELP (J, 6) = HELP (J, 2) / HELP (J, 4) 
               HELP (J, 7) = HELP (J, 3) / HELP (J, 4) 
  120       END DO 
!                                                                       
!     Compute new areas at  I-2, I-1, I and I+1, and store              
!                                                                       
            DO 130, J = I - 2, I + 1 
               sklprd = HELP (J, 5) * HELP (J + 1, 5) + HELP (J, 6)     &
               * HELP (J + 1, 6) + HELP (J, 7) * HELP (J + 1, 7)        
               HELP (J, 8) = DSQRT (1.0D0 - sklprd * sklprd) 
  130       END DO 
!                                                                       
!     Increase point count                                              
!                                                                       
            N = N + 1 
            IMAX = IMAX + 1 
         ENDIF 
!                                                                       
!     Set index for next point                                          
!                                                                       
         I = I + 1 
         IF (I.LE.IMAX) GOTO 75 
!                                                                       
      ENDIF 
!                                                                       
      IF (IMARK.EQ.1) THEN 
!                                                                       
!     Curve is closed                                                   
!                                                                       
!     Compute additional chordal data and store                         
!                                                                       
         DO 140, I = N - 2, N - 1 
            DO 150, K = 1, 7 
               HELP (I - N, K) = HELP (I, K) 
  150       END DO 
  140    END DO 
         DO 160, I = 0, 1 
            DO 170, K = 1, 7 
               HELP (I + N, K) = HELP (I, K) 
  170       END DO 
  160    END DO 
      ELSE 
!                                                                       
!     Curve not closed                                                  
!                                                                       
         DO 180, K = 1, 3 
            HELP ( - 2, K) = 3.D0 * HELP (0, K) - 2.D0 * HELP (1, K) 
            HELP ( - 1, K) = 2.D0 * HELP (0, K) - HELP (1, K) 
  180    END DO 
!                                                                       
!     Compute additional chordal data at -2, -1, N and N+1              
!                                                                       
         DO 190, I = - 2, - 1 
            HELP (I, 4) = DSQRT (HELP (I, 1) * HELP (I, 1) + HELP (I, 2)&
            * HELP (I, 2) + HELP (I, 3) * HELP (I, 3) )                 
            IF (HELP (I, 4) .GT.EPS) THEN 
               HELP (I, 5) = HELP (I, 1) / HELP (I, 4) 
               HELP (I, 6) = HELP (I, 2) / HELP (I, 4) 
               HELP (I, 7) = HELP (I, 3) / HELP (I, 4) 
            ELSE 
               HELP (I, 5) = 0.0D0 
               HELP (I, 6) = 0.0D0 
               HELP (I, 7) = 0.0D0 
            ENDIF 
  190    END DO 
         DO 200, K = 1, 3 
            HELP (N, K) = 2.D0 * HELP (N - 1, K) - HELP (N - 2, K) 
            HELP (N + 1, K) = 3.D0 * HELP (N - 1, K) - 2.D0 * HELP (N - &
            2, K)                                                       
  200    END DO 
         DO 210, I = N, N + 1 
            HELP (I, 4) = DSQRT (HELP (I, 1) * HELP (I, 1) + HELP (I, 2)&
            * HELP (I, 2) + HELP (I, 3) * HELP (I, 3) )                 
            IF (HELP (I, 4) .GT.EPS) THEN 
               HELP (I, 5) = HELP (I, 1) / HELP (I, 4) 
               HELP (I, 6) = HELP (I, 2) / HELP (I, 4) 
               HELP (I, 7) = HELP (I, 3) / HELP (I, 4) 
            ELSE 
               HELP (I, 5) = 0.0D0 
               HELP (I, 6) = 0.0D0 
               HELP (I, 7) = 0.0D0 
            ENDIF 
  210    END DO 
      ENDIF 
!                                                                       
!     Compute new areas for  -2, -1, N-1 and N                          
!                                                                       
      DO 220, I = - 2, - 1 
         sklprd = HELP (I, 5) * HELP (I + 1, 5) + HELP (I, 6) * HELP (I &
         + 1, 6) + HELP (I, 7) * HELP (I + 1, 7)                        
         HELP (I, 8) = DSQRT (1.0D0 - sklprd * sklprd) 
  220 END DO 
      DO 230, I = N - 1, N 
         sklprd = HELP (I, 5) * HELP (I + 1, 5) + HELP (I, 6) * HELP (I &
         + 1, 6) + HELP (I, 7) * HELP (I + 1, 7)                        
         HELP (I, 8) = DSQRT (1.0D0 - sklprd * sklprd) 
  230 END DO 
!                                                                       
!     Compute the left and right hand unit tangent vectors and store in 
!     columns 9 to 11 and 12 to 14 of HELP, respectively                
!                                                                       
      DO 240, I = 0, N 
         XL = HELP (I - 2, 8) 
         XR = HELP (I, 8) 
         IF (XL + XR.GT.EPS) THEN 
            ALPHA = XL / (XL + XR) 
            DO 250, K = 1, 3 
               HELP (I, K + 8) = HELP (I - 1, K) + ALPHA * (HELP (I, K) &
               - HELP (I - 1, K) )                                      
  250       END DO 
            HELPT = DSQRT (HELP (I, 9) * HELP (I, 9) + HELP (I, 10)     &
            * HELP (I, 10) + HELP (I, 11) * HELP (I, 11) )              
            HELP (I, 9) = HELP (I, 9) / HELPT 
            HELP (I, 10) = HELP (I, 10) / HELPT 
            HELP (I, 11) = HELP (I, 11) / HELPT 
            HELP (I, 12) = HELP (I, 9) 
            HELP (I, 13) = HELP (I, 10) 
            HELP (I, 14) = HELP (I, 11) 
         ELSE 
            DO 260, K = 5, 7 
               HELP (I, K + 4) = HELP (I - 1, K) 
               HELP (I, K + 7) = HELP (I, K) 
  260       END DO 
         ENDIF 
  240 END DO 
!                                                                       
!     Compute the lengths of the parameter intervals                    
!                                                                       
      DO 270, I = 0, N - 1 
         TS1 = HELP (I, 12) + HELP (I + 1, 9) 
         TS2 = HELP (I, 13) + HELP (I + 1, 10) 
         TS3 = HELP (I, 14) + HELP (I + 1, 11) 
         A1 = 16.0D0 - (TS1 * TS1 + TS2 * TS2 + TS3 * TS3) 
         B1 = 6.0D0 * (HELP (I, 1) * TS1 + HELP (I, 2) * TS2 + HELP (I, &
         3) * TS3)                                                      
         C1 = 36.0D0 * HELP (I, 4) * HELP (I, 4) 
         T (I) = ( - B1 + DSQRT (B1 * B1 + A1 * C1) ) / A1 
         HELPT = 1.0D0 / T (I) 
!                                                                       
!     Compute the spline coefficient vectors                            
!                                                                       
         DO 280, K = 1, 3 
            B (I, K) = HELP (I, K + 11) 
            C (I, K) = (3.0D0 * HELPT * HELP (I, K) - 2.0D0 * B (I, K)  &
            - HELP (I + 1, K + 8) ) * HELPT                             
            D (I, K) = (B (I, K) + HELP (I + 1, K + 8) - 2.0D0 * HELPT *&
            HELP (I, K) ) * HELPT * HELPT                               
  280    END DO 
  270 END DO 
!                                                                       
      RETURN 
      END SUBROUTINE RENN3D                         
