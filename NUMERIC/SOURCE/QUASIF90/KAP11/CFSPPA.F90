![  {Parametric Cubic Fitting Splines}                                  
![  {Parametric Cubic Fitting Splines}*)                                
      SUBROUTINE CFSPPA (N, XN, FN, WX, WY, T, MT, IB, ALPHA, BETA, MW, &
      AX, BX, CX, DX, AY, BY, CY, DY, AUXF, IERR)                       
!                                                                       
!*****************************************************************      
!                                                                *      
!  CFSPPA computes the coefficients AX(I), BX(I), CX(I), DX(I),  *      
!  AY(I), BY(I), CY(I), DY(I) for I=0, 1, ..., N-1 of a          *      
!  parametric cubic fitting spline for various end point con-    *      
!  ditions. The end point conditions are prescribed via the      *      
!  parameter IB.                                                 *      
!  The parametric spline S with parameters T(I) for I=0, 1, ..., *      
!  N, has two component functions SX and SY of the following form*      
!                                                                *      
!  SX := SX(T) = AX(I) + BX(I)(T-T(I)) + CX(I)(T-T(I))**2 +      *      
!                                      + DX(I)(T-T(I))**3        *      
!                                                                *      
!  SY := SY(T) = AY(I) + BY(I)(T-T(I)) + CY(I)(T-T(I))**2 +      *      
!                                      + DY(I)(T-T(I))**3        *      
!                                                                *      
!  for T in the interval [T(I),T(I+1)], I=0, 1, ..., N-1.        *      
!                                                                *      
!  SX and SY are nonparametric cubic splines.                    *      
!                                                                *      
!                                                                *      
!  ASSUMPTIONS:    1.         N > 4     , for IB = 1, 2 or 3     *      
!  ============               N > 5     , for IB =  4            *      
!                  2.      T(I) < T(I+1), I=0, 1, ..., N-1       *      
!                  3.     WX(I) > 0.0   ] I=0, 1, ..., N         *      
!                  3.     WY(I) > 0.0   ]                        *      
!                  4.     WX(0) = WX(N) ] for IB = 4             *      
!                         WY(0) = WY(N) ]                        *      
!                  5.     FN(0) = FN(N) , for IB = 4             *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N  :  Index of the last node                                  *      
!  XN :  vector XN(0:N); XN(I) is the Ith node, I = 0, ..., N    *      
!  FN :  vector FN(0:N); FN(I) is the data at the node XN(I)     *      
!  WX :  vector WX(0:N); weight for the node XN(I), I=0, ..., N  *      
!  WY :  vector WXY0:N); weight for the data point FN(I),        *      
!        I = 0, 1, ..., N                                        *      
!  T  :  vector T(0:N); the parameter value corresponding to     *      
!        XN(I),FN(I)                                             *      
!  MT :  label for the assignment of the parameter T(I)          *      
!        MT = 0 : The user will prescribe the parameter values   *      
!                 T(I) for I=0, 1, ..., N.                       *      
!        MT = 1 : The parameter values shall be computed in sub- *      
!                 routine PSPPV from the chordal lengths.        *      
!        MT = 2 : The parameter values shall be computed in sub- *      
!                 routine PSPPV from the arc length.             *      
!                                                                *      
!  IB : determines the end point condition:                      *      
!       IB = 1: first end point derivative with respect to the   *      
!               parameter prescribed                             *      
!       IB = 2: second end point derivative with respect to the  *      
!               parameter prescribed                             *      
!       IB = 3: first end point derivative DY/DX prescribed      *      
!       IB = 4: periodic spline                                  *      
!                                                                *      
!  ALPHA :  vector ALPHA(1:2)                                    *      
!  BETA  :  vector BETA(1:2)                                     *      
!           for IB = 1 or 2 : first or second derivative wrt     *      
!                             parameter:                         *      
!                             ALPHA(1)=SX(IB)(T(0))              *      
!                             ALPHA(2)=SY(IB)(T(0))              *      
!                             BETA(1) =SX(IB)(T(N))              *      
!                             BETA(2) =SY(IB)(T(N))              *      
!           for IB = 3 :  first end point derivative DY/DX       *      
!                           ALPHA(1) = DY/DX (XN(0))             *      
!                           BETA(1)  = DY/DX (XN(N))             *      
!                           ALPHA(2) : not needed                *      
!                           BETA(2)  : not needed                *      
!                If the magnitude of ALPHA(1) or BETA(1) exceeds *      
!                1.E10, then the corresponding tangent vector is *      
!                computed as follows:                            *      
!                 . .                                            *      
!                (X,Y) = (0,SIGN(1,FN(1)-FN(0))  (left end point)*      
!                 . .                                            *      
!                (X,Y) = (0,SIGN(1,FN(N)-FN(N-1)) (right end     *      
!                                                        point)  *      
!                                                                *      
!           for IB = 4 : not needed                              *      
!                                                                *      
!  (A natural parametric fitting spline can be obtained for IB=2 *      
!   and ALPHA(1) = ALPHA(2) = BETA(1) = BETA(2) = 0.0.)          *      
!                                                                *      
!  MW : label for preassignment of weights WX(I),WY(I), I=0,...,N*      
!       MW < 1: weights not preassigned,                         *      
!               CFSPPA uses  WX(I) = WY(I) = 1.0                 *      
!       MW = 1: the user preassigns the weights WY(I),           *      
!               CFSPPA uses WX(I) = WY(I)                        *      
!       MW > 1: the user preassigns both WX(I) and WY(I)         *      
!                                                                *      
!                                                                *      
!  AUXILIARY VARIABLES:                                          *      
!  ====================                                          *      
!  AUXF :  vector AUXF(1:14*N-10)                                *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  AX :  vector  AX(0:N) ]   The elements in positions 0 to N-1  *      
!  BX :  vector  BX(0:N) ]   are the coefficients of the com-    *      
!  CX :  vector  CX(0:N) ]   component spline function SX        *      
!  DX :  vector  DX(0:N) ]                                       *      
!                                                                *      
!  AY :  vector  AY(0:N) ]   The elements in positions 0 to N-1  *      
!  BY :  vector  BY(0:N) ]   are the coefficients of the com-    *      
!  CY :  vector  CY(0:N) ]   component spline function SY        *      
!  DY :  vector  DY(0:N) ]                                       *      
!                            The elements in position N are      *      
!                            auxiliary variables                 *      
!  IERR :  error parameter                                       *      
!          =  0 :  All is o.k                                    *      
!          = -1 :  N < 5  (for IB = 1, 2, 3)                     *      
!                  N < 6  (for IB =  4 )                         *      
!          = -2 :  IB < 1  or  IB > 4                            *      
!          = -3 :  Inadmissable weight WX or WY                  *      
!          = -4 :  parameter values T(I) not ordered monotonical-*      
!                  ly, i.e., T(I) >= T(I+1) for some I=0,..., N-1*      
!          = -5 :  IB = 4 and FN(0) not equal to FN(N) or        *      
!                  WX(0) not equal to WX(N) or                   *      
!                  WY(0) not equal to WY(N)                      *      
!          =  1 :  error in FDISY or NCYFSY (system matrix       *      
!                  numerically singular)                         *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: CFSP1D, CFSP2D, CFSPPE, PSPPV           *      
!                                                                *      
!                                                                *      
!  Reference : Engeln-MÅllges, G.; Reutter, F., [ENGE87].        *      
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
      DOUBLEPRECISION XN (0:N), FN (0:N), WX (0:N), WY (0:N), T (0:N),  &
      AX (0:N), BX (0:N), CX (0:N), DX (0:N), AY (0:N), BY (0:N),       &
      CY (0:N), DY (0:N), ALPHA (2), BETA (2), AUXF (1:14 * N - 10)     
!                                                                       
!-----Check assumptions                                                 
!                                                                       
      IERR = - 1 
      IF (N.LT.5) RETURN 
      IF (IB.LT.1.OR.IB.GT.4) THEN 
         IERR = - 2 
         RETURN 
      ENDIF 
      IERR = - 3 
!                                                                       
!-----Check or compute the weights                                      
!                                                                       
      IF (MW.GT.1) THEN 
!                                                                       
!       Check the prescribed weights                                    
!                                                                       
         DO 10 I = 0, N, 1 
            IF (WX (I) .LE.0.0D0.OR.WY (I) .LE.0.0D0) RETURN 
   10    END DO 
         MREP = 1 
      ELSE 
         IF (MW.EQ.1) THEN 
!                                                                       
!         Check the prescribed weights WY                               
!                                                                       
            DO 20 I = 0, N, 1 
               IF (WY (I) .LE.0.0D0) RETURN 
   20       END DO 
         ELSE 
!                                                                       
!         If not prescribed, assign 1.0 to all weights WY               
!                                                                       
            DO 30 I = 0, N, 1 
               WY (I) = 1.0D0 
   30       END DO 
         ENDIF 
!                                                                       
!       Set the weights  WX equal to the weights WY                     
!                                                                       
         DO 40 I = 0, N, 1 
            WX (I) = WY (I) 
   40    END DO 
         MREP = 2 
      ENDIF 
!                                                                       
!-----Compute and/or check the parameter values                         
!                                                                       
      IF (MT.GT.0) THEN 
!                                                                       
!       Compute the parameter values in subroutine PSPPV                
!                                                                       
         CALL PSPPV (N, XN, FN, T, MT, IERR) 
         IF (IERR.NE.0) THEN 
            IERR = - 4 
            RETURN 
         ENDIF 
      ELSE 
!                                                                       
!       Check the prescribed parameter values                           
!                                                                       
         IERR = - 4 
         DO 50 I = 0, N - 1, 1 
            IF (T (I + 1) .LE.T (I) ) RETURN 
   50    END DO 
      ENDIF 
!                                                                       
!-----Compute the spline coefficients                                   
!                                                                       
      IF (IB.EQ.1) THEN 
         CALL CFSP1D (N, T, XN, WX, ALPHA (1), BETA (1), 1, AX, BX, CX, &
         DX, AUXF (1), AUXF (N + 1), AUXF (2 * N + 1), AUXF (3 * N + 1),&
         AUXF (4 * N), AUXF (5 * N - 1), AUXF (6 * N - 2), IERR)        
         IF (IERR.NE.0) RETURN 
         CALL CFSP1D (N, T, FN, WY, ALPHA (2), BETA (2), MREP, AY, BY,  &
         CY, DY, AUXF (1), AUXF (N + 1), AUXF (2 * N + 1), AUXF (3 * N +&
         1), AUXF (4 * N), AUXF (5 * N - 1), AUXF (6 * N - 2), IERR)    
      ELSEIF (IB.EQ.2) THEN 
         CALL CFSP2D (N, T, XN, WX, ALPHA (1), BETA (1), 1, AX, BX, CX, &
         DX, AUXF (1), AUXF (N + 1), AUXF (2 * N + 1), AUXF (3 * N + 1),&
         AUXF (4 * N), AUXF (5 * N - 1), AUXF (6 * N - 2), IERR)        
         IF (IERR.NE.0) RETURN 
         CALL CFSP2D (N, T, FN, WY, ALPHA (2), BETA (2), MREP, AY, BY,  &
         CY, DY, AUXF (1), AUXF (N + 1), AUXF (2 * N + 1), AUXF (3 * N +&
         1), AUXF (4 * N), AUXF (5 * N - 1), AUXF (6 * N - 2), IERR)    
      ELSEIF (IB.EQ.3) THEN 
!                                                                       
!       Compute the tangent vectors of the first derivatives            
!       for SX and SY                                                   
!                                                                       
         UB = 1.0D10 
         IF (DABS (ALPHA (1) ) .GE.UB) THEN 
            ALPHAX = 0.0D0 
            ALPHAY = DSIGN (1.0D0, FN (1) - FN (0) ) 
         ELSE 
            ROOT = DSQRT (1.0D0 / (1.0D0 + ALPHA (1) * ALPHA (1) ) ) 
            ALPHAX = DSIGN (ROOT, XN (1) - XN (0) ) 
            ALPHAY = ALPHAX * ALPHA (1) 
         ENDIF 
         IF (DABS (BETA (1) ) .GE.UB) THEN 
            BETAX = 0.0D0 
            BETAY = DSIGN (1.0D0, FN (N) - FN (N - 1) ) 
         ELSE 
            ROOT = DSQRT (1.0D0 / (1.0D0 + BETA (1) * BETA (1) ) ) 
            BETAX = DSIGN (ROOT, XN (N) - XN (N - 1) ) 
            BETAY = BETAX * BETA (1) 
         ENDIF 
         CALL CFSP1D (N, T, XN, WX, ALPHAX, BETAX, 1, AX, BX, CX, DX,   &
         AUXF (1), AUXF (N + 1), AUXF (2 * N + 1), AUXF (3 * N + 1),    &
         AUXF (4 * N), AUXF (5 * N - 1), AUXF (6 * N - 2), IERR)        
         IF (IERR.NE.0) RETURN 
         CALL CFSP1D (N, T, FN, WY, ALPHAY, BETAY, MREP, AY, BY, CY, DY,&
         AUXF (1), AUXF (N + 1), AUXF (2 * N + 1), AUXF (3 * N + 1),    &
         AUXF (4 * N), AUXF (5 * N - 1), AUXF (6 * N - 2), IERR)        
      ELSE 
         IERR = - 1 
         IF (N.LT.6) RETURN 
         CALL CFSPPE (N, T, XN, WX, 1, AX, BX, CX, DX, AUXF (1),        &
         AUXF (N + 2), AUXF (2 * N + 3), AUXF (3 * N + 4), AUXF (4 * N +&
         5), AUXF (5 * N + 5), IERR)                                    
         IF (IERR.NE.0) RETURN 
         CALL CFSPPE (N, T, FN, WY, MREP, AY, BY, CY, DY, AUXF (1),     &
         AUXF (N + 2), AUXF (2 * N + 3), AUXF (3 * N + 4), AUXF (4 * N +&
         5), AUXF (5 * N + 5), IERR)                                    
      ENDIF 
      RETURN 
      END SUBROUTINE CFSPPA                         
