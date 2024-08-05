![             Splines}*)                                               
      SUBROUTINE ISPLPA (N, XN, FN, T, MT, IB, ALPHA, BETA, BX, CX, DX, &
      BY, CY, DY, DUMMY, IERR)                                          
!                                                                       
!*****************************************************************      
!                                                                *      
!  ISPLPA computes the coefficients BX(I), CX(I), DX(I), BY(I),  *      
!  CY(I), DY(I) for I=0,1,...,N-1 of a parametric cubic inter-   *      
!  polating spline for various end point conditions.             *      
!  The end point conditions are specified via the parameter IB.  *      
!  The parametric two-component function depends on the parameter*      
!  T(I) for  I=0,1,...,N and has two components SX and SY of the *      
!  the following form:                                           *      
!                                                                *      
!  SX := SX(T) = XN(I) + BX(I)(T-T(I)) + CX(I)(T-T(I))**2 +      *      
!                                      + DX(I)(T-T(I))**3        *      
!                                                                *      
!  SY := SY(T) = FN(I) + BY(I)(T-T(I)) + CY(I)(T-T(I))**2 +      *      
!                                      + DY(I)(T-T(I))**3        *      
!                                                                *      
!  for T in the interval [T(I),T(I+1)], I=0,1,...N-1.            *      
!                                                                *      
!  SX and SY each are nonparametric cubic splines.               *      
!                                                                *      
!                                                                *      
!  ASSUMPTIONS:    1.         N > 2                              *      
!  ============    2.      T(I) < T(I+1), I=0,1,...,N-1          *      
!                  3.     FN(0) = FN(N) , for IB = 4             *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N  :  Index of the last node                                  *      
!  XN :  vector XN(0:N); the nodes XN(I), I = 0,1,..,N           *      
!  FN :  vector FN(0:N); the functional values FN(I) = FN(XN(I)) *      
!  T  :  vector T(0:N); paramete value for XN(I) and FN(I)       *      
!  MT :  indicates origin of the curve parameter T(I):           *      
!        MT = 0 :  The user specifies the parameter values T(I), *      
!                  I=0,1,...,N                                   *      
!        MT = 1 :  Parameter values not prescribed. Their values *      
!                  are computed in SUBROUTINE PSPPV from the     *      
!                  chordal distances.                            *      
!        MT = 2 :  Parameter values not prescribed. Their values *      
!                  are computed in SUBROUTINE PSPPV from the     *      
!                  arc length.                                   *      
!                                                                *      
!  IB : describes the end point conditions                       *      
!       IB = 1: first derivative w.r.t. parameter is given       *      
!       IB = 2: second derivative w.r.t. parameter is given      *      
!       IB = 3: first derivative DY/DX is given                  *      
!       IB = 4: periodic spline                                  *      
!                                                                *      
!  ALPHA : ] 2-vectors ..(1:2)                                   *      
!  BETA  : ]                                                     *      
!              if IB = 1 : ] IBth derivative w.r.t. parameter    *      
!                 IB = 2 : ]       ALPHA(1)=SX(IB)(T(0))         *      
!                                  ALPHA(2)=SY(IB)(T(0))         *      
!                                  BETA(1) =SX(IB)(T(N))         *      
!                                  BETA(2) =SY(IB)(T(N))         *      
!              if IB = 3 :  first derivative DY/DX               *      
!                           ALPHA(1) = DY/DX (XN(0))             *      
!                           BETA(1)  = DY/DX (XN(N))             *      
!                           ALPHA(2) : not used                  *      
!                           BETA(2)  : not used                  *      
!               If the magnitude of ALPHA(1) or BETA(1) exceeds  *      
!               1.E10, then we compute the corresponding tangent *      
!               vector as follows:                               *      
!                . .                                             *      
!               (X,Y) = (0,DSIGN(1,FN(1)-FN(0))  (left end point)*      
!                . .                                             *      
!               (X,Y) = (0,DSIGN(1,FN(N)-FN(N-1))                *      
!                                               (right end point)*      
!                                                                *      
!              if IB = 4 : not used                              *      
!                                                                *      
!  (A natural parametric interpolating spline is obtained for    *      
!   IB=2 and ALPHA(1)=ALPHA(2)=BETA(1)=BETA(2)=0.0)              *      
!                                                                *      
!                                                                *      
!  AUXILIARY VBARIABLE:                                          *      
!  ====================                                          *      
!  DUMMY :  vector DUMMY(1:5*N+1)                                *      
!                                                                *      
!                                                                *      
!  AUSGABEPARAMETER:                                             *      
!  =================                                             *      
!  XN :  ] N+1-vectors ..(0:N);                                  *      
!  BX :  ] their elements 0, ..., N-1 describe the component     *      
!  CX :  ] function SX.                                          *      
!  DX :  ]                                                       *      
!                                                                *      
!  FN :  ] N+1-vectors ..(0:N);                                  *      
!  BY :  ] their elements 0, ..., N-1 describe the component     *      
!  CY :  ] function SY.                                          *      
!  DY :  ]                                                       *      
!          The elements BX(N), CX(N), DX(N), BY(N), CY(N) and    *      
!          DY(N) are used for auxiliary purposes.                *      
!  IERR :  error parameter                                       *      
!          =  0 :  All is ok                                     *      
!          = -1 :  N < 3                                         *      
!          = -2 :  IB < 1  or  IB > 4                            *      
!          = -3 :  Parameter values T(I) non monotonic:          *      
!                  T(I) > or = T(I+1) for some I=0,1,...,N-1     *      
!          = -4 :  For IB = 4 we have FN(0) not equal to FN(N)   *      
!          =  1 :  crash in TRDSY or CYTSY, system matrix numer- *      
!                  ically singular.                              *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Subroutines required: ISPL1D, ISPL2D, ISPLPE, PSPPV           *      
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
      DOUBLEPRECISION XN (0:N), FN (0:N), T (0:N), BX (0:N), CX (0:N),  &
      DX (0:N), BY (0:N), CY (0:N), DY (0:N), ALPHA (2), BETA (2),      &
      DUMMY (1:5 * N + 1)                                               
!                                                                       
!-----Check input data                                                  
!                                                                       
      IERR = - 1 
      IF (N.LT.3) RETURN 
      IF (IB.LT.1.OR.IB.GT.4) THEN 
         IERR = - 2 
         RETURN 
      ENDIF 
!                                                                       
!-----Find and check parameter values                                   
!                                                                       
      IF (MT.GT.0) THEN 
!                                                                       
!       Find parameter values in SUBROUTINE PSPPV                       
!                                                                       
         CALL PSPPV (N, XN, FN, T, MT, IERR) 
         IF (IERR.NE.0) THEN 
            IERR = - 3 
            RETURN 
         ENDIF 
      ELSE 
!                                                                       
!       Check the given parameter values                                
!                                                                       
         IERR = - 3 
         DO 20 I = 0, N - 1, 1 
            IF (T (I + 1) .LE.T (I) ) RETURN 
   20    END DO 
      ENDIF 
!                                                                       
!-----Compute the spline coefficients CX(I) and CY(I). Note that        
!     on second or subsequent calls, the system matrix need not         
!     be decomposed anew.                                               
!                                                                       
      IF (IB.EQ.1) THEN 
         CALL ISPL1D (N, T, XN, ALPHA (1), BETA (1), 1, BX, CX, DX,     &
         DUMMY (1), DUMMY (N + 1), DUMMY (2 * N), DUMMY (3 * N - 1),    &
         IERR)                                                          
         IF (IERR.NE.0) RETURN 
         CALL ISPL1D (N, T, FN, ALPHA (2), BETA (2), 2, BY, CY, DY,     &
         DUMMY (1), DUMMY (N + 1), DUMMY (2 * N), DUMMY (3 * N - 1),    &
         IERR)                                                          
      ELSEIF (IB.EQ.2) THEN 
         CALL ISPL2D (N, T, XN, ALPHA (1), BETA (1), 1, BX, CX, DX,     &
         DUMMY (1), DUMMY (N + 1), DUMMY (2 * N), DUMMY (3 * N - 1),    &
         IERR)                                                          
         IF (IERR.NE.0) RETURN 
         CALL ISPL2D (N, T, FN, ALPHA (2), BETA (2), 2, BY, CY, DY,     &
         DUMMY (1), DUMMY (N + 1), DUMMY (2 * N), DUMMY (3 * N - 1),    &
         IERR)                                                          
      ELSEIF (IB.EQ.3) THEN 
!                                                                       
!       Compute the tangent vectors using the first                     
!       derivatives of SX and SY                                        
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
         CALL ISPL1D (N, T, XN, ALPHAX, BETAX, 1, BX, CX, DX, DUMMY (1),&
         DUMMY (N + 1), DUMMY (2 * N), DUMMY (3 * N - 1), IERR)         
         IF (IERR.NE.0) RETURN 
         CALL ISPL1D (N, T, FN, ALPHAY, BETAY, 2, BY, CY, DY, DUMMY (1),&
         DUMMY (N + 1), DUMMY (2 * N), DUMMY (3 * N - 1), IERR)         
      ELSE 
         CALL ISPLPE (N, T, XN, 1, BX, CX, DX, DUMMY (1), DUMMY (N + 2),&
         DUMMY (2 * N + 2), DUMMY (3 * N + 2), DUMMY (4 * N + 2),       &
         IERR)                                                          
         IF (IERR.NE.0) RETURN 
         CALL ISPLPE (N, T, FN, 2, BY, CY, DY, DUMMY (1), DUMMY (N + 2),&
         DUMMY (2 * N + 2), DUMMY (3 * N + 2), DUMMY (4 * N + 2),       &
         IERR)                                                          
      ENDIF 
      RETURN 
      END SUBROUTINE ISPLPA                         
