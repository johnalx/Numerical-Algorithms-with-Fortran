      SUBROUTINE FRKFSY (A, DA, N, Y, DES, H, HMX, ABSERR, RELERR, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  A system of ordinary differential equations of 1st order is   *      
!  integrated by applying the RUNGE-KUTTA-FEHLBERG method        *      
!  [ of order O(H**5) ] with estimates for the local error and   *      
!  step size control.                                            *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  A      : starting value for the integration interval          *      
!  DA     : length of the integration interval;                  *      
!           DA may be < 0.0 if we want to integrate to the left. *      
!  N      : number of equations; N < 11.                         *      
!  Y      : vector Y(1:N); the initial values at A               *      
!  DES    : SUBROUTINE, that describes the system of differential*      
!           equations, given in the following form:              *      
!                      SUBROUTINE  DES(X,Y,YS)                   *      
!                      X : independent variable                  *      
!                      Y : vector of dependent variables         *      
!                      YS: vector YS(I)=DY(I)/DX of derivatives  *      
!                          at X, I=1,...,N                       *      
!                      Y and YS are dimensioned as DOUBLE        *      
!                      PRECISION Y(1), YS(1), however, they may  *      
!                      be to be used as a vector of length N.    *      
!           example :  SUBROUTINE DES(X,Y,YS)                    *      
!                      DOUBLE PRECISION Y(1),YS(1)               *      
!                      YS(1)=Y(1)                                *      
!                      YS(2)=-Y(2)                               *      
!                      RETURN                                    *      
!                      END                                       *      
!  H      : initial step size; if H is provided unrealistically, *      
!           H is modified internally; H may be negative if       *      
!           DA < 0.0.                                            *      
!  HMX    : upper bound for the step size magnitude used during  *      
!           calculation. HMX > 0.0                               *      
!  ABSERR :] bounds for the acceptable local error, relative to  *      
!  RELERR :] the current step size. If the following holds for   *      
!         :] each component of the computed solution Y(I)        *      
!                  ABS ( estimate of the local error) .LE.       *      
!                      ABS(H)*(RELERR*ABS(Y(I))+ABSERR),         *      
!           then the solution is accepted in the current step.   *      
!           If ABSERR = 0.0, we test for the relative error;     *      
!           If RELERR = 0.0, we test for the absolute error.     *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  A      : last x value for which a solution was successfully   *      
!           determined. Normally the following will hold:        *      
!           A on output = A on input + DA.                       *      
!  Y      : computed solution vector at  A on output             *      
!  H      : optimal step size, which was used for the last step. *      
!  IERR   : = 1, everything o.k.; solution found at A + DA.      *      
!           = 2, after 3000 calls of SUBROUTINE DES we stop with-*      
!                out having reached the endpoint A+DA. If com-   *      
!                putations are to be continued, call FRKFSY again*      
!                with unchanged parameters.                      *      
!           = 3, false input data; i.e.                          *      
!                ABSERR.LT.0.0     or    RELERR.LT.0.0     or    *      
!                ABSERR + RELERR = 0.0  or  HMX.LE.0.0: Return.  *      
!           = 4, the optimal step size cannot be achieved for the*      
!                computer.     RETURN                            *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required : MACHPD                                 *      
!                                                                *      
!                                                                *      
!  sources : SHAMPINE/ALLEN, see [SHAM73].                       *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Richard Reuter                                     *      
!  date     : 02.09.1983                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION YT (10), T (10), R (10), K1 (10), K2 (10) 
      DOUBLEPRECISION K3 (10), K4 (10), K5 (10), K6 (10) 
      DOUBLEPRECISION Y (N) 
!                                                                       
!     determine machine constant                                        
!                                                                       
      FMACHP = 1.0D0 
    2 FMACHP = 0.5D0 * FMACHP 
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 2 
      FMACHP = FMACHP * 2.0D0 
!                                                                       
!     check the input data                                              
!                                                                       
      IERR = 3 
      IF (RELERR.LT.0.0D0.OR.ABSERR.LT.0.0D0.OR.RELERR +                &
      ABSERR.EQ.0.0D0.OR.HMX.LE.0.0D0) RETURN                           
!                                                                       
      IERR = 4 
      B = A + DA 
      IF (DABS (DA) .LE.13.0D0 * FMACHP * DMAX1 (DABS (A), DABS (B) ) ) &
      RETURN                                                            
!                                                                       
      HMAX = DMIN1 (HMX, DABS (DA) ) 
      IF (DABS (H) .LE.13.0D0 * FMACHP * DABS (A) ) H = HMAX 
!                                                                       
!     Initialize counter for calls of SUBROUTINE DES                    
!                                                                       
      LFD = 0 
      IAD = 0 
!                                                                       
!     H is bounded by HMAX and is chosen so that the                    
!     endpoint B is reached, if possible.                               
!                                                                       
    3 H = DSIGN (DMIN1 (DABS (H), HMAX), DA) 
      IF (DABS (B - A) .LE.1.25D0 * DABS (H) ) THEN 
         HF = H 
!                                                                       
!        if IAD=1 and H=B-A acceptable, we stop after                   
!        the next integration step.                                     
!                                                                       
         IAD = 1 
         H = B - A 
      ENDIF 
!                                                                       
!     an integration step is executed                                   
!                                                                       
      CALL DES (A, Y, K1) 
      LFD = LFD+1 
    5 CONTINUE 
      X = 0.25D0 * H 
      DO 6 I = 1, N 
         YT (I) = Y (I) + X * K1 (I) 
    6 END DO 
      X = A + X 
      CALL DES (X, YT, K2) 
      DO 7 I = 1, N 
         YT (I) = Y (I) + H * (K1 (I) * (3.0D0 / 32.0D0) + K2 (I)       &
         * (9.0D0 / 32.0D0) )                                           
    7 END DO 
      X = A + H * (3.0D0 / 8.0D0) 
      CALL DES (X, YT, K3) 
      DO 8 I = 1, N 
         YT (I) = Y (I) + H * (K1 (I) * (1932.0D0 / 2197.0D0) - K2 (I)  &
         * (7200.0D0 / 2197.0D0) + K3 (I) * (7296.0D0 / 2197.0D0) )     
    8 END DO 
      X = A + H * (12.0D0 / 13.0D0) 
      CALL DES (X, YT, K4) 
      DO 9 I = 1, N 
         YT (I) = Y (I) + H * (K1 (I) * (439.0D0 / 216.0D0) - 8.0D0 *   &
         K2 (I) + K3 (I) * (3680.0D0 / 513.0D0) - K4 (I) * (845.0D0 /   &
         4104.0D0) )                                                    
    9 END DO 
      X = A + H 
      CALL DES (X, YT, K5) 
      DO 10 I = 1, N 
         YT (I) = Y (I) + H * ( - K1 (I) * (8.0D0 / 27.0D0) + 2.0D0 *   &
         K2 (I) - K3 (I) * (3544.0D0 / 2565.0D0) + K4 (I) * (1859.0D0 / &
         4104.0D0) - K5 (I) * (11.0D0 / 40.0D0) )                       
   10 END DO 
      X = A + 0.5D0 * H 
      CALL DES (X, YT, K6) 
      DO 11 I = 1, N 
         T (I) = K1 (I) * (25.0D0 / 216.0D0) + K3 (I) * (1408.0D0 /     &
         2565.0D0) + K4 (I) * (2197.0D0 / 4104.0D0) - K5 (I) * 0.20D0   
         YT (I) = Y (I) + H * T (I) 
   11 END DO 
!                                                                       
!     YT(I) now represents the latest result of this pass.              
!     Determine R(I), the estimate of the local                         
!     error, relative to the current step size.                         
!                                                                       
      DO 12 I = 1, N 
         R (I) = K1 (I) / 360.0D0 - K3 (I) * (128.0D0 / 4275.0D0)       &
         - K4 (I) * (2197.0D0 / 75240.0D0) + K5 (I) / 50.0D0 + K6 (I)   &
         * (2.0D0 / 55.0D0)                                             
   12 END DO 
!                                                                       
!     Check accuracy                                                    
!                                                                       
      QUOT = 0.0D0 
      DO 13 I = 1, N 
         TR = DABS (R (I) ) / (RELERR * DABS (YT (I) ) + ABSERR) 
         QUOT = DMAX1 (QUOT, TR) 
   13 END DO 
!                                                                       
!     If  QOUOT.LE.1.0   ==> integration step is accepted               
!                                                                       
      IF (QUOT.LE.1.0D0) THEN 
!                                                                       
!        result is accepted                                             
!                                                                       
         DO 14 I = 1, N 
            Y (I) = YT (I) 
   14    END DO 
         A = A + H 
!                                                                       
!        if A=B ,  RETURN                                               
!                                                                       
         IF (IAD.EQ.1) THEN 
            IERR = 1 
            H = HF 
            RETURN 
         ENDIF 
!                                                                       
!        prepare next step                                              
!                                                                       
         QUOT = DMAX1 (QUOT, 6.5536D-4) 
      ENDIF 
      QUOT = DMIN1 (QUOT, 4096.0D0) 
      H = 0.8D0 * H / DSQRT (DSQRT (QUOT) ) 
!                                                                       
!     We just achieved that H was increased by at most a factor of 5,   
!     or alternatively, that it was decreased by a factor of 10         
!     at most                                                           
!                                                                       
      IF (DABS (H) .LE.13.0D0 * FMACHP * DABS (A) ) THEN 
         IERR = 4 
         RETURN 
      ENDIF 
      LFD = LFD+5 
      IF (LFD.GE.2995) THEN 
         IERR = 2 
         RETURN 
      ENDIF 
      IF (QUOT.LE.1.0D0) THEN 
!                                                                       
!        the step was successful. Continue with another step.           
!                                                                       
         GOTO 3 
      ELSE 
!                                                                       
!        the step is repeated for a smaller H.                          
!                                                                       
         IAD = 0 
         GOTO 5 
      ENDIF 
      END SUBROUTINE FRKFSY                         
