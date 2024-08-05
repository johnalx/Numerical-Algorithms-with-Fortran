      SUBROUTINE DESABM (X, Y, FCT, N, XEND, H, HMAX, ABSERR, RELERR,   &
      IND, IERR, AUXF, LDF)                                             
!                                                                       
!*****************************************************************      
!                                                                *      
!  Numerical solution of a system of ordinary differential       *      
!  equations                                                     *      
!                                                                *      
!      Y' = F(X,Y) with initial condition Y(X0) = Y0             *      
!                                                                *      
!  using the predictor-corrector method by ADAMS-BASHFORTH-      *      
!  MOULTON. The starting values required are produced by the     *      
!  RUNGE-KUTTA method with the same order as the A-B-M method.   *      
!  An automatic step size control is used, that operates on the  *      
!  principle of doubling or halving of the step size with sub-   *      
!  sequent error estimation.                                     *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!                                                                *      
!  X     : starting x value.                                     *      
!  Y     : DOUBLE PRECISION vector Y(1:N), the initial values    *      
!          at location X.                                        *      
!  FCT   : SUBROUTINE, that evaluates the right hand side of the *      
!          differential equation. It has to be of the following  *      
!          form:                                                 *      
!                                                                *      
!             SUBROUTINE FCT(X,Y,N,F)                            *      
!             DOUBLE PRECISION Y(N), F(N)                        *      
!                  .                                             *      
!                  .                                             *      
!                  .                                             *      
!                                                                *      
!          In the calling program it must be declared as         *      
!          EXTERNAL. Its variables have the following meaning:   *      
!             X   independent variable.                          *      
!             Y   dependent variable.                            *      
!             F   function value of the right hand side of the   *      
!                 system of differential equations Y'=F(X,Y).    *      
!  N     : number of differential equations in the system.       *      
!  XEND  : location where the solution is desired; may be        *      
!          smaller than X.                                       *      
!  H     : step size for the next step, normally it is           *      
!          determined by the program.                            *      
!  HMAX  : maximum step size allowed, it has to be positive.     *      
!  ABSERR:)                                                      *      
!  RELERR:) error parameters, which have to be >= 0. A mixed test*      
!           is used:                                             *      
!              ABS (LOCAL ERROR) .LE. ABS(Y)*RELERR + ABSERR.    *      
!           Thus, if RELERR = 0 is chosen, this corresponds to a *      
!           test for the absolute error. If ABSERR = 0 is chosen *      
!           this corresponds to a test for the relative error.   *      
!           RELERR and ABSERR should be chosen to exceed ten     *      
!           times the machine constant. If this is not the case, *      
!           they are automatically set equal to this value.      *      
!  IND   : indicator, which must be set equal to one on the      *      
!          first call.                                           *      
!  IERR  : error parameter, on first call it has to be set equal *      
!          to  0.                                                *      
!  AUXF  : DOUBLE PRECISION auxiliary array AUXF(1:LDF,12),      *      
!          LDF .GE. N, where LDF is as defined in the calling    *      
!          program.                                              *      
!  LDF   : leading dimension of AUXF.                            *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!                                                                *      
!  X     : x value reached during last integration;              *      
!          normally X < XEND.                                    *      
!  Y     : value for the solution at location X.                 *      
!  H     : step size used last, should remain unchanged for the  *      
!          next step. If it is changed, IERR has to be reset to  *      
!          zero.                                                 *      
!  IND   : after the first call set equal to 2, should not be    *      
!          changed.                                              *      
!  IERR  : error parameter.                                      *      
!          = 0: everything o.k., after resetting of XEND         *      
!               DESABM can be called again.                      *      
!          = 1: everything o.k., after resetting of XEND         *      
!               DESABM can be called again.                      *      
!               the starting values for the A-B-M-method are     *      
!               known, the R-K set up procedure is omitted.      *      
!          = 2: after ICTMAX function evaluations the procedure  *      
!               did not reach XEND, a repeated call without any  *      
!               change of the parameters may be successful.      *      
!               (otherwise try to increase the error parameters).*      
!               ICTMAX is defined locally as DATA with a preset  *      
!               value of 500. Local changes can be made inside   *      
!               the program.                                     *      
!          = 3: the step size is less than eight times the       *      
!               machine constant at an x value. Before further   *      
!               calls, H and the error parameters need to be     *      
!               increased.                                       *      
!          = 4: ABSERR or RELERR is negative or both             *      
!               ABSERR = 0.0 and RELERR = 0.0                    *      
!          = 5: XEND equals X.                                   *      
!          = 6: HMAX is negative.                                *      
!                                                                *      
!  NOTE:                                                         *      
!  =====                                                         *      
!                                                                *      
!  The point XEND normally is not reached by the program, but due*      
!  to the step size control a point X with XEND - H < X < XEND.  *      
!  If the solution is required at XEND, DESABM has to be called  *      
!  again with H = XEND - X and IERR = 0.                         *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: MACHPD                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Jobst Hoffmann                                     *      
!  date     : 07.25.1988 (revised 11.11.1994)                    *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION X, XEND, H, HMAX, ABSERR, RELERR 
      DIMENSION Y (N), AUXF (LDF, 12) 
      DIMENSION ALPHA (4), P (4), BETA (4), BETDER (4), FERR (2) 
!                                                                       
      SAVE ALPHA, BETA, BETDER, P, FERR, FMACHP, X0 
!                                                                       
      DATA ICTMAX / 500 / 
!.                                                                      
!     Initialize the constants. This is performed                       
!     during the first call of the program only,                        
!     afterwards it is omitted.                                         
!.                                                                      
      IF (IND.NE.1) GOTO 10 
!.                                                                      
!     machine constant                                                  
!.                                                                      
      FMACHP = 1.0D0 
    5 FMACHP = 0.50D0 * FMACHP 
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 5 
      FMACHP = FMACHP * 2.0D0 
!.                                                                      
!     RUNGE-KUTTA coefficients                                          
!.                                                                      
      ALPHA (1) = 1.0D0 / 6.0D0 
      ALPHA (2) = 1.0D0 / 3.0D0 
      ALPHA (3) = 1.0D0 / 3.0D0 
      ALPHA (4) = 1.0D0 / 6.0D0 
!                                                                       
      P (1) = 0.0D0 
      P (2) = 0.50D0 
      P (3) = 0.50D0 
      P (4) = 1.0D0 
!.                                                                      
!     ADAMS-BASHFORTH coefficients                                      
!.                                                                      
      BETA (1) = - 9.0D0 
      BETA (2) = 37.0D0 
      BETA (3) = - 59.0D0 
      BETA (4) = 55.0D0 
!.                                                                      
!     ADAMS-MOULTON coefficients                                        
!.                                                                      
      BETDER (1) = 1.0D0 
      BETDER (2) = - 5.0D0 
      BETDER (3) = 19.0D0 
      BETDER (4) = 9.0D0 
!.                                                                      
!     factors for the error estimate (extrapolation)                    
!.                                                                      
      FERR (1) = 1.0D0 / 80.0D0 
      FERR (2) = - 19.0D0 / 270.0D0 
!                                                                       
      IND = 2 
!.                                                                      
!     Input error check, Initialize variables                           
!.                                                                      
   10 IF (ABSERR.LT.0.0.OR.RELERR.LT.0.0.OR.ABSERR + RELERR.LE.0.0)     &
      THEN                                                              
         IERR = 4 
         RETURN 
      ELSEIF (XEND.EQ.X) THEN 
         IERR = 5 
         RETURN 
      ELSEIF (HMAX.LE.0.0D0) THEN 
         IERR = 6 
         RETURN 
      ENDIF 
      IF (ABSERR.LT.10.0D0 * FMACHP) ABSERR = 10.0D0 * FMACHP 
      IF (RELERR.LT.10.0D0 * FMACHP) RELERR = 10.0D0 * FMACHP 
      ICOUNT = 0 
      IF (IERR.EQ.2) GOTO 170 
      IF (IERR.EQ.1) GOTO 210 
      IERR = 0 
      CALL FCT (X, Y (1), N, AUXF (1, 2) ) 
      ICOUNT = ICOUNT + 1 
!.                                                                      
!     store the starting step size                                      
!.                                                                      
      HLOC = H 
!.                                                                      
!     determine the current step size                                   
!.                                                                      
   20 IMETH = 1 
      IF (HMAX.LT.H) H = HMAX 
      ABSH = DABS (H) 
      H = DSIGN (DMIN1 (ABSH, DABS (XEND-X) / 3.0D0), XEND-X) 
      ABSH = DABS (H) 
      IF (ABSH.LE.8.0D0 * FMACHP * DABS (X) ) THEN 
         H = HLOC 
         IERR = 0 
         RETURN 
      ENDIF 
!.                                                                      
!     RUNGE-KUTTA steps with error estimates, in order                  
!     to determine the starting values for the multi-step               
!     method.                                                           
!     If  ISTART = 1:  step size  H.                                    
!     If  ISTART = 0:  step size 3H.                                    
!.                                                                      
   30 ISTART = 1 
   40 CONTINUE 
!.                                                                      
!     store the starting values                                         
!.                                                                      
      X0 = X 
      XDUMMY = X 
      DO 50 I = 1, N 
         AUXF (I, 6) = Y (I) 
         AUXF (I, 7) = Y (I) 
   50 END DO 
!.                                                                      
!     Perform the three steps                                           
!.                                                                      
      DO 120 JBEG = 1, 3 
         DO 60 I = 1, N 
            AUXF (I, 9) = 0.0D0 
   60    END DO 
         IRK = 1 
   70    CALL FCT (X0, AUXF (1, 7), N, AUXF (1, 8) ) 
         DO 80 I = 1, N 
            AUXF (I, 9) = AUXF (I, 9) + ALPHA (IRK) * AUXF (I, 8) 
   80    END DO 
         IF (IRK.EQ.4) GOTO 100 
         IRK = IRK + 1 
         X0 = XDUMMY + P (IRK) * H 
         DO 90 I = 1, N 
            AUXF (I, 7) = AUXF (I, 6) + P (IRK) * H * AUXF (I, 8) 
   90    END DO 
         GOTO 70 
  100    IF (ISTART.EQ.0) GOTO 130 
         XDUMMY = XDUMMY + H 
         X0 = XDUMMY 
         DO 110 I = 1, N 
            AUXF (I, 6) = AUXF (I, 6) + H * AUXF (I, 9) 
            AUXF (I, 7) = AUXF (I, 6) 
            IF (JBEG.EQ.3) AUXF (I, 11) = AUXF (I, 6) 
  110    END DO 
         CALL FCT (X0, AUXF (1, 6), N, AUXF (1, JBEG + 2) ) 
  120 END DO 
!.                                                                      
!     After three steps with step size H, perform one step              
!     with step size 3H                                                 
!.                                                                      
      ISTART = 1 - ISTART 
      HDUMMY = H 
      H = H + H + H 
      GOTO 40 
  130 ICOUNT = ICOUNT + 13 
      DO 140 I = 1, N 
         AUXF (I, 10) = AUXF (I, 6) + H * AUXF (I, 9) 
  140 END DO 
      X0 = X + H 
!.                                                                      
!     continue calculations with the original step size                 
!.                                                                      
      H = HDUMMY 
  150 IF (ICOUNT.GE.ICTMAX) THEN 
!.                                                                      
!        return with message about excessive function calls             
!.                                                                      
         X = X0 
         DO 160 I = 1, N 
            Y (I) = AUXF (I, 11) 
  160    END DO 
         IERR = 2 
         RETURN 
      ENDIF 
!.                                                                      
!     error estimate                                                    
!.                                                                      
  170 DMAX = 0.0D0 
      YMAX = 0.0D0 
      DO 180 I = 1, N 
         AUXF (I, 12) = FERR (IMETH) * (AUXF (I, 11) - AUXF (I, 10) ) 
         DMAX = DMAX1 (DMAX, DABS (AUXF (I, 12) ) ) 
         YMAX = DMAX1 (YMAX, DABS (AUXF (I, 11) ) ) 
  180 END DO 
      IF (DMAX.GE.RELERR * YMAX + ABSERR) THEN 
!.                                                                      
!        unsuccessful run, the step is repeated with half the step size 
!.                                                                      
         H = 0.5D0 * H 
         IF (DABS (H) .GT.8.0D0 * FMACHP * DABS (X0) ) GOTO 30 
         IERR = 3 
         RETURN 
      ENDIF 
!.                                                                      
!     the step was successful, calculation is continued with the        
!     old step size                                                     
!.                                                                      
      X = X0 
      DO 190 I = 1, N 
         AUXF (I, 11) = AUXF (I, 11) + AUXF (I, 12) 
         Y (I) = AUXF (I, 11) 
  190 END DO 
      CALL FCT (X0, Y (1), N, AUXF (1, 5) ) 
      ICOUNT = ICOUNT + 1 
      IF (DMAX.LE.0.02D0 * (RELERR * YMAX + ABSERR) ) THEN 
!.                                                                      
!        The precision obtained is too high, the step size is           
!        doubled, continue with the RUNGE-KUTTA starting procedure      
!.                                                                      
         H = 2.0D0 * H 
         HLOC = DMAX1 (H, HLOC) 
         DO 200 I = 1, N 
            AUXF (I, 2) = AUXF (I, 5) 
  200    END DO 
         IF (H.GT.0..AND.X0.LT.XEND.OR.H.LT.0..AND.X0.GT.XEND) GOTO 20 
         IERR = 0 
         RETURN 
      ENDIF 
      IF (H.GT.0..AND.X0 + H.GE.XEND.OR.H.LT.0..AND.X0 + H.LE.XEND)     &
      THEN                                                              
!.                                                                      
!        We are at the end of the interval. For further calls IERR is se
!        1, so that on repeated calls the calculations can be continued 
!        with ADAMS-BASHFORTH-MOULTON steps                             
!.                                                                      
         IF (ABSH.LE.8.0D0 * FMACHP * DABS (X0) ) THEN 
            H = HLOC 
            IERR = 0 
            RETURN 
         ENDIF 
         IERR = 1 
         RETURN 
      ENDIF 
!.                                                                      
!     start of the ADAMS-BASHFORTH-MOULTON steps                        
!.                                                                      
  210 IMETH = 2 
      DO 230 J = 1, 4 
         DO 220 I = 1, N 
            AUXF (I, J) = AUXF (I, J + 1) 
  220    END DO 
  230 END DO 
!.                                                                      
!     predictor-step                                                    
!.                                                                      
      DO 240 I = 1, N 
         AUXF (I, 8) = 0.0D0 
  240 END DO 
      DO 260 J = 1, 4 
         DO 250 I = 1, N 
            AUXF (I, 8) = AUXF (I, 8) + BETA (J) * AUXF (I, J) 
  250    END DO 
  260 END DO 
      DO 270 I = 1, N 
         AUXF (I, 10) = AUXF (I, 11) + H * AUXF (I, 8) / 24.0D0 
  270 END DO 
!.                                                                      
!     new node                                                          
!.                                                                      
      X0 = X0 + H 
      CALL FCT (X0, AUXF (1, 10), N, AUXF (1, 5) ) 
!.                                                                      
!     corrector step                                                    
!.                                                                      
      DO 280 I = 1, N 
         AUXF (I, 8) = 0.0D0 
  280 END DO 
      DO 300 J = 1, 4 
         DO 290 I = 1, N 
            AUXF (I, 8) = AUXF (I, 8) + BETDER (J) * AUXF (I, J + 1) 
  290    END DO 
  300 END DO 
      DO 310 I = 1, N 
         AUXF (I, 11) = AUXF (I, 11) + H * AUXF (I, 8) / 24.0D0 
  310 END DO 
!.                                                                      
!     function value of the corrector                                   
!.                                                                      
      CALL FCT (X0, AUXF (1, 11), N, AUXF (1, 5) ) 
      ICOUNT = ICOUNT + 2 
      GOTO 150 
      END SUBROUTINE DESABM                         
