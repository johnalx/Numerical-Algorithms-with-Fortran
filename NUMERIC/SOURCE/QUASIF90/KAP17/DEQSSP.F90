![KA{P 17}{Initial Value Problems}                                      
![        {Initial Value Problems for Ordinary Differential             
![         Equations}*)                                                 
      SUBROUTINE DEQSSP (X, Y, FCT, XEND, H, HMAX, ABSERR, RELERR, IND, &
      IERR, AUXF)                                                       
!                                                                       
!*****************************************************************      
!                                                                *      
!  Numerical solution of an ordinary differential equation       *      
!                                                                *      
!     Y' = F(X,Y) with given initial value  Y(X0) = Y0           *      
!                                                                *      
!  by a user specified one-step method.                          *      
!  Any one of the following methods may be chosen:               *      
!  1) EULER-CAUCHY polygonal method,                             *      
!  2) method of HEUN, and                                        *      
!  3) classical RUNGE-KUTTA method.                              *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!                                                                *      
!  X     : starting value for X.                                 *      
!  Y     : initial value of the solution Y at X.                 *      
!  FCT   : SUBROUTINE that evaluates the right hand side of the  *      
!          differential equation. It has to be formatted as      *      
!          follows:                                               *     
!                                                                *      
!             SUBROUTINE FCT(X,Y,F)                              *      
!             DOUBLE PRECISION X,Y,F                             *      
!                  .                                             *      
!                  .                                             *      
!                  .                                             *      
!                                                                *      
!          In the calling program it must be declared as         *      
!          EXTERNAL. Its variables have the following meaning:   *      
!             X   independent variable X.                        *      
!             Y   dependent variable Y.                          *      
!             F   function value of the right hand side F(X,Y)   *      
!                 of the differential equation Y'=F(X,Y).        *      
!  XEND  : value for X where the solution is desired; it may be  *      
!          to the left of the starting value for X.              *      
!  H     : step size for the next step. It is usually determined *      
!          by the program.                                       *      
!  HMAX  : maximum step size allowed, HMAX must be positive.     *      
!  ABSERR: )                                                     *      
!  RELERR: ) error parameters, which must be >= 0. The subroutine*      
!            performs a mixed test:                              *      
!              ABS (LOCAL ERROR) .LE. ABS(Y)*RELERR + ABSERR.    *      
!           Thus, if RELERR is chosen as zero, this tests for the*      
!           absolute error. If ABSERR is chosen as zero, this    *      
!           tests for the relative error.                        *      
!           RELERR and ABSERR should be chosen at least ten times*      
!           larger than the machine constant. If this is not     *      
!           the case the program defaults to these settings.     *      
!  IND   : INTEGER vector IND(1:3).                              *      
!          IND(1) indicates whether DEQSSP has already been      *      
!          called:                                               *      
!          = 0: DEQSSP was not been called during this program   *      
!               run.                                             *      
!          = 1: DEQSSP has been called during this run.          *      
!          IND(2) indicates the method used:                     *      
!          = 1: EULER-CAUCHY poygonal method.                    *      
!          = 2: method of HEUN.                                  *      
!          = 3: classical RUNGE-KUTTA method.                    *      
!          IND(3) indicates whether integration is acceptable    *      
!          beyond XEND (this can be used e.g. for discontinuities*      
!          of the right hand side F of the differential equation)*      
!          = 0: no restrictions.                                 *      
!          = 1: the differential equation may only be integrated *      
!               up to XEND.                                      *      
!  AUXF : REAL auxiliary vector AUXF(0:4,3)                      *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!                                                                *      
!  X     : X value for which solution is desired (should equal   *      
!          input XEND)                                           *      
!  Y     : value for the solution at X.                          *      
!  H     : terminal step size used. It should remain unchanged   *      
!          for the next call.                                    *      
!  IND(1): after the first call set to 2, it should not be       *      
!          changed.                                              *      
!  IERR  : error parameter:                                      *      
!          = 0: everything o.k. After re-setting of XEND, DEQSSP *      
!               can be called again.                             *      
!          = 1: everything o.k. After re-setting of XEND, DEQSSP *      
!               can be called again. The solution has already    *      
!               been determined up to the x value AUXF(IPTR,1)   *      
!               with AUXF(IPTR,2) as its y value. For XEND .LT.  *      
!               AUXF(IPTR,1) further solution values will be     *      
!               determined by interpolation.                     *      
!          = 2: after 500 function evaluations the procedure did *      
!               not reach XEND. A repeated call without any      *      
!               change of the parameters may be successful.      *      
!               (otherwise, try to increase the error parameters)*      
!          = 3: the step size is less than eight times the       *      
!               machine constant at an x value. To continue with *      
!               another call, H and the error parameters must be *      
!               increased.                                       *      
!          = 4: ABSERR or RELERR is negative or both             *      
!               ABSERR = 0.0 and RELERR = 0.0.                   *      
!          = 5: XEND is equal to X.                              *      
!          = 6: HMAX is negative.                                *      
!          = 7: number for method is erroneous.                  *      
!          = 8: IND(3) is erroneous.                             *      
!                                                                *      
!  NOTES:                                                        *      
!  ======                                                        *      
!                                                                *      
!  1) The program uses step size control by doubling or halving  *      
!     the step size. Thus, it independently finds the step sizes *      
!     for which the required error bounds can be maintained.     *      
!  2) If the distance remaining towards XEND is smaller than     *      
!     the current, automatically determined step size, the       *      
!     solution at XEND is determined by interpolation.           *      
!     The interpolation is effected by use of a circular buffer. *      
!     This buffer has to contain a sufficient number of entries  *      
!     before it can be used to interpolate with a sufficiently   *      
!     high order. Thus it may happen that the solution is        *      
!     determined much farther than to XEND in the first call.    *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: MACHPD                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Jobst Hoffmann                                     *      
!  date     : 04.21.1990                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      INTEGER IND (3), IBUF (3), ICTFV (3), ICTV (3) 
!                                                                       
      DIMENSION ERRHAT (3), AUXF (0:4, 3), ERREST (3) 
!                                                                       
      LOGICAL ERRV, FULL 
!                                                                       
      SAVE IPTR, FULL, LEN, FMACHP 
!                                                                       
      DATA IBEG / 1 / METH / 2 / IEND / 3 / 
      DATA ICTMAX / 500 / IBUF / 2, 3, 5 / ICTFV / 1, 3, 7 / ICTV / 2,  &
      5, 11 /                                                           
      DATA ERRHAT / 0.2D0, 0.1D0, 0.02D0 / ERREST / 1.0D0, 3.0D0,       &
      15.0D0 /                                                          
!                                                                       
      IF (IND (IBEG) .EQ.1) GOTO 10 
!                                                                       
!     determine the machine constant                                    
!.                                                                      
      FMACHP = 1.0D0 
    5 FMACHP = 0.5D0 * FMACHP 
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 5 
      FMACHP = FMACHP * 2.0D0 
!                                                                       
      IERR = 0 
      IND (IBEG) = 1 
!.                                                                      
!     determine the length of the circular buffer needed for the        
!     interpolation and initialize pointer to the first element         
!     of this buffer                                                    
!.                                                                      
      LEN = IBUF (IND (METH) ) 
      IPTR = 0 
      FULL = .FALSE. 
   10 IF (IERR.EQ.1) THEN 
!.                                                                      
!        the solution is known up to the point (X,Y) =                  
!        (AUXF(IPTR,1), AUXF(IPTR,2)). Further output                   
!        values may be produced by interpolation.                       
!.                                                                      
         X0 = X 
         Y0 = Y 
         X1 = AUXF (IPTR, 1) 
         Y1 = AUXF (IPTR, 2) 
         GOTO 170 
      ELSEIF (IERR.EQ.2) THEN 
!.                                                                      
!        repeat call after excessive number of function calls           
!.                                                                      
         ICOUNT = 0 
         GOTO 120 
      ENDIF 
!.                                                                      
!     Input data error check, initialize variables                      
!.                                                                      
      IF (ABSERR.LE.0.0D0.AND.RELERR.LE.0.0D0) THEN 
         IERR = 4 
         RETURN 
      ELSEIF (XEND.EQ.X) THEN 
         IERR = 5 
         RETURN 
      ELSEIF (HMAX.LE.0.0D0) THEN 
         IERR = 6 
         RETURN 
      ELSEIF (IND (METH) .LT.1.OR.IND (METH) .GT.3) THEN 
         IERR = 7 
         RETURN 
      ELSEIF (IND (IEND) .LT.0.OR.IND (IEND) .GT.1) THEN 
         IERR = 8 
         RETURN 
      ENDIF 
      ICOUNT = 0 
!.                                                                      
!     start computation, store the starting values in AUXF              
!.                                                                      
      HLOC = H 
      X0 = X 
      Y0 = Y 
      AUXF (IPTR, 1) = X 
      AUXF (IPTR, 2) = Y 
!.                                                                      
!     determine step size                                               
!.                                                                      
      IF (DABS (H) .GT.HMAX) H = HMAX 
      IF (IND (IEND) .EQ.0) THEN 
         H = DSIGN (H, XEND-X) 
      ELSE 
         H = DSIGN (DMIN1 (DABS (H), DABS (XEND-X) * 0.5D0), XEND-X) 
      ENDIF 
      IF (DABS (H) .LE.8.0D0 * FMACHP * DABS (X) ) THEN 
         H = HLOC 
         IERR = 3 
         RETURN 
      ENDIF 
!.                                                                      
!     switch between the three allowed methods                          
!.                                                                      
      GOTO (20, 50, 80) IND (METH) 
!.                                                                      
!     *** EULER-CAUCHY  polygonal method ***                            
!.                                                                      
   20 CONTINUE 
   30 CALL FCT (X0, Y0, F0) 
      XH = X0 + H 
      X1 = XH + H 
      Y2 = Y0 + 2.0D0 * H * F0 
   40 YH = Y0 + H * F0 
      CALL FCT (XH, YH, F1) 
      Y1 = YH + H * F1 
      GOTO 110 
!.                                                                      
!     *** method of HEUN ***                                            
!.                                                                      
   50 CONTINUE 
   60 CALL FCT (X0, Y0, F0) 
      XH = X0 + H 
      X1 = XH + H 
      Y2 = Y0 + 2.0D0 * H * F0 
      CALL FCT (X1, Y2, F2) 
      Y2 = Y0 + H * (F0 + F2) 
   70 YH = Y0 + H * F0 
      CALL FCT (XH, YH, FH) 
      YH = Y0 + 0.5D0 * H * (F0 + FH) 
      CALL FCT (XH, YH, FH) 
      Y1 = YH + H * FH 
      CALL FCT (X1, Y1, F1) 
      Y1 = YH + 0.5D0 * H * (FH + F1) 
      GOTO 110 
!.                                                                      
!     *** classical RUNGE-KUTTA method ***                              
!.                                                                      
   80 CONTINUE 
   90 CALL FCT (X0, Y0, F0) 
      XH = X0 + H 
      X1 = XH + H 
      RK12 = 2.0D0 * H * F0 
      CALL FCT (XH, Y0 + 0.5D0 * RK12, F2) 
      RK22 = 2.0D0 * H * F2 
      CALL FCT (XH, Y0 + 0.5D0 * RK22, F2) 
      RK32 = 2.0D0 * H * F2 
      CALL FCT (X1, Y0 + RK32, F2) 
      RK42 = 2.0D0 * H * F2 
      Y2 = Y0 + (RK12 + 2.0D0 * (RK22 + RK32) + RK42) / 6.0D0 
  100 RK11 = H * F0 
      CALL FCT (X0 + 0.5D0 * H, Y0 + 0.5D0 * RK11, FH) 
      RK21 = H * FH 
      CALL FCT (X0 + 0.5D0 * H, Y0 + 0.5D0 * RK21, FH) 
      RK31 = H * FH 
      CALL FCT (XH, Y0 + RK31, FH) 
      RK41 = H * FH 
      YH = Y0 + (RK11 + 2.0D0 * (RK21 + RK31) + RK41) / 6.0D0 
      CALL FCT (XH, YH, F1) 
      RK11 = H * F1 
      CALL FCT (XH + 0.5D0 * H, YH + 0.5D0 * RK11, F1) 
      RK21 = H * F1 
      CALL FCT (XH + 0.5D0 * H, YH + 0.5D0 * RK21, F1) 
      RK31 = H * F1 
      CALL FCT (X1, YH + RK31, F1) 
      RK41 = H * F1 
      Y1 = YH + (RK11 + 2.0D0 * (RK21 + RK31) + RK41) / 6.0D0 
!.                                                                      
!     count the number of function evaluations                          
!.                                                                      
  110 ICOUNT = ICOUNT + ICTV (IND (METH) ) 
      IF (ERRV) THEN 
         ERRV = .FALSE. 
         ICOUNT = ICOUNT + ICTFV (IND (METH) ) 
      ENDIF 
      IF (ICOUNT.GT.ICTMAX) THEN 
         IERR = 2 
         X = AUXF (IPTR, 1) 
         Y = AUXF (IPTR, 2) 
         RETURN 
      ENDIF 
!.                                                                      
!     error estimate                                                    
!.                                                                      
  120 D = DABS (Y1 - Y2) / ERREST (IND (METH) ) 
      IF (D.GT.ABSERR + DABS (Y1) * RELERR) THEN 
         ERRV = .TRUE. 
!.                                                                      
!        computation was unsuccessful, recompute with                   
!        half the step size                                             
!.                                                                      
         H = 0.5D0 * H 
         IF (DABS (H) .LE.8.0D0 * FMACHP * DABS (X) ) THEN 
            IERR = 3 
            X = AUXF (IPTR, 1) 
            Y = AUXF (IPTR, 2) 
            RETURN 
         ENDIF 
         X1 = XH 
         XH = X0 + H 
         Y2 = YH 
         GOTO (40, 70, 100) IND (METH) 
      ELSE 
!.                                                                      
!        computation was successful, compute an improved                
!        estimate ...                                                   
!.                                                                      
         GOTO (130, 140, 150) IND (METH) 
!.                                                                      
!        ... for EULER-CAUCHY ...                                       
!.                                                                      
  130    Y1 = 2.0D0 * Y1 - Y2 
         GOTO 160 
!.                                                                      
!        ... for HEUN ...                                               
!.                                                                      
  140    Y1 = (4.0D0 * Y1 - Y2) / 3.0D0 
         GOTO 160 
!.                                                                      
!        ... for RUNGE-KUTTA ...                                        
!.                                                                      
  150    Y1 = (16.0D0 * Y1 - Y2) / 15.0D0 
!.                                                                      
!        ... and the new values are stored in the buffer AUXF           
!.                                                                      
  160    IF (IPTR + 2.EQ.LEN) FULL = .TRUE. 
         IPTR = MOD (IPTR + 1, LEN) 
         AUXF (IPTR, 1) = X1 
         AUXF (IPTR, 2) = Y1 
      ENDIF 
!.                                                                      
!     if the solution varies only slightly, the next step may be        
!     executed with double the previous step size, provided this        
!     does not exceed the maximal allowable step size                   
!.                                                                      
      IF (D.LT.ERRHAT (IND (METH) ) * (ABSERR + DABS (Y1) * RELERR) ) H &
      = DSIGN (DMIN1 (DABS (2.0D0 * H), HMAX), H)                       
!.                                                                      
!     if XEND was not been reached, another step will be                
!     executed. Otherwise ...                                           
!.                                                                      
  170 IF ( (X1 - XEND) * DSIGN (1.0D0, H) ) 180, 190, 200 
  180 IF (IND (IEND) .EQ.1) THEN 
         HLOC = H 
         H = DSIGN (DMIN1 (DABS (H), DABS (XEND-X1) * 0.5D0), XEND-X1) 
      ENDIF 
      X0 = X1 
      Y0 = Y1 
      GOTO (30, 60, 90) IND (METH) 
!.                                                                      
!     ... return to calling program ...                                 
!.                                                                      
  190 X = X1 
      Y = Y1 
      IERR = 0 
      RETURN 
!.                                                                      
!     ... or the desired solution is obtained by interpolation          
!     in case that the steps went beyond XEND.                          
!     If there are not sufficiently many values in the buffer           
!     we perform another step.                                          
!.                                                                      
  200 IF (.NOT.FULL) THEN 
         X0 = X1 
         Y0 = Y1 
         GOTO (30, 60, 90) IND (METH) 
      ENDIF 
      DO 210 I = 0, LEN - 1 
         AUXF (I, 3) = AUXF (I, 2) 
  210 END DO 
      DO 230 I = 1, LEN - 1 
         DO 220 K = LEN - 1, I, - 1 
            XDIF = AUXF (K, 1) - AUXF (K - I, 1) 
            AUXF (K, 3) = (AUXF (K, 3) - AUXF (K - 1, 3) ) / XDIF 
  220    END DO 
  230 END DO 
      Y = AUXF (LEN - 1, 3) 
      DO 240 I = LEN - 2, 0, - 1 
         Y = Y * (XEND-AUXF (I, 1) ) + AUXF (I, 3) 
  240 END DO 
      X = XEND 
      IERR = 1 
      RETURN 
!                                                                       
      END SUBROUTINE DEQSSP                         
