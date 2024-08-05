![  {Bulirsch--Stoer--Gragg Extrapolation}                              
![  {Bulirsch--Stoer--Gragg Extrapolation}*)                            
      SUBROUTINE DESEXT (X, Y, FCT, N, XEND, H, HMAX, ABSERR, RELERR,   &
      IERR, AUXF, LDF)                                                  
!                                                                       
!*****************************************************************      
!                                                                *      
!  Numerical solution of a system of ordinary differential       *      
!  equations                                                     *      
!                                                                *      
!      Y' = F(X,Y)  with  given initial condition Y(X0) = Y0     *      
!                                                                *      
!  using the extrapolations method of BULIRSCH-STOER.            *      
!  The maximal order of the extrapolation is determined depending*      
!  the machine constant; the step size control lollows [HALL76]. *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!                                                                *      
!  X     : initial x value.                                      *      
!  Y     : REAL vector Y(1:N), the initial value at X.           *      
!  FCT   : SUBROUTINE, that evaluates the right hand side of the *      
!          system of differential equations. It has to be of the *      
!          following form:                                       *      
!                                                                *      
!             SUBROUTINE FCT(X,Y,F,N)                            *      
!             DOUBLE PRECISION X, Y(N), F(N)                     *      
!                  .                                             *      
!                  .                                             *      
!                  .                                             *      
!                                                                *      
!          In the calling program it has to be defined as        *      
!          EXTERNAL. In it the variables have the following      *      
!          meaning:                                              *      
!             X   independent variable.                          *      
!             Y   dependent variable.                            *      
!             F   function value of the right hand side of the   *      
!                 system of differential equations Y'=F(X,Y)     *      
!             N   number of differential equations.              *      
!  N     : number of differential equations.                     *      
!  XEND  : location where the solution is desired; may be        *      
!          smaller than X.                                       *      
!  H     : step size for the next step. Normally it is deter-    *      
!          mined by the program.                                 *      
!  HMAX  : maximally allowed step size, it has to be positivs.   *      
!  ABSERR:)                                                      *      
!  RELERR:) error parameters, which have to be >= 0. A mixed     *      
!           test is performed:                                   *      
!              ABS (LOCAL ERROR) .LE. ABS(Y)*RELERR + ABSERR.    *      
!           Thus if RELERR = 0 is chosen, this corresponds to a  *      
!           test for the absolute error. If ABSERR = 0 is chosen,*      
!           this corresponds to a test for the relative error.   *      
!           RELERR and ABSERR should be chosen larger than ten   *      
!           times the machine constant. If this is not the case, *      
!           they are automatically set equal to this value.      *      
!  IERR  : error parameter, on first call it has to be set equal *      
!          to 0.                                                 *      
!  AUXF  : REAL auxiliary array AUXF(1:LDF,12), LDF .GE. N, with *      
!          LDF as defined in the calling program.                *      
!  LDF   : leading dimension of AUXF.                            *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!                                                                *      
!  X     : x value that was reached after last during integration*      
!          Normally  X = XEND.                                   *      
!  Y     : solution value at X.                                  *      
!  H     : step size used last, this should remain unchanged for *      
!          the next step. If it is changed IERR has to be set to *      
!          equal zero.                                           *      
!  IERR  : error parameter.                                      *      
!          = 0: everything o.k., after resetting of XEND,        *      
!               DESEXT can be called again.                      *      
!          = 1: after 700 function evaluations the procedure did *      
!               not reach XEND. A repeated call without any      *      
!               change of the parameters may be successful.      *      
!               (otherwise try with increased error parameters)  *      
!          = 2: the step size is below four times the machine    *      
!               constant for an x value. Before further calls H  *      
!               and the error parameters must be increased.      *      
!          = 3: ABSERR or RELERR is negative or both             *      
!               ABSERR = 0.0 and RELERR = 0.0.                   *      
!          = 4: XEND is equal to X.                              *      
!          = 5: HMAX is negative.                                *      
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
      DIMENSION Y (N), AUXF (LDF, 12), BUFOL (12) 
      LOGICAL IND 
!                                                                       
      SAVE IEXTMX, FMACHP, IND 
!                                                                       
      DATA IND / .TRUE. /, ICTMAX / 700 / 
      DATA BUFOL / 2.0D0, 4.0D0, 6.0D0, 8.0D0, 12.0D0, 16.0D0, 24.0D0,  &
      32.0D0, 48.0D0, 64.0D0, 96.0D0, 128.0D0 /                         
      IF (IND) THEN 
!.                                                                      
!     determine machine constant                                        
!.                                                                      
         FMACHP = 1.0D0 
    5    FMACHP = 0.5D0 * FMACHP 
         IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 5 
         FMACHP = FMACHP * 2.0D0 
!.                                                                      
!     determine the maximal level of extrapolation. Due to              
!     increasing rounding errors effects, extrapolation should not be   
!     performed at a higher level.                                      
!.                                                                      
         IEXTMX = INT ( - DLOG (FMACHP) / DLOG (2.0D0) / 7.0D0 + 0.5D0) 
         IND = .FALSE. 
      ENDIF 
!.                                                                      
!     Checking input data, initializing variables.                      
!.                                                                      
      IF (ABSERR.LT.0.0D0.OR.RELERR.LT.0.0D0.OR.ABSERR +                &
      RELERR.LE.0.0D0) THEN                                             
         IERR = 3 
         RETURN 
      ELSEIF (XEND.EQ.X) THEN 
         IERR = 4 
         RETURN 
      ELSEIF (HMAX.LE.0.0D0) THEN 
         IERR = 5 
         RETURN 
      ENDIF 
      ICOUNT = 0 
      IF (IERR.EQ.1) GOTO 150 
      ILINE = 0 
      YMAX = 0.0D0 
!.                                                                      
!     determin the first step size                                      
!.                                                                      
   20 ABSH = DABS (H) 
      IF (HMAX.LT.ABSH) ABSH = HMAX 
      HLOC = DSIGN (DMIN1 (ABSH, DABS (XEND-X) ), XEND-X) 
      IF (DABS (HLOC) .LE.FMACHP * DABS (X) ) THEN 
         IERR = 0 
         RETURN 
      ENDIF 
!.                                                                      
!     determine the step size for the extrapolation.                    
!.                                                                      
   30 ILINE = ILINE+1 
      H0 = HLOC / BUFOL (ILINE) 
!.                                                                      
!     EULER step, store the initial values.                             
!.                                                                      
      X0 = X 
      DO 40 I = 1, N 
         AUXF (I, 9) = Y (I) 
   40 END DO 
      CALL FCT (X0, AUXF (1, 9), AUXF (1, 12), N) 
      DO 50 I = 1, N 
         AUXF (I, 10) = AUXF (I, 9) + H0 * AUXF (I, 12) 
   50 END DO 
      X0 = X0 + H0 
      CALL FCT (X0, AUXF (1, 10), AUXF (1, 12), N) 
!.                                                                      
!     the midpoint rule is applied                                      
!     BUFOL(ILINE)-1 - times.                                           
!.                                                                      
      IBUFO = INT (BUFOL (ILINE) ) - 1 
      DO 90 NLOOPS = 1, IBUFO 
         DO 60 I = 1, N 
            AUXF (I, 11) = AUXF (I, 9) + 2.0D0 * H0 * AUXF (I, 12) 
   60    END DO 
         X0 = X0 + H0 
         CALL FCT (X0, AUXF (1, 11), AUXF (1, 12), N) 
!.                                                                      
!     alter storage for the next step.                                  
!.                                                                      
         DO 80 J = 9, 11 
            DO 70 I = 1, N 
               AUXF (I, J) = AUXF (I, J + 1) 
   70       END DO 
   80    END DO 
   90 END DO 
!.                                                                      
!     stabilize using the trapezoidal rule.                             
!.                                                                      
      CALL FCT (X0, AUXF (1, 10), AUXF (1, 12), N) 
      DO 100 I = 1, N 
         AUXF (I, ILINE) = 0.5D0 * (AUXF (I, 10) + AUXF (I, 9) + H0 *   &
         AUXF (I, 12) )                                                 
  100 END DO 
!.                                                                      
!     at least two values are required for the extrapolation.           
!.                                                                      
      IF (ILINE.EQ.1) GOTO 30 
!.                                                                      
!     extrapolation.                                                    
!.                                                                      
      MINZ = MIN0 (ILINE, IEXTMX) 
      DO 120 ICOLUM = 2, MINZ 
         IDUMMY = MIN0 (11 - ICOLUM, ILINE-ICOLUM + 2) 
         DMAX = 0.0D0 
         DO 110 I = 1, N 
            AUXF (I, IDUMMY - 1) = AUXF (I, IDUMMY) + (AUXF (I, IDUMMY) &
            - AUXF (I, IDUMMY - 1) ) / ( (BUFOL (ILINE) / BUFOL (IDUMMY &
            - 1) ) **2 - 1)                                             
            YMAX = DMAX1 (YMAX, DABS (AUXF (I, IDUMMY - 1) ) ) 
            DMAX = DMAX1 (DMAX, DABS (AUXF (I, IDUMMY - 1) - AUXF (I,   &
            IDUMMY) ) )                                                 
  110    END DO 
         IF (DMAX.LT.RELERR * YMAX + ABSERR) GOTO 180 
  120 END DO 
      ICOUNT = ICOUNT + INT (BUFOL (ILINE) ) + 1 
      IF (ICOUNT.GE.ICTMAX) THEN 
!.                                                                      
!        return with message for excessive functional evaluations,      
!        store the current data.                                        
!.                                                                      
         IERR = 1 
         AUXF (1, 10) = X 
         DO 130 I = 1, N 
            AUXF (I, 9) = Y (I) 
  130    END DO 
         X = X0 
         DO 140 I = 1, N 
            Y (I) = AUXF (I, IDUMMY - 1) 
  140    END DO 
         RETURN 
      ENDIF 
      GOTO 170 
!.                                                                      
!     entry point on repeat call following return due to                
!     excessive functional evaluations.                                 
!.                                                                      
  150 X = AUXF (1, 10) 
      DO 160 I = 1, N 
         Y (I) = AUXF (I, 9) 
  160 END DO 
  170 IF (ILINE.LT.8) GOTO 30 
!.                                                                      
!     despite a completed extrapolation scheme the                      
!     required accuracyprecision was not achieved,                      
!     calculations are repeated with a smaller step size.               
!.                                                                      
      H = 0.9D0 * 0.6D0** (ILINE-IEXTMX) * H 
      IF (DABS (H) .LE.4.0D0 * FMACHP * DABS (X0) ) THEN 
         IERR = 2 
         RETURN 
      ENDIF 
      ILINE = 0 
      GOTO 30 
!.                                                                      
!     the required accuracy was achieved.                               
!.                                                                      
  180 X = X + HLOC 
      DO 190 I = 1, N 
         Y (I) = AUXF (I, IDUMMY - 1) 
  190 END DO 
!.                                                                      
!     the next step size is determined.                                 
!.                                                                      
      H = 0.9D0 * 0.6D0** (ICOLUM - IEXTMX) * H 
      ILINE = 0 
      GOTO 20 
      END SUBROUTINE DESEXT                         
