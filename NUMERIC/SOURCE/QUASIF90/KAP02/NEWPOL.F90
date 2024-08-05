      SUBROUTINE NEWPOL (A, N, X0, MAXIT, ABSERR, RELERR, X, PN, NUMIT, &
      IERR, WORK)                                                       
!                                                                       
!*****************************************************************      
!                                                                *      
!  This SUBROUTINE determines an approximation X for a zero of   *      
!  the algebraic polynomial                                      *      
!       PN(X)=A(0)+A(1)*X+A(2)*X**2+...+A(N)*X**N                *      
!  with real coefficients A(I) using Newton's method for simple  *      
!  roots. The iteration function is:                             *      
!              PHI(X):=X-PN(X)/PNDER(X).                         *      
!  Starting with X0, the description of the parameters is        *      
!  identical to that of the SUBROUTINE NEWPSZ, as is the organi- *      
!  zation of this program with the exception that the functional *      
!  evaluations and those for the derivative are determined       *      
!  by a Horner scheme (SUBROUTINE NEWPOH).                       *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  A        : (n+1) - vector A(0:N) containing the               *      
!             coefficients A(I) of PN.                           *      
!  N        : degree of the polynomial PN.                       *      
!  X0       : starting value for the iteration.                  *      
!  MAXIT    : maximum number of iteration steps (MAXIT >= 1).    *      
!  ABSERR   : ) error parameters. Both have to be >= 0.0, while  *      
!  RELERR   : ) their sum must be > 0.0. The following mixed     *      
!               break-off criterion is used:                     *      
!                    ABS(X1-X2) <= ABS(X2)*RELERR+ABSERR.        *      
!               Thus if RELERR=0.0 is chosen, this is a test for *      
!               the absolute error; if ABSERR=0.0 is chosen, this*      
!               is a test for the relative error.                *      
!               The values entered for ABSERR and RELERR are     *      
!               accepted unchanged by the program if they both   *      
!               exceed four times the machine constant, or if one*      
!               is zero, if the other exceeds four times the     *      
!               machine constant.                                *      
!               If this is not the case, both error bounds, or   *      
!               one of them are set to that value internally.    *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  ABSERR   : ) error parameters actually used.                  *      
!  RELERR   : )                                                  *      
!  X        : last approximate value for the zero X              *      
!             of PN(X).                                          *      
!  PN       : functional value PN(X) at the last approximate X.  *      
!  NUMIT    : number of iteration steps executed.                *      
!  IERR     : = 0, faulty input parameters.                      *      
!             = 1, zero X found; break-off criterion for the     *      
!                  difference between the last two approximations*      
!                  has been met.                                 *      
!             = 2, zero X found.                                 *      
!                  ABS(PN(X)) < 4 * machine constant.            *      
!             = 3, starting value X0 already is a zero of PN:    *      
!                         PN(X0)=0.0 (machine zero).             *      
!             = 4, zero has not been found. The maximum number   *      
!                  of iterations has been reached without        *      
!                  meeting one of the break-off criteria.        *      
!                  Possibly the error parameters have been chosen*      
!                  too small or the starting value X0 is not     *      
!                  close enough to a zero of PN.                 *      
!                                                                *      
!  AUXILIARY VECTOR:                                             *      
!  =================                                             *      
!  WORK     : (n+1) - vector WORK(0:N).                          *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: NEWPOH, MACHPD                          *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Gisela Engeln-Muellges                           *      
!  date       : 09.23.1985                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
!  declarations.                                                        
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
      DOUBLEPRECISION A (0:N), WORK (0:N) 
!                                                                       
      NUMIT = 0 
!                                                                       
!  test of input parameters.                                            
!                                                                       
      IF (ABSERR.LT.0.0D0.OR.RELERR.LT.0.0D0.OR.ABSERR +                &
      RELERR.LE.0.0D0.OR.MAXIT.LT.1) THEN                               
         IERR = 0 
         RETURN 
      ENDIF 
      X = X0 
      CALL NEWPOH (A, N, X, PN, PNDER, WORK) 
      IF (PN.EQ.0.0D0) THEN 
         IERR = 3 
         RETURN 
      ENDIF 
      FMACHP = 1.0D0 
   10 FMACHP = 0.5D0 * FMACHP 
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 10 
      FMACHP = 2.0D0 * FMACHP 
      DUMMY = 4.0D0 * FMACHP 
      IF (RELERR.EQ.0.0D0) THEN 
         IF (ABSERR.LT.DUMMY) ABSERR = DUMMY 
      ELSEIF (ABSERR.EQ.0.0D0) THEN 
         IF (RELERR.LT.DUMMY) RELERR = DUMMY 
      ELSE 
         IF (ABSERR.LT.DUMMY) ABSERR = DUMMY 
         IF (RELERR.LT.DUMMY) RELERR = DUMMY 
      ENDIF 
!                                                                       
!  iteration loop.                                                      
!                                                                       
   20 IF (DABS (PN) .LT.DUMMY) THEN 
         IERR = 2 
         RETURN 
      ELSE 
         IF (DABS (PNDER) .LT.DUMMY) THEN 
            PNDER = DSIGN (1.0D-6, PNDER) 
         ENDIF 
         NUMIT = NUMIT + 1 
         DIFF = PN / PNDER 
         X = X - DIFF 
         CALL NEWPOH (A, N, X, PN, PNDER, WORK) 
         IF (DABS (DIFF) .LE.DABS (X) * RELERR + ABSERR) THEN 
            IERR = 1 
            RETURN 
         ELSEIF (NUMIT.GE.MAXIT) THEN 
            IERR = 4 
            RETURN 
         ENDIF 
      ENDIF 
      GOTO 20 
      END SUBROUTINE NEWPOL                         
!                                                                       
!                                                                       
      SUBROUTINE NEWPOH (A, N, X, PN, PNDER, WORK) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This SUBROUTINE uses the Horner scheme in order to calculate  *      
!  the functional value PN(X) and the 1st derivative PNDER(X) of *      
!  an Nth degree algebraic polynomial with real coefficients     *      
!        PN(X)=A(0)+A(1)*X+A(2)*X**2+...+A(N)*X**N               *      
!  at X.                                                         *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  A        : (n+1) - vector A(0:N) containing the               *      
!             coefficients A(I) of PN.                           *      
!  N        : degree of the polynomial PN.                       *      
!  X        : value X for which the functional value or the 1st  *      
!             derivative value of the polynomial will be         *      
!             calculated.                                        *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  PN       : functional value of the polynomial at X.           *      
!  PNDER    : 1st derivative of the polynomial at X.             *      
!                                                                *      
!                                                                *      
!  AUXILIARY VECTOR:                                             *      
!  =================                                             *      
!  WORK     : (n+1) - vector WORK(0:N).                          *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Gisela Engeln-Muellges                           *      
!  date       : 09.23.1985                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION A (0:N), WORK (0:N) 
      DO 10 I = 0, N, 1 
         WORK (I) = A (I) 
   10 END DO 
      DO 20 K = 0, 1 
         S = 0.0D0 
         DO 20 I = N, K, - 1 
            S = S * X + WORK (I) 
            WORK (I) = S 
   20 CONTINUE 
      PN = WORK (0) 
      PNDER = WORK (1) 
      RETURN 
      END SUBROUTINE NEWPOH                         
