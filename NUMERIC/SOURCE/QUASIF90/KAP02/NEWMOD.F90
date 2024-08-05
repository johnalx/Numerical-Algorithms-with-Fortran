![  {Newton's Method for Multiple Zeros}                                
![  {Newton's Method for Multiple Zeros;                                
![   a Modified Newton's Method}*)                                      
      SUBROUTINE NEWMOD (FCT, FDER, F2DER, X0, MAXIT, LUN, ABSERR,      &
      RELERR, X, F, J, NUMIT, IERR)                                     
!                                                                       
!*****************************************************************      
!                                                                *      
!  This SUBROUTINE determines an approximation X for a J-fold    *      
!  zero of the function FCT by applying the modified Newton      *      
!  method for the starting value X0.                             *      
!  The iteration function is :                                   *      
!             PHI:=X-J(X)*FCT(X)/FDER(X), where                  *      
!          J(X)=1./(1.-FCT(X)*F2DER(X)/FDER(X)**2).              *      
!  The values J(X) converge towards the order J of the zero X.   *      
!  The method converges quadratically. However, its efficiency   *      
!  is only E=1.26, since each iteration step requires three      *      
!  functional evaluations.                                       *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  FCT      : function whose zero is to be determined.           *      
!             It has the form                                    *      
!                 DOUBLE PRECISION FUNCTION FCT(X)               *      
!             and must be declared as EXTERNAL by the            *      
!             calling program (or as INTRINSIC if a FORTRAN      *      
!             standard function is used).                        *      
!  FDER     : 1st derivative of the function FCT. It has the     *      
!             form                                               *      
!                 DOUBLE PRECISION FUNCTION FDER(X)              *      
!             and has to be declared as EXTERNAL by the          *      
!             calling program (or as INTRINSIC if a FORTRAN      *      
!             standard function is used).                        *      
!  F2DER    : 2nd derivative of the function FCT. It has the     *      
!             form                                               *      
!                 DOUBLE PRECISION FUNCTION F2DER(X)             *      
!             and has to be declared as EXTERNAL in the          *      
!             calling program (or as INTRINSIC if a FORTRAN      *      
!             standard function is used).                        *      
!  X0       : starting value for the iteration.                  *      
!  MAXIT    : maximum number of iteration steps (MAXIT >= 1).    *      
!  LUN      : > 0, tape unit onto which the iteration values X   *      
!             and their corresponding orders J(X) are stored.    *      
!  ABSERR   : ) error bounds that both must be >= 0.0, while     *      
!  RELERR   : ) their sum has to be > 0.0. The following mixed   *      
!               test will be used inside the subroutine:         *      
!                  ABS(X1-X2) <= ABS(X2)*RELERR+ABSERR.          *      
!               Thus if RELERR was chosen as 0.0, this tests for *      
!               the absolute error; if ABSERR=0.0 is chosen this *      
!               tests for the relative error.                    *      
!               The values chosen for ABSERR and RELERR are      *      
!               used unaltered by the program if they both       *      
!               exceed the machine constant by a factor of five, *      
!               unless one of their values is zero, in which     *      
!               case the other must exceed the machine constant  *      
!               by a factor of five.                             *      
!               If this is not the case, then both, or one of    *      
!               the error bounds are set to five times the       *      
!               machine constant internally.                     *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  ABSERR   : ) error parameters actually used.                  *      
!  RELERR   : )                                                  *      
!  X        : last approximate value for the zero of FCT.        *      
!  F        : function value of FCT at the last                  *      
!             approximate value X.                               *      
!  J        : order of the zero X.                               *      
!             (see remark below)                                 *      
!  NUMIT    : number of iteration steps performed.               *      
!  IERR     : = 0, input parameter are faulty                    *      
!             = 1, zero found; break-off criterion for the       *      
!                  difference of the last two approxiations      *      
!                  has been met.                                 *      
!             = 2, zero X found.                                 *      
!                  ABS(FCT(X)) < 5 * machine constant.           *      
!             = 3, starting value X0 already is a zero           *      
!                  of FCT(X0)=0.0 (machine zero).                *      
!             = 4, zero not found. The maximum number of         *      
!                  iterations reached without meeting the        *      
!                  break-off criterion.                          *      
!                  The error bounds might have been chosen too   *      
!                  small or the starting value X0 is not close   *      
!                  enough to a zero of FCT.                      *      
!                                                                *      
!                                                                *      
!                                                                *      
!  NOTE  : The quotients FCT/FDER and FCT*F2DER/FDER**2 become   *      
!          indeterminant in a neighborhood of the (multiple)     *      
!          zero. Thus, even for small error bounds the zero can  *      
!          only be determined with limited precision.            *      
!          In some rare test cases the value for the order J of  *      
!          the zero was computed incorrectly due to the inde-    *      
!          terminate expression for J(X). Therefore the results  *      
!          should always be inspected on TAPE unit # LUN and     *      
!          interpreted according to [ENGE87], p.51.              *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: MACHPD                                  *      
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
!                                                                       
!  initializing the iteration step counter NUMIT.                       
!                                                                       
      NUMIT = 0 
!                                                                       
!  test for validity of the input parameters ABSERR, RELERR and MAXIT.  
!                                                                       
      IF (ABSERR.LT.0.0D0.OR.RELERR.LT.0.0D0.OR.ABSERR +                &
      RELERR.LE.0.0D0.OR.MAXIT.LT.1) THEN                               
         IERR = 0 
         RETURN 
      ENDIF 
!                                                                       
!  test whether X0 already is already a zero of FCT: FCT(X0)=0.0        
!  (machine zero).                                                      
!                                                                       
      X = X0 
      F = FCT (X) 
      IF (F.EQ.0.0D0) THEN 
         IERR = 3 
         RETURN 
      ENDIF 
!                                                                       
!  compute the machine constant and, if necessary,                      
!  adjust the error bounds.                                             
!                                                                       
      FMACHP = 1.0D0 
   10 FMACHP = 0.5D0 * FMACHP 
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 10 
      FMACHP = 2.0D0 * FMACHP 
      DUMMY = 5.0D0 * FMACHP 
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
      WRITE (LUN, 900) 
   20 IF (DABS (F) .LT.DUMMY) THEN 
         IERR = 2 
         WRITE (LUN, 910) NUMIT, X 
         WRITE (LUN, * ) 
         RETURN 
      ELSE 
         FD = FDER (X) 
         F2D = F2DER (X) 
         IF (DABS (FD) .LT.DUMMY) THEN 
            FD = DSIGN (1.0D-6, FD) 
         ENDIF 
         XJ = 1.0D0 / (1.0D0 - F * F2D / (FD * FD) ) 
         WRITE (LUN, 910) NUMIT, X, XJ 
         DIFF = XJ * F / FD 
         X = X - DIFF 
         NUMIT = NUMIT + 1 
         F = FCT (X) 
         J = INT (XJ + 0.5D0) 
!                                                                       
!  test whether the break-off criterion is met.                         
!                                                                       
         IF (DABS (DIFF) .LT.DABS (X) * RELERR + ABSERR) THEN 
            IERR = 1 
            WRITE (LUN, 910) NUMIT, X 
            WRITE (LUN, * ) 
            RETURN 
         ELSEIF (NUMIT.GE.MAXIT) THEN 
            IERR = 4 
            WRITE (LUN, 910) NUMIT, X 
            WRITE (LUN, * ) 
            RETURN 
         ENDIF 
      ENDIF 
      GOTO 20 
  900 FORMAT(1X,'Results of the iteration:',/,1X,25('*'),//,3X,         &
     &       'I',1X,'I',9X,'X(I)',9X,'I',7X,'J(X(I))',/,1X,             &
     &       49('='),/,5X,'I',22X,'I')                                  
  910 FORMAT(1X,I3,2(' I ',D20.14)) 
      END SUBROUTINE NEWMOD                         
