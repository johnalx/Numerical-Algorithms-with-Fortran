![KA{P 2}{Nonlinear Equations in One Variable}                          
![       {Nonlinear Equations in One Variable}*)                        
![  {The Newton Method for Simple Roots}                                
![  {The Newton Method for Simple Roots}*)                              
      SUBROUTINE NEWPSZ (FCT, FDER, X0, MAXIT, ABSERR, RELERR, X, F,    &
      NUMIT, IERR)                                                      
!                                                                       
!*****************************************************************      
!                                                                *      
!  This SUBROUTINE determines an approximate value X for a zero  *      
!  of the function FCT with known derivative FDER by applying    *      
!  Newton's method for simple zeros. The iteration function is:  *      
!                 PHI(X):=X-FCT(X)/FDER(X).                      *      
!  For a simple zero the order of convergence is P=2 with        *      
!  efficiency E=1.414... If FCT is an algebraic polynomial the   *      
!  SUBROUTINE NEWPOL should be applied instead.                  *      
!  For multiple zeros the modified Newton method                 *      
!  (SUBROUTINE NEWMOD) is recommended. Although the order of     *      
!  convergence is slightly lower for methods that do not use     *      
!  derivatives such as SUBROUTINE ZERORF, PEGASU and others,     *      
!  they are preferable to Newton's method. Their efficiency is   *      
!  generally better since they only require one functional       *      
!  evaluation per iteration step.                                *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  FCT      : function for which a zero is to be found.          *      
!             It has the form                                    *      
!                 DOUBLE PRECISION FUNCTION FCT(X)               *      
!             and has to be declared as EXTERNAL within the      *      
!             calling program (or as INTRINSIC if a FORTRAN      *      
!             standard function is used).                        *      
!  FDER     : 1st derivative of the function FCT.                *      
!             It is declared as                                  *      
!                 DOUBLE PRECISION FUNCTION FDER(X)              *      
!             and has to be defined as EXTERNAL within the       *      
!             calling program (or as INTRINSIC, if a FORTRAN     *      
!             standard function is used).                        *      
!  X0       : starting value for the iteration.                  *      
!  MAXIT    : maximum number of iteration steps (MAXIT >= 1).    *      
!  ABSERR   : ) error parameters. Both have be set to >= 0.0     *      
!  RELERR   : ) their sum has to be > 0.0. The following mixed   *      
!               test is used:                                    *      
!                   ABS(X1-X2) <= ABS(X2)*RELERR+ABSERR.         *      
!               Thus if RELERR=0.0, we test for the absolute     *      
!               error; if ABSERR=0.0 is chosed, then the         *      
!               relative error is tested for.                    *      
!               The values entered for ABSERR and RELERR are     *      
!               accepted unchanged by the program if they both   *      
!               exceed four times the machine constant, or in    *      
!               case one of them equals zero, the other exceeds  *      
!               four times the machine constant.                 *      
!               If this is not the case for both or one of the   *      
!               values, they are set to that value internally.   *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  ABSERR   : ) error parameters actually used.                  *      
!  RELERR   : )                                                  *      
!  X        : last approximate value for the zero of FCT         *      
!  F        : function value of FCT at the last approximate      *      
!             value X.                                           *      
!  NUMIT    : number of iteration steps performed. The number    *      
!             of functional evaluations is two times NUMIT.      *      
!  IERR     : = 0, input parameter are not as specified.         *      
!             = 1, zero found; break-off criterion for the       *      
!                  difference between the last two approximations*      
!                  has been met.                                 *      
!             = 2, zero X found.                                 *      
!                      ABS(FCT(X)) < 4 * machine constant.       *      
!             = 3, starting value X0 is a zero of FCT:           *      
!                      FCT(X0)=0.0 (machine zero).               *      
!             = 4, zero not found. The maximum number of         *      
!                  iterations was reached without meeting        *      
!                  the break-off criterion. Possibly the         *      
!                  error parameters were chosen too small or     *      
!                  the starting value X0 was not close enough to *      
!                  a zero position of FCT.                       *      
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
!  testing the validity of the input parameters ABSERR, RELERR and MAXIT
!                                                                       
      IF (ABSERR.LT.0.0D0.OR.RELERR.LT.0.0D0.OR.ABSERR +                &
      RELERR.LE.0.0D0.OR.MAXIT.LT.1) THEN                               
         IERR = 0 
         RETURN 
      ENDIF 
!                                                                       
!  test whether X0 already is a zero of FCT:                            
!        FCT(X0)=0.0 (machine zero).                                    
!                                                                       
      X = X0 
      F = FCT (X) 
      IF (F.EQ.0.0D0) THEN 
         IERR = 3 
         RETURN 
      ENDIF 
!                                                                       
!  computation of the machine constant FMACHP.                          
!                                                                       
      FMACHP = 1.0D0 
   10 FMACHP = 0.5D0 * FMACHP 
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 10 
      FMACHP = 2.0D0 * FMACHP 
!                                                                       
!  in case ABSERR and/or RELERR were chosen too small, they are now set 
!  to four times the machine constant.                                  
!                                                                       
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
!  iteration loop for calculating a new approximate value X for         
!  the desired zero. First it we check whether F=FCT(X) is less         
!  than four times the machine constant.                                
!  If this is true the program is ended by setting IERR=2.              
!                                                                       
   20 IF (DABS (F) .LT.DUMMY) THEN 
         IERR = 2 
         RETURN 
      ELSE 
!                                                                       
!  calculation of the 1st derivative FD of the function FCT at X;       
!  if this turns out to be zero, the derivative will be set to 1.E-6    
!  in order to prevent division by zero.                                
!                                                                       
         FD = FDER (X) 
         IF (DABS (FD) .LT.DUMMY) THEN 
            FD = DSIGN (1.0D-6, FD) 
         ENDIF 
!                                                                       
!  increase the iteration step counter by 1 and calculate               
!  a new approximation for the zero and its corresponding               
!  functional value.                                                    
!                                                                       
         NUMIT = NUMIT + 1 
         DIFF = F / FD 
         X = X - DIFF 
         F = FCT (X) 
!                                                                       
!  test whether the break-off criterion is met for the last two         
!  approximations or whether the maximum number MAXIT of                
!  iterations has been reached.                                         
!                                                                       
         IF (DABS (DIFF) .LE.DABS (X) * RELERR + ABSERR) THEN 
            IERR = 1 
            RETURN 
         ELSEIF (NUMIT.GE.MAXIT) THEN 
            IERR = 4 
            RETURN 
         ENDIF 
      ENDIF 
      GOTO 20 
      END SUBROUTINE NEWPSZ                         
