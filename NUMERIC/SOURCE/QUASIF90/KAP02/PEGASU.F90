      SUBROUTINE PEGASU (FCT, A, B, ABSERR, RELERR, MAXIT, XSI, X1, X2, &
      NUMIT, IERR)                                                      
!                                                                       
!*****************************************************************      
!                                                                *      
!  We assume that the given continuous function FCT(X) satisfies:*      
!          FCT(A) * FCT(B) < 0.0  on the interval [A,B].         *      
!  Thus there is at least one zero of odd order in [A,B].        *      
!  The SUBROUTINE PEGASU determines one of these zeros by        *      
!  applying the Pegasus-method.                                  *      
!  The method always converges provided the following condition  *      
!  is met:                                                       *      
!               FCT(A) * FCT(B) < 0.0                            *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  FCT   : function for which a zero is to be determined.        *      
!          It is declared as                                     *      
!              DOUBLE PRECISION FUNCTION FCT(X)                  *      
!          and has to be defined as EXTERNAL within the          *      
!          calling program (or as INTRINSIC if a FORTRAN         *      
!          standard function is used).                           *      
!  A,B   : endpoints of the interval containing a zero.          *      
!  ABSERR: ) error parameters. ABSERR >= 0.0 and RELERR >= 0.0.  *      
!  RELERR: ) Their sum has to be > 0.0. The following mixed test *      
!            is used as the break-off criterion:                 *      
!                 ABS(X1-X2) <= ABS(X2)*RELERR+ABSERR.           *      
!            If RELERR=0.0, this tests the absolute error,       *      
!            if ABSERR=0.0, this tests for the relative error.   *      
!            The values entered for ABSERR and RELERR are accep- *      
!            ted unchanged by the program if they both exceed    *      
!            four times the machine constant.                    *      
!            Or, if one value is zero and the other excceds      *      
!            four times the machine constan, they are accepted   *      
!            as well. If this is not the case both or one of the *      
!            values are internally set to this value.            *      
!  MAXIT : maximum number of iteration steps.                    *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  ABSERR: ) error parameters actually used.                     *      
!  RELERR: )                                                     *      
!  XSI   : zero or approximate value for the zero                *      
!          of the function FCT.                                  *      
!  X1,X2 : the last two iterates, so that [X1,X2] contains a     *      
!          zero of FCT.                                          *      
!  NUMIT : number of iteration steps executed.                   *      
!  IERR  : =-2, ABSERR or RELERR is negative or both equal zero, *      
!               or MAXIT < 1.                                    *      
!          =-1, the condition FCT(A)*FCT(B) < 0.0 is not met.    *      
!          = 0, A or B already are a zero of FCT.                *      
!          = 1, XSI is a zero with ABS(FCT(XSI)) < 4* machine    *      
!               constant.                                        *      
!          = 2, break-off criterion has been met, XSI:=X2 if     *      
!                   ABS(FCT(X2)) <= ABS(FCT(X1)),                *      
!               otherwise XSI:=X1.                               *      
!               The absolute error of the computed zero is less  *      
!               than or equal to  ABS(X1-X2).                    *      
!          = 3, the maximum number of iteration steps was        *      
!               reached without meeting the break-off criterion. *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: MACHPD                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Gisela Engeln-Muellges                           *      
!  date       : 09.02.1985                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
!  initializing the parameters NUMIT, X1, X2.                           
!                                                                       
      NUMIT = 0 
      X1 = A 
      X2 = B 
!                                                                       
!  calculating the machine constant FMACHP.                             
!                                                                       
      FMACHP = 1.0D0 
   10 FMACHP = 0.5D0 * FMACHP 
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 10 
      FMACHP = 2.0D0 * FMACHP 
!                                                                       
!  calculating the functional values at the endpoints A and B.          
!                                                                       
      F1 = FCT (X1) 
      F2 = FCT (X2) 
!                                                                       
!  testing for alternate signs of FCT at A and B:  FCT(A)*FCT(B) < 0.0 .
!                                                                       
      IF (F1 * F2.GT.0.0D0) THEN 
         IERR = - 1 
         RETURN 
      ELSEIF (F1 * F2.EQ.0.0D0) THEN 
         IERR = 0 
         RETURN 
      ENDIF 
!                                                                       
!  testing for validity of the error bounds and MAXIT.                  
!                                                                       
      IF (ABSERR.GE.0.0D0.AND.RELERR.GE.0.0D0.AND.ABSERR +              &
      RELERR.GT.0.0D0.AND.MAXIT.GE.1) GOTO 20                           
      IERR = - 2 
      RETURN 
   20 DUMMY = 4.0D0 * FMACHP 
      IF (RELERR.EQ.0.0D0) THEN 
         IF (ABSERR.LT.DUMMY) ABSERR = DUMMY 
      ELSEIF (ABSERR.EQ.0.0D0) THEN 
         IF (RELERR.LT.DUMMY) RELERR = DUMMY 
      ELSE 
         IF (ABSERR.LT.DUMMY) ABSERR = DUMMY 
         IF (RELERR.LT.DUMMY) RELERR = DUMMY 
      ENDIF 
!                                                                       
!  executing the Pegasus-method.                                        
!                                                                       
      DO 30 I = 1, MAXIT, 1 
         NUMIT = I 
!                                                                       
!     testing whether the value of F2 is less than four times           
!     the machine constant. If this is the case, X2 is accepted         
!     as a zero of FCT with IERR=1.                                     
!                                                                       
         IF (DABS (F2) .LT.4.0D0 * FMACHP) THEN 
            XSI = X2 
            IERR = 1 
            RETURN 
!                                                                       
!     testing for the break-off criterion.                              
!                                                                       
         ELSEIF (DABS (X2 - X1) .LE.DABS (X2) * RELERR + ABSERR) THEN 
            XSI = X2 
            IF (DABS (F1) .LT.DABS (F2) ) XSI = X1 
            IERR = 2 
            RETURN 
         ELSE 
!                                                                       
!     calculating the secant slope.                                     
!                                                                       
            S12 = (F2 - F1) / (X2 - X1) 
!                                                                       
!     calculating the secant intercept X3 with the x-axis.              
!                                                                       
            X3 = X2 - F2 / S12 
!                                                                       
!     calculating a new functional value at X3.                         
!                                                                       
            F3 = FCT (X3) 
!                                                                       
!     definition of the endpoints of a smaller inclusion interval.      
!                                                                       
            IF (F2 * F3.LE.0.0D0) THEN 
               X1 = X2 
               F1 = F2 
            ELSE 
               F1 = F1 * F2 / (F2 + F3) 
            ENDIF 
            X2 = X3 
            F2 = F3 
         ENDIF 
   30 END DO 
      IERR = 3 
      RETURN 
      END SUBROUTINE PEGASU                         
