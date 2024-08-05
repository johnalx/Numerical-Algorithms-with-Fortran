![  {The King and the Anderson--Bj"orck--King Method}                   
![  {The King and the Anderson--Bj"orck--King Methods,                  
![   the Illinois Method}*)                                             
      SUBROUTINE ZERORF (FCT, ABSERR, RELERR, FAB, MAXIT, IMETH, IEXTRA,&
      DELX, X1, X2, X3, F1, F2, F3, NUMIT, IHELP1, IHELP2, INCL, IERR)  
!                                                                       
!*****************************************************************      
!                                                                *      
!  The SUBROUTINE ZERORF is the governing program for several    *      
!  subroutines that determine zeros of a continuous real-valued  *      
!  function FCT. The methods used are primarily intended for     *      
!  simple zeros and zeros of odd order.                          *      
!  This SUBROUTINE adapts its approach depending on whether      *      
!  the starting values enclose a zero of the function FCT or not.*      
!  If a zero is enclosed, an iteration sequence is constructed   *      
!  using one of several iteration methods: Pegasus,              *      
!  Anderson/Bjoerck, or King. These methods ensure that inclusion*      
!  is preserved. All methods used work without derivatives and   *      
!  each possesses a convergence order of P > 1.6.                *      
!  If no two starting estimates are known that enclose a         *      
!  functional zero, then two identical starting values X1, X2    *      
!  may be specified. The SUBROUTINE ZERORF will construct a      *      
!  second starting value X2 by adding DELX to the starting value *      
!  X1. DELX has to be chosen by the user.                        *      
!  If there is still no enclosure, linear extrapolation is tried,*      
!  or, following the first secant step, quadratic extrapolation  *      
!  using the tangent line of an interpolating parabola through   *      
!  three different interpolation points is used. If enclosure    *      
!  is finally achieved, the program continues with one of the    *      
!  three special iterative methods as chosen by the user.        *      
!  NOTE: Extensive tests have shown that the number of iteration *      
!        steps required per zero for the same degree of          *      
!        accuracy is the largest for the Pegasus-method          *      
!        (SUBROUTINE PEG) and by far the smallest for the        *      
!        combination of the Anderson/Bjoerck and King            *      
!        methods (SUBROUTINE ANDBJK).                            *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  FCT      : real-valued function for which a zero is to be     *      
!             determined. It is declared as                      *      
!                 DOUBLE PRECISION FUNCTION FCT(X)               *      
!             and has to be defined as EXTERNAL within the       *      
!             calling program (or as INTRINSIC if a FORTRAN      *      
!             standard function is used).                        *      
!  ABSERR   : ) error bounds both of which have to be >= 0.0.    *      
!  RELERR   : ) Their sum has to be > 0.0. The following mixed   *      
!               test is used as a break-off criterion:           *      
!                   ABS(X1-X2) <= ABS(X2)*RELERR+ABSERR.         *      
!               Thus if RELERR=0.0 is chosen, this tests for the *      
!               absolute error, if ABSERR=0.0, this tests for the*      
!               relative error.                                  *      
!               The values entered for ABSERR and RELERR are     *      
!               accepted unchanged by the program if they        *      
!               both exceed four times the machine constant, or, *      
!               if one is zero, then the other has to exceed     *      
!               four times the machine constant. If this is not  *      
!               the case, then both or one of the bounds is set  *      
!               internally to this value.                        *      
!  FAB      : break-off criterion constant for the functional    *      
!             value at the last approximation of the zero. The   *      
!             value for FAB can be chosen between zero and four  *      
!             times the machine constant. If it is chosen        *      
!             negative or larger than four times the machine     *      
!             constant, it is set to that value internally.      *      
!  MAXIT    : maximum number of functional evaluations.          *      
!             (this equals the number of iteration steps)        *      
!  IMETH    : = 1, use the Pegasus-method.                       *      
!             = 2, use the King-method.                          *      
!             = 3, use the Anderson/Bjoerck-method.              *      
!             = 4, use the method of Anderson/Bjoerck and King.  *      
!  IEXTRA   : = 0, quadratic extrapolation shall be allowed.     *      
!             = 1, only linear extrapolation is acceptable (if   *      
!                  e.g. a double root is expected or if setting  *      
!                  IEXTRA=0 has resulted in IERR=0).             *      
!             If the starting values do not enclose a zero, i.e.,*      
!             if ( F1*F2 > 0.0, one may predetermine via IEXTRA, *      
!             whether only linear extrapolation should be tried  *      
!             or whether quadratic extrapolation is also         *      
!             acceptable. If no double zero is expected, then    *      
!             IEXTRA=0 should be used first.                     *      
!  DELX     : is used for altering a starting value in case two  *      
!             identical starting values were entered. A meaning- *      
!             ful choice for DELX would be:                      *      
!                     DELX=1E-6 or DELX=1E-8.                    *      
!  X1,X2    : starting values for the iteration; if only one     *      
!             starting value is known, setting X2=X1 is          *      
!             acceptable. By adding DELX to one of the two       *      
!             identical starting values, the program will create *      
!             two different starting values and proceed.         *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  ABSERR   : ) error bounds actually used.                      *      
!  RELERR   : )                                                  *      
!  X1,X2,X3 : approximate values for the desired zero.           *      
!             (see IERR).                                        *      
!  F1,F2,F3 : functional values FCT(X1), FCT(X2), FCT(X3).       *      
!  NUMIT    : number of functional evaluations performed.        *      
!  INCL     : = 0, starting values with no enclosure.            *      
!             = 1, starting values with known enclosure.         *      
!  IERR     : = 0, zero was not found.                           *      
!             = 1, zero lies between X1 and X2  (of the two      *      
!                  values the one with the smaller absolute      *      
!                  functional value should be chosen as the      *      
!                  zero).                                        *      
!                  The absolute error of the computed zero is    *      
!                  smaller than or equal to ABS(X1-X2).          *      
!             = 2, X2 is a zero of FCT: F2=0.0 (machine zero).   *      
!             = 3, X3 is a zero with ABS(F3) < 4 * machine       *      
!                  constant.                                     *      
!             = 4, zero of FCT is at X2 (enclosure interval      *      
!                  not definable).                               *      
!             = 5, maximum number MAXIT of functional evaluations*      
!                  exceeded.                                     *      
!             = 6, ABSERR or RELERR are negative, or both are    *      
!                  equal to zero, or MAXIT < 1.                  *      
!                                                                *      
!                                                                *      
!  AUXILIARY PARAMETERS:                                         *      
!  =====================                                         *      
!  IHELP1,IHELP2 : internal variables.                           *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: EXCHG, PEG, PEGK, ANDBJ, ANDBJK,        *      
!                        SUB1, SUB2, MACHPD                      *      
!                                                                *      
!                                                                *      
!  sources: 1. Anderson/Bjoerck, see [ANDE73].                   *      
!           2. Dowell/Jarrat, see [DOWE71], [DOWE72].            *      
!           3. King, see [KING73].                               *      
!           4. unpublished manuscript by R. Wodicka,             *      
!              RWTH Aachen.                                      *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Gisela Engeln-Muellges                           *      
!  date       : 08.26.1985                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
!  if two identical starting values are initially given (X1=X2),        
!  a second starting value is internally generated.                     
!                                                                       
      IF (X1.EQ.X2) THEN 
         X1 = X1 + DELX 
      ENDIF 
!                                                                       
!  initializing the parameters IERR and NUMIT.                          
!                                                                       
      IERR = 1 
      NUMIT = 2 
!                                                                       
!  calculation of the machine constant FMACHP.                          
!                                                                       
      FMACHP = 1.0D0 
   10 FMACHP = 0.5D0 * FMACHP 
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 10 
      FMACHP = 2.0D0 * FMACHP 
!                                                                       
!  testing the validity of the error bounds and MAXIT.                  
!                                                                       
      IF (ABSERR.GE.0.0D0.AND.RELERR.GE.0.0D0.AND.ABSERR +              &
      RELERR.GT.0.0D0.AND.MAXIT.GE.1) GOTO 20                           
      IERR = 6 
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
!  testing the validity of the break-off parameter FAB for the          
!  functional value at the approximate zero.                            
!                                                                       
      IF (FAB.LT.0.0D0.OR.FAB.GT.DUMMY) FAB = DUMMY 
!                                                                       
!  calculating the functional values at the starting points.            
!                                                                       
      F1 = FCT (X1) 
      F2 = FCT (X2) 
!                                                                       
!  labelling the zeros, so that ABS(F2) <= ABS(F1) holds.               
!                                                                       
      IF (DABS (F2) .GT.DABS (F1) ) THEN 
         CALL EXCHG (X1, X2, F1, F2) 
      ENDIF 
!                                                                       
!  test whether X2 already is a zero of FCT,                            
!  in which case IERR=2.                                                
!                                                                       
      IF (F2.EQ.0.0D0) THEN 
         IERR = 2 
         RETURN 
      ENDIF 
!                                                                       
!  test whether the starting values X1, X2 enclose a zero of FCT.       
!  If so, the internal variables are set to INCLUD=1, IHELP2=1,         
!  otherwise INCLUD=0.                                                  
!  IHELP2=1 indicates that the next iteration step is to be a secant ste
!                                                                       
      IF (F1 * F2.GT.0.0D0) THEN 
         INCLUD = 0 
         F3 = F2 
      ELSE 
         INCLUD = 1 
         IHELP2 = 1 
      ENDIF 
      INCL = INCLUD 
!                                                                       
!  iteration loop to find a new approximate value X3. We account for    
!  the fact here whether we have inclusion or not. The new              
!  approximate value X3 is taken as the x-intercept of the straight     
!  line connecting (X2, F2) and (X1, F1). Unless a secant step is       
!  performed, F1 is changed: in case of enclosure according to          
!  the specified method, otherwise by using quadratic extrapolation     
!  (see above remarks).                                                 
!                                                                       
   30 IF (INCLUD.EQ.0) THEN 
         F1DF2 = F1 / F2 
         IF (F1DF2.GT.1.0D0) THEN 
!                                                                       
!     check whether linear extrapolation is the only acceptable         
!     method or whether quadratic extrapolation is also allowed.        
!                                                                       
            IF (IEXTRA.EQ.0) THEN 
               IF ( (F1DF2 - F1 / F3) .GT.1.0D0) THEN 
                  G = 1.0D0 - F2 / F3 
                  F1 = G * F1 
               ENDIF 
            ENDIF 
         ELSE 
!                                                                       
!     now F1/F2 <= 1.0. If moreover F1/F2 < 0.0, then we have enclosure.
!                                                                       
            IF (F1DF2.LT.0.0D0) THEN 
               INCLUD = 1 
               IF (DABS (X1 - X2) .LE.DABS (X2) * RELERR + ABSERR) THEN 
                  IERR = 1 
                  RETURN 
               ELSE 
                  IHELP2 = 1 
               ENDIF 
            ELSE 
!                                                                       
!     if there is no enclosure, then a zero cannot be found             
!     because of ABS(F1) <= ABS(F2).                                    
!                                                                       
               IERR = 0 
               RETURN 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
!  calculation of the scaling factor Q for X3=X2+Q(X1-X2).              
!                                                                       
      Q = F2 / (F2 - F1) 
!                                                                       
!  calculation of the new approximate value.                            
!                                                                       
      X3 = X2 + Q * (X1 - X2) 
!                                                                       
!  testing whether the new approximate value X3 differs from both       
!  X1 and X2. If this is not the case, an alternate value               
!  X3NEW is calculated for X3. If there is no enclosure and             
!  X2=X3=X3NEW, then the program is stopped with setting                
!  IERR=4 (zero of FCT is at X2).                                       
!                                                                       
      IF (INCLUD.EQ.0) THEN 
         IF (X2.EQ.X3) THEN 
            X3NEW = X2 + (X2 - X1) / 9.0D0 
            IF (X2.EQ.X3NEW) THEN 
               IERR = 4 
               RETURN 
            ELSE 
               X3 = X3NEW 
            ENDIF 
         ELSE 
            IF (X2.EQ.X3) THEN 
               X3NEW = X2 + (X1 - X2) / 3.0D0 
               IF (X2.EQ.X3NEW) THEN 
                  IERR = 1 
                  RETURN 
               ELSE 
   40             Q = 2.0D0 * Q 
                  X3NEW = X2 + Q * (X1 - X2) 
                  IF (X3NEW.EQ.X2) GOTO 40 
                  IF (X3NEW.EQ.X1) THEN 
                     IERR = 1 
                     RETURN 
                  ELSE 
                     X3 = X3NEW 
                  ENDIF 
               ENDIF 
            ELSE 
               IF (X1.EQ.X3) THEN 
                  X3NEW = X1 + (X2 - X1) / 3.0D0 
                  IF (X3NEW.EQ.X1) THEN 
                     IERR = 1 
                     RETURN 
                  ELSE 
                     Q = F1 / (F1 - F2) 
   50                Q = 2.0D0 * Q 
                     X3NEW = X1 + Q * (X2 - X1) 
                     IF (X3NEW.EQ.X1) GOTO 50 
                     IF (X3NEW.EQ.X2) THEN 
                        IERR = 1 
                        RETURN 
                     ELSE 
                        X3 = X3NEW 
                     ENDIF 
                  ENDIF 
               ENDIF 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
!  now X3 differs from both X1 and X2, and F3 has been calculated.      
!                                                                       
      F3 = FCT (X3) 
!                                                                       
!  increasing the counter NUMIT.                                        
!                                                                       
      NUMIT = NUMIT + 1 
!                                                                       
!  test whether |F3| is less than four times the                        
!  machine constant. In this case, X3 is a zero with IERR=3.            
!                                                                       
      IF (DABS (F3) .LE.FAB) THEN 
         IERR = 3 
         RETURN 
      ENDIF 
!                                                                       
!  test whether the maximum number of functional evaluations            
!  allowed has been exceeded.                                           
!                                                                       
      IF (NUMIT.GE.MAXIT) THEN 
         IERR = 5 
         RETURN 
      ENDIF 
!                                                                       
!  the approximate values and their functional values                   
!  are relabelled by applying the SUBROUTINE EXCHG.                     
!  If F2*F3 > 0.0, the auxiliary variable IHELP1 is set to 0;           
!  if F2*F3 < 0.0, it is set to IHELP1=1.                               
!                                                                       
      IF (INCLUD.EQ.0) THEN 
         CALL EXCHG (X1, X2, F1, F2) 
      ELSE 
         IF (F2 * F3.LT.0.0D0) THEN 
            CALL EXCHG (X1, X2, F1, F2) 
            IHELP1 = 1 
         ELSE 
            IHELP1 = 0 
         ENDIF 
      ENDIF 
      CALL EXCHG (X2, X3, F2, F3) 
!                                                                       
!  if we have enclosure, we check whether the break-off criterion       
!  holds for the new enclosure interval. If the break-off condition     
!  is met, the iteration is stopped with IERR=1.                        
!  Otherwise the user specified method is used to determine F1.         
!                                                                       
      IF (INCLUD.EQ.1) THEN 
         IF (DABS (X1 - X2) .GT.DABS (X2) * RELERR + ABSERR) THEN 
            IF (IMETH.EQ.1) THEN 
               CALL PEG (IHELP1, F1, F2, F3) 
            ELSEIF (IMETH.EQ.2) THEN 
               CALL PEGK (IHELP1, IHELP2, F1, F2, F3) 
            ELSEIF (IMETH.EQ.3) THEN 
               CALL ANDBJ (IHELP1, F1, F2, F3) 
            ELSE 
               CALL ANDBJK (IHELP1, IHELP2, F1, F2, F3) 
            ENDIF 
         ELSE 
            IERR = 1 
            RETURN 
         ENDIF 
      ENDIF 
      GOTO 30 
      END SUBROUTINE ZERORF                         
!                                                                       
!                                                                       
      SUBROUTINE PEG (IHELP1, F1, F2, F3) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This SUBROUTINE uses the Pegasus-method to calculate a new    *      
!  functional value F1 for the governing program ZERORF.         *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  IHELP1   : auxiliary variable, as specified by the            *      
!             governing program ZERORF.                          *      
!  F1,F2,F3 : functional values FCT(X1), FCT(X2), FCT(X3).       *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  F1       : new functional value at X1.                        *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: SUB1                                    *      
!                                                                *      
!                                                                *      
!  sources: Dowell/Jarrat, see at [DOWE71], [DOWE72].            *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Gisela Engeln-Muellges                           *      
!  date       : 08.26.1985                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      IF (IHELP1.EQ.0) THEN 
         CALL SUB1 (F1, F2, F3) 
      ENDIF 
      RETURN 
      END SUBROUTINE PEG                            
!                                                                       
!                                                                       
!                                                                       
      SUBROUTINE PEGK (IHELP1, IHELP2, F1, F2, F3) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This SUBROUTINE uses an improved version of the Pegasus       *      
!  method by King for calculating a new functional value F1      *      
!  for use in the governing program ZERORF.                      *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  IHELP1   : ) auxiliary variables that are specified by the    *      
!  IHELP2   : ) governing program ZERORF.                        *      
!  F1,F2,F3 : functional values FCT(X1), FCT(X2), FCT(X3).       *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETER:                                             *      
!  =================                                             *      
!  F1       : new functional value at X1.                        *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: SUB1                                    *      
!                                                                *      
!                                                                *      
!  sources: method of King, see at [KING73].                     *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Gisela Engeln-Muellges                           *      
!  date       : 08.26.1985                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      IF (IHELP2.EQ.1) THEN 
         IHELP2 = 0 
         CALL SUB1 (F1, F2, F3) 
      ELSEIF (IHELP1.EQ.0) THEN 
         CALL SUB1 (F1, F2, F3) 
      ELSE 
         IHELP2 = 1 
      ENDIF 
      RETURN 
      END SUBROUTINE PEGK                           
!                                                                       
!                                                                       
!                                                                       
      SUBROUTINE ANDBJ (IHELP1, F1, F2, F3) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This SUBROUTINE uses the Anderson / Bjoerck method to calcu-  *      
!  late a new functional value F1 for the governing program      *      
!  ZERORF.                                                       *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  IHELP1   : auxiliary variable that is determined in the       *      
!             governing program ZERORF.                          *      
!  F1,F2,F3 : functional values FCT(X1), FCT(X2), FCT(X3).       *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETER:                                             *      
!  =================                                             *      
!  F1       : new functional value at X1.                        *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: SUB2                                    *      
!                                                                *      
!                                                                *      
!  sources: method of Anderson/Bjoerck, see at [ANDE73].        *       
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Gisela Engeln-Muellges                           *      
!  date       : 08.26.1985                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      IF (IHELP1.EQ.0) THEN 
         CALL SUB2 (F1, F2, F3) 
      ENDIF 
      RETURN 
      END SUBROUTINE ANDBJ                          
!                                                                       
!                                                                       
      SUBROUTINE ANDBJK (IHELP1, IHELP2, F1, F2, F3) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This SUBROUTINE uses a combination of the Anderson/Bjoerck    *      
!  and King methods suggested by King to calculate a new func-   *      
!  tional value F1 for use in the governing program ZERORF.      *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  IHELP1   : ) auxiliary variables set by the governing         *      
!  IHELP2   : ) program ZERORF.                                  *      
!  F1,F2,F3 : functional values FCT(X1), FCT(X2), FCT(X3).       *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETER:                                             *      
!  =================                                             *      
!  F1       : new functional value at X1.                        *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: SUB2                                    *      
!                                                                *      
!                                                                *      
!  sources: 1. method of Anderson/Bjoerck, see at [ANDE73].      *      
!           2. method of King, see at [KING73].                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Gisela Engeln-Muellges                           *      
!  date       : 08.26.1985                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      IF (IHELP2.EQ.1) THEN 
         IHELP2 = 0 
         CALL SUB2 (F1, F2, F3) 
      ELSEIF (IHELP1.EQ.0) THEN 
         CALL SUB2 (F1, F2, F3) 
      ELSE 
         IHELP2 = 1 
      ENDIF 
      RETURN 
      END SUBROUTINE ANDBJK                         
!                                                                       
!                                                                       
      SUBROUTINE EXCHG (X, Y, FX, FY) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This SUBROUTINE exchanges X and Y, and FX and FY.             *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DUMMY = X 
      X = Y 
      Y = DUMMY 
      DUMMY = FX 
      FX = FY 
      FY = DUMMY 
      RETURN 
      END SUBROUTINE EXCHG                          
!                                                                       
!                                                                       
      SUBROUTINE SUB1 (F1, F2, F3) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Auxiliary routine for subroutines PEG and PEGK.               *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      G = F3 / (F2 + F3) 
      F1 = G * F1 
      RETURN 
      END SUBROUTINE SUB1                           
!                                                                       
!                                                                       
      SUBROUTINE SUB2 (F1, F2, F3) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Auxiliary routine for subroutines ANDBJ and ANDBJK.           *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      G = 1.0D0 - F2 / F3 
      IF (G.LE.0.0D0) THEN 
         G = 0.5D0 
      ENDIF 
      F1 = G * F1 
      RETURN 
      END SUBROUTINE SUB2                           
