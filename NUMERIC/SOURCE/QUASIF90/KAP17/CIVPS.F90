      SUBROUTINE CIVPS (X, H, BETA, ABSERR, RELERR, N, FSAL, M, DES, Y, &
      EPS, XZI, QG, COEFF, MAXSTP, IERR)                                
!                                                                       
!*****************************************************************      
!                                                                *      
! This program solves a system of at most 12 ordinary            *      
! differential equations of first order by using a RUNGE-KUTTA   *      
! embedding formula over the interval of integration             *      
! I= [X0,BETA]. The step size is automatically controlled.       *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! X       : DOUBLE PRECISION initial value for the integration:  *      
!           X=X0                                                 *      
! H       : DOUBLE PRECISION initial step size                   *      
! BETA    : DOUBLE PRECISION endpoint X=BETA at which we want to *      
!           find the solution                                    *      
! ABSERR  : DOUBLE PRECISION error bound for the absolute error  *      
!           (ABSERR >= 0). If ABSERR=0, then only the relative   *      
!           error is checked.                                    *      
! RELERR  : DOUBLE PRECISION error bound for the relative error  *      
!           (RELERR >= 0). If RELERR=0, then only the absolute   *      
!           error is checked.                                    *      
! N       : number of differential equations in the system,      *      
!           or the size of Y:   0 < N < 13                       *      
! FSAL    : (LOGICAL) variable indicating whether the method     *      
!           FSAL (First Same As Last) is used by the RUNGE-KUTTA *      
!           embedding formula                                    *      
! M       : level of the embedding formula, also used for        *      
!           dimensioning COEFF.                                  *      
! DES     : SUBROUTINE DES must be declared as EXTERNAL in the   *      
!           calling program. DES describes the system of         *      
!           differential equations and must have the following   *      
!           form:                                                *      
!                  SUBROUTINE DES(X,Y,N,YPUNKT)                  *      
!                  DOUBLE PRECISION Y(N),YPUNKT(N),X             *      
!                  YPUNKT(1)=....                                *      
!                  YPUNKT(2)=....                                *      
!                         .                                      *      
!                         .                                      *      
!                         .                                      *      
!                   YPUNKT(N)=....                               *      
!                   RETURN                                       *      
!                   END                                          *      
! Y       : DOUBLE PRECISION vector Y(1:N), the solution at X=X0 *      
! EPS     : DOUBLE PRECISION 100 * machine constant              *      
! XZI     : DOUBLE PRECISION largest representable number for    *      
!           testing for OVERFLOW                                 *      
! QG      : global error order of the low order RUNGE-KUTTA      *      
!           method in use                                        *      
! COEFF   : 2-dim. DOUBLE PRECISION array COEFF(1:16,1:M) which  *      
!           contains the coefficients of the formula             *      
! MAXSTP  : maximal number of allowed integration steps          *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! X       : DOUBLE PRECISION value for X, where the integration  *      
!           has stopped (normally X=BETA)                        *      
! H       : DOUBLE PRECISION last step size used                 *      
! Y       : DOUBLE PRECISION solution vector Y(1:N) for X        *      
! IERR    : error parameter:                                     *      
!               IERR=0     all is ok                             *      
!               IERR=-1    the desired relative accuracy is less *      
!                          than 100 times the machine constant   *      
!                          in certain parts of the integration   *      
!                          interval. In these regions we compute *      
!                          with 100 times the machine constant as*      
!                          an absolute error bound.              *      
!               IERR=-2    the nunber of maximally allowed steps *      
!                          has been reached.                     *      
!               IERR=-20   OVERFLOW, the program stops.          *      
!               IERR=-30   the computed step size is too small.  *      
!                          The program stops.                    *      
!                                                                *      
!                                                                *      
! LOCAL VARIABLES:                                               *      
! ================                                               *      
! LSTSTP  : (LOGICAL)  LSTSTP=FALSE  continue integration        *      
!                      LSTSTP=TRUE   stop integration            *      
! ISTEP   : loop variable                                        *      
! J       : loop variable                                        *      
! YHILO   : 2-dim. DOUBLE PRECISION array YHILO(1:12,1:2)        *      
!           see SUBROUTINE RKSTEP or STEP32                      *      
! K       : 2-dim. DOUBLE PRECISION array K(1:12,1:16)           *      
!           see SUBROUTINE RKSTEP or STEP32                      *      
! YDIFF   : DOUBLE PRECISION auxiliary vector YDIFF(1:12)        *      
! NOSTEP  : LOGICAL variable, see SUBROUTINE RKSTEP or STEP32    *      
! FSALHP  : LOGICAL variable, auxiliary variable for FSAL        *      
! DELTA   : DOUBLE PRECISION estimate of the local error         *      
! EPSLON  : DOUBLE PRECISION tolerance for the maximal local     *      
!           error                                                *      
! S       : DOUBLE PRECISION auxiliary variable                  *      
! TEMP    : DOUBLE PRECISION auxiliary variable, the last step   *      
!           size                                                 *      
! EXPO    : DOUBLE PRECISION variable                            *      
! XDUMMY  : DOUBLE PRECISION variable                            *      
! XEND    : DOUBLE PRECISION variable for checking the endpoint  *      
!           of the interval                                      *      
! JERR    : error parameter of the SUBROUTINE RKSTEP or STEP32   *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: RKSTEP, STEP32, VMNORM                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : Volker Krger                                      *      
!  Date     : 28.04.1993                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
! Declarations                                                          
!                                                                       
      EXTERNAL DES 
      DOUBLEPRECISION Y (N), YHILO (12, 2), K (12, 16), YDIFF (12),     &
      COEFF (16, M)                                                     
      DOUBLEPRECISION X, H, TEMP, EXPO, BETA, ABSERR, RELERR, VMNORM,   &
      QG, DELTA, EPS, XZI, EPSLON, S, XDUMMY, XEND                      
      LOGICAL NOSTEP, FSAL, FSALHP, LSTSTP 
!                                                                       
! Initialize FSALHP, NOSTEP and LSTSTP                                  
!                                                                       
      FSALHP = .FALSE. 
      NOSTEP = .FALSE. 
      LSTSTP = .FALSE. 
!                                                                       
! Compute EXPO                                                          
!                                                                       
      EXPO = 1.0D0 / QG 
!                                                                       
! Integrate:                                                            
!   the maximal number of steps is MAXSTP                               
!                                                                       
      DO 20 ISTEP = 1, MAXSTP 
         IF (M.GT.3) THEN 
!                                                                       
!   SUBROUTINE RKSTEP performs one integration                          
!                                                                       
            CALL RKSTEP (X, H, Y, N, M, K, DES, YHILO, COEFF, NOSTEP,   &
            FSALHP, XZI, JERR)                                          
         ELSE 
!                                                                       
!   SUBROUTINE STEP32 performs one integration                          
!                                                                       
            CALL STEP32 (X, H, Y, N, K, DES, YHILO, NOSTEP, XZI, JERR) 
         ENDIF 
!                                                                       
!   after the first step we use the feature FSAL                        
!                                                                       
         FSALHP = FSAL 
!                                                                       
!                                                                       
!   If OVERFLOW is encounterd, return                                   
!   to calling program                                                  
!                                                                       
         DO 200 J = 1, N 
            IF (JERR.EQ.1.OR.DABS (YHILO (J, 1) ) .GT.XZI.OR.DABS (     &
            YHILO (J, 2) ) .GT.XZI) THEN                                
               IERR = - 20 
               RETURN 
            ENDIF 
  200    END DO 
!                                                                       
!   Determine step size for next step                                   
!                                                                       
         DO 30 J = 1, N 
            YDIFF (J) = YHILO (J, 1) - YHILO (J, 2) 
   30    END DO 
         DELTA = VMNORM (YDIFF, N) 
         EPSLON = ABSERR + RELERR * VMNORM (YHILO (1, 1), N) 
         IF (EPSLON.LT.EPS) THEN 
            EPSLON = EPS 
            IERR = - 1 
         ENDIF 
         EPSLON = DABS (H) * EPSLON 
         IF (DELTA.LT.EPS.OR. (DELTA.LE.1.0D+00.AND.XZI *               &
         DELTA.LT.EPSLON) ) THEN                                        
!                                                                       
!     Check prevent division by zero or OVERFLOW                        
!                                                                       
            S = 2.0D0 
         ELSE 
            S = (EPSLON / DELTA) **EXPO 
         ENDIF 
!                                                                       
!   S is less than 1                                                    
!                                                                       
         IF (S.LT.1.0D0) THEN 
            NOSTEP = .TRUE. 
!                                                                       
!   the new step size is found at least half the size                   
!   of the old one. If the step size becomes too small,                 
!   return to calling program                                           
!                                                                       
            H = H * DMAX1 (0.5D0, S) 
            IF (DABS (H) .LT.EPS) THEN 
               IERR = - 30 
               RETURN 
            ENDIF 
!                                                                       
!   If we repeat the last integration we adjust LSTSTP accordingly      
!                                                                       
            IF (LSTSTP) LSTSTP = .FALSE. 
!                                                                       
!   S is larger than 1                                                  
!                                                                       
         ELSE 
!                                                                       
!   Update X and Y                                                      
!                                                                       
            X = X + H 
            DO 10 J = 1, N 
               Y (J) = YHILO (J, 1) 
   10       END DO 
!                                                                       
!   the last integration was successful, the last                       
!   step size TEMP is used for H,                                       
!   return to calling program                                           
!                                                                       
            IF (LSTSTP) THEN 
               H = TEMP 
               RETURN 
!                                                                       
!   Prepare for new integration step                                    
!                                                                       
            ELSE 
               NOSTEP = .FALSE. 
!                                                                       
!   the new step size is chosen at most twice as large                  
!   as the last one                                                     
!                                                                       
               H = H * DMIN1 (2.0D0, S) 
!                                                                       
!   Check whether we have reached the endpoint of                       
!   the interval                                                        
!                                                                       
               XDUMMY = X + H 
               XEND = BETA - 0.1D0 * H 
               IF ( (H.LT.0.0D+00) .EQV. (XDUMMY.LT.XEND) ) THEN 
!                                                                       
!   Initialize LSTSTP once more.                                        
!   Store the last step size H in TEMP, and                             
!   recompute H                                                         
!                                                                       
                  LSTSTP = .TRUE. 
                  TEMP = H 
                  H = BETA - X 
               ENDIF 
            ENDIF 
         ENDIF 
   20 END DO 
!                                                                       
! the integration is stopped when the maximal allowable                 
! number of integration steps has been reached                          
!                                                                       
      IERR = - 2 
      RETURN 
      END SUBROUTINE CIVPS                          
