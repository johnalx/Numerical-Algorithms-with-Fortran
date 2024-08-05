      SUBROUTINE HULL (X, H, BETA, ABSERR, RELERR, N, FSAL, M, DES, Y,  &
      EPS, XZI, QG, COEFF, MAXSTP, IERR)                                
!                                                                       
!*****************************************************************      
!                                                                *      
! This program solves a system of at most 12 ordinary            *      
! differential equations of first order by using a RUNGE-KUTTA   *      
! embedding formula over the interval of integration             *      
! I= [X0,BETA]. The step size is automatically controlled using  *      
! the method of HULL.                                            *      
!                                                                *      
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
      QG, DELTA, EPS, XZI, EPSLON, XDUMMY, XEND                         
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
      EXPO = 1.0D0 / (QG + 1.0D0) 
!                                                                       
! Integrate:                                                            
!   the maximal number of steps is MAXSTP                               
!                                                                       
      DO 20 ISTEP = 1, MAXSTP 
         IF (M.GT.3) THEN 
!                                                                       
!   Integrate using the SUBROUTINE RKSTEP                               
!                                                                       
            CALL RKSTEP (X, H, Y, N, M, K, DES, YHILO, COEFF, NOSTEP,   &
            FSALHP, XZI, JERR)                                          
         ELSE 
!                                                                       
!   Integrate using the SUBROUTINE STEP32                               
!                                                                       
            CALL STEP32 (X, H, Y, N, K, DES, YHILO, NOSTEP, XZI, JERR) 
         ENDIF 
!                                                                       
!   after the first step we use the feature FSAL                        
!                                                                       
         FSALHP = FSAL 
!                                                                       
!   If we encounter OVERFLOW:                                           
!   return to calling program                                           
!                                                                       
         DO 200 J = 1, N 
            IF (JERR.EQ.1.OR.DABS (YHILO (J, 1) ) .GT.XZI.OR.DABS (     &
            YHILO (J, 2) ) .GT.XZI) THEN                                
               IERR = - 20 
               RETURN 
            ENDIF 
  200    END DO 
!                                                                       
!   Compute step size for next step                                     
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
!                                                                       
!   DELTA exxceeds the tolerance EPSLON                                 
!                                                                       
         IF (DELTA.GT.EPSLON) THEN 
            NOSTEP = .TRUE. 
!                                                                       
!   the new step size is limited to one fourth of the old step size.    
!   If it then falls below the machine constant, we return to the       
!   calling program                                                     
!                                                                       
            H = H * DMAX1 (2.5D-01, 9.0D-01 * (EPSLON / DELTA) **EXPO) 
            IF (DABS (H) .LT.EPS) THEN 
               IERR = - 30 
               RETURN 
            ENDIF 
!                                                                       
!   Adjust LSTSTP if the last integration is repeated                   
!                                                                       
            IF (LSTSTP) LSTSTP = .FALSE. 
!                                                                       
!   DELTA becomes smaller than the tolerance EPSLON                     
!                                                                       
         ELSE 
!                                                                       
!   Reinitialize X and Y                                                
!                                                                       
            X = X + H 
            DO 10 J = 1, N 
               Y (J) = YHILO (J, 1) 
   10       END DO 
!                                                                       
!   the last integration was successful, store the last step            
!   size TEMP in H, go back to calling program                          
!                                                                       
            IF (LSTSTP) THEN 
               H = TEMP 
               RETURN 
!                                                                       
!   prepare for a new integration                                       
!                                                                       
            ELSE 
               NOSTEP = .FALSE. 
!                                                                       
!   the new step size is bounded by 4 times the old one                 
!                                                                       
               IF (DELTA.LT.EPS.OR. (DELTA.LE.1.0D+00.AND.XZI *         &
               DELTA.LT.EPSLON) ) THEN                                  
!                                                                       
!     Check prevent division by zero or OVERFLOW                        
!                                                                       
                  H = 4.0D+00 * H 
               ELSE 
                  H = H * DMIN1 (4.0D+00, 9.0D-01 * (EPSLON / DELTA) ** &
                  EXPO)                                                 
!                                                                       
!   If the step size becomes too small, return to calling program       
!                                                                       
                  IF (DABS (H) .LT.EPS) THEN 
                     IERR = - 30 
                     RETURN 
                  ENDIF 
               ENDIF 
!                                                                       
!   Check whether the end of the interval has been reached              
!                                                                       
               XDUMMY = X + H 
               XEND = BETA - 0.1D0 * H 
               IF ( (H.LT.0.0D+00) .EQV. (XDUMMY.LT.XEND) ) THEN 
!                                                                       
!   Adjust LSTSTP accordingly, store the last step size TEMP in H.      
!   Compute H for the last integration step                             
!                                                                       
                  LSTSTP = .TRUE. 
                  TEMP = H 
                  H = BETA - X 
               ENDIF 
            ENDIF 
         ENDIF 
   20 END DO 
!                                                                       
! Integration is stopped after the maximally allowed                    
! number of integration steps                                           
!                                                                       
      IERR = - 2 
      RETURN 
      END SUBROUTINE HULL                           
