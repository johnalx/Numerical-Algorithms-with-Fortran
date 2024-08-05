      SUBROUTINE IVP (XK, HK, YK, N, DES, XE, EPSABS, EPSREL, INDEX,    &
      NMAX, NUSED, IERR)                                                
!                                                                       
!*****************************************************************      
!                                                                *      
!  For an approximation YK of the solution Y at XK of a system   *      
!  of ordinary differential equations of 1st order               *      
!                    Y' = F(X,Y),                                *      
!  this program computes an approximation for the solution Y     *      
!  at XE.                                                        *      
!  We use step size control in such a way that the error of the  *      
!  computed approximation falls either absolutely or relatively  *      
!  within the given error bounds EPSABS or EPSREL.               *      
!  In case of the Prince-Dormand embedding formula we also check *      
!  for stiffness of the system.                                  *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
! XK     - initial value of the independent variable X           *      
! HK     - proposed step size for the next step                  *      
! YK     - vector YK(1:N); value of the solution of the          *      
!          differential equation at XK                           *      
! N      - number of differential equations ( 1 <= N <= 20 )     *      
! DES    - right-hand side of the differential equation which    *      
!          must be provided as a SUBROUTINE in the form:         *      
!             SUBROUTINE DES (X, Y, F, N)                        *      
!          (starting with:  DOUBLE PRECISION Y(N), F(N), etc.).  *      
!          Here F represents the value of the differential       *      
!          equation's right-hand side at (X,Y). (DES has to be   *      
!          defined as EXTERNAL in the calling program)           *      
! XE     - location where the solution is desired; XE may not    *      
!          be chosen smaller than XK.                            *      
! EPSABS - error bound for the absolute accuracy of the desired  *      
!          solution. EPSABS has to be >= 0; if EPSABS = 0, only  *      
!          the relative accuracy is considered.                  *      
! EPSREL - error bound for the relative accuracy of the desired  *      
!          solution. EPSREL has to be >= 0; if EPSREL = 0,       *      
!          only the absolute accuracy is considered.             *      
! INDEX  - chooses the embedding formula with step size control: *      
!             = 0: RUNGE-KUTTA method  2nd/3rd order             *      
!             =-1: Prince-Dormand embedding formula of 4th/5th   *      
!                  order (checking for stiffness of the system   *      
!                  of the DEs here, see output IERR)             *      
!            else: formula of ENGLAND  4th/5th order             *      
! NMAX   - upper limit for the number of function evaluations    *      
!          allowed for the right-hand side F.                    *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
! XK     - location that was reached during last integration.    *      
!          If IERR = 0 normally XK is equal to XE                *      
! HK     - local step size used last ( it should remain un-      *      
!          changed for the next step )                           *      
! YK     - approximate value of the solution Y at the new        *      
!          location XK                                           *      
! NUSED  - number of function evaluations actually used          *      
! IERR   - error parameter:                                      *      
!          = 0: everything o.k.                                  *      
!               If Prince-Dormand was specified, the sytsem was  *      
!               not found to be stiff.                           *      
!          = 1: both error bounds  EPS...  are too small         *      
!               (relative to the machine constant)               *      
!          = 2: XE <= XK       (wrt machine constant)            *      
!          = 3: step size  HK <= 0   (wrt machine constant)      *      
!          = 4: N > 20  or  N <= 0                               *      
!          = 5: NUSED > NMAX:  the number of allowed functional  *      
!               evaluations is insufficient to determine an      *      
!               adequate approximate solution with the required  *      
!               accuracy; on termination, XK and HK contain      *      
!               the current values                               *      
!          =-1: The computations terminated ok, but the Prince-  *      
!               Dormand formula has detected possible stiffness. *      
!          =-2: The computations have terminated ok, but the     *      
!               Prince-Dormand formula has recognized the system *      
!               as stiff using two criteria: we recommend to use *      
!               a method suited for stiff DEs instead.           *      
!          =-3: NUSED > NMAX:  The number of allowed functional  *      
!               evaluations does not suffice, see IERR = 5.      *      
!               Moreover, the Prince-Dormand formula seems to    *      
!               indicate that the system is stiff; we recommend  *      
!               to use a suitable stiff DE solver.               *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: DVNORM, RUKU23, ENGL45, MACHPD, PRDO45  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author      : Klaus Niederdrenk                               *      
!  Date        : 11.01.1985 / 4.1.1995                           *      
!  Source code : FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION YK (N) 
      DOUBLEPRECISION Y (20), YT (20), Y00 (20) 
      LOGICAL IEND, LEPSI, STEIF1, STEIF2, STEIFA 
      EXTERNAL DES 
      SAVE FMACHP, EPS1, EPS2, LEPSI 
      DATA LEPSI / .TRUE. /, Y00 / 20 * 0.0D0 / 
!                                                                       
!** FMACHP is the machine constant                                      
!** (i.e., the smallest positive machine number for which  1 + FMACHP > 
!** EPS1 bounds the admissable initial step size HK below.              
!** EPS2 is used to determine a level for negligable quantities.        
!                                                                       
      IF (LEPSI) THEN 
         FMACHP = 1.0D0 
   10    FMACHP = 0.5D0 * FMACHP 
         IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 10 
         FMACHP = 2.0D0 * FMACHP 
         EPS1 = FMACHP**0.75D0 
         EPS2 = 100.0D0 * FMACHP 
         LEPSI = .FALSE. 
      ENDIF 
!                                                                       
!** Preassign local variables                                           
!                                                                       
      SG = DSIGN (1.0D0, XE) 
      XEND = (1.0D0 - SG * EPS2) * XE 
      IERR = 0 
      NUSED = 0 
      IEND = .FALSE. 
      STEIF1 = .FALSE. 
      STEIF2 = .FALSE. 
      STEIFA = .FALSE. 
      IANZ = 0 
!                                                                       
!** Check input parameters                                              
!                                                                       
      YMAX = DVNORM (YK, Y00, N) 
      IF (EPSABS.LE.EPS2 * YMAX.AND.EPSREL.LE.EPS2) THEN 
         IERR = 1 
      ELSEIF (XEND.LT.XK) THEN 
         IERR = 2 
      ELSEIF (HK.LT.EPS2 * DABS (XK) ) THEN 
         IERR = 3 
      ELSEIF (N.LE.0.OR.N.GT.20) THEN 
         IERR = 4 
      ENDIF 
      IF (IERR.NE.0) RETURN 
!                                                                       
!**********  C O N T R O L L I N G   a l g o r i t h  **********        
!                                                                       
      IF (XK + HK.GT.XEND) THEN 
         HK = XE-XK 
         DUMMY = HK 
         IEND = .TRUE. 
      ENDIF 
!                                                                       
!** Integrate on the interval *[*XK, XE*]* in suitable steps            
!                                                                       
   50 CONTINUE 
!                                                                       
!** Call the desired one step method                                    
!                                                                       
      IF (INDEX.EQ.0) THEN 
         CALL RUKU23 (XK, HK, YK, N, DES, Y, YT) 
         NUSED = NUSED+3 
      ELSEIF (INDEX.EQ. - 1) THEN 
         CALL PRDO45 (XK, HK, YK, N, DES, Y, YT, STEIF1, STEIFA) 
         NUSED = NUSED+7 
         IF (STEIFA) THEN 
            IANZ = IANZ + 1 
            IF (IANZ.GE.3) STEIF2 = .TRUE. 
         ELSE 
            IANZ = 0 
         ENDIF 
      ELSE 
         CALL ENGL45 (XK, HK, YK, N, DES, Y, YT) 
         NUSED = NUSED+6 
      ENDIF 
!                                                                       
      DIFF = DVNORM (Y, YT, N) 
!                                                                       
      IF (DIFF.LT.EPS2) THEN 
         S = 2.0D0 
      ELSE 
         YMAX = DVNORM (YT, Y00, N) 
         S = DSQRT (HK * (EPSABS + EPSREL * YMAX) / DIFF) 
         IF (INDEX.NE.0) S = DSQRT (S) 
      ENDIF 
!                                                                       
      IF (S.GT.1.0D0) THEN 
!                                                                       
!** accept the performed integration with step size HK                  
!                                                                       
         DO 60 I = 1, N 
            YK (I) = YT (I) 
   60    END DO 
!                                                                       
         XK = XK + HK 
!                                                                       
!** if then endpoint XE has been reached or if more than the            
!** allowable function evaluations were used: go back                   
!                                                                       
   70    IF (IEND) THEN 
            HK = DUMMY 
            IF (INDEX.EQ. - 1) THEN 
               IF (STEIF1.OR.STEIF2) IERR = - 1 
               IF (STEIF1.AND.STEIF2) IERR = - 2 
            ENDIF 
            RETURN 
         ELSEIF (NUSED.GT.NMAX) THEN 
            IERR = 5 
            IF (INDEX.EQ. - 1.AND. (STEIF1.OR.STEIF2) ) IERR = - 3 
            RETURN 
         ENDIF 
!                                                                       
!** increase the step size for the next step maximally by a factor of tw
!                                                                       
         HK = HK * DMIN1 (2.0D0, 0.98D0 * S) 
!                                                                       
         IF ( (XK + HK) .GE.XEND) THEN 
            DUMMY = HK 
            HK = XE-XK 
            IEND = .TRUE. 
!                                                                       
!** if very close to XE: go back                                        
!                                                                       
            IF (HK.LT.EPS1 * DABS (XE) ) GOTO 70 
         ENDIF 
      ELSE 
!                                                                       
!** the previous step is unaccaptable, the step size HK must be decrease
!** at most it must be halved                                           
!                                                                       
         HK = HK * DMAX1 (0.5D0, 0.98D0 * S) 
         IEND = .FALSE. 
      ENDIF 
!                                                                       
      GOTO 50 
!                                                                       
      END SUBROUTINE IVP                            
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION DVNORM (F1, F2, N) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This program finds the maximum norm of the difference F1 - F2 *      
!  of two vectors F1 and F2 of length N.                         *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  required subroutines: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author      : Klaus Niederdrenk                               *      
!  Date        : 11.01.1985                                      *      
!  Source code : FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION F1 (N), F2 (N) 
      DVNORM = 0.0D0 
      DO 10 I = 1, N 
         DVNORM = DMAX1 (DVNORM, DABS (F1 (I) - F2 (I) ) ) 
   10 END DO 
      RETURN 
      END FUNCTION DVNORM                           
!                                                                       
!                                                                       
      SUBROUTINE RUKU23 (X, H, Y, N, DES, Y2, Y3) 
!                                                                       
!*****************************************************************      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Starting with an approximation Y at X, this program computes  *      
!  approximations Y2 and Y3 at X + H using the RUNGE-KUTTA       *      
!  embedding formula of 2nd and 3rd order for a system of        *      
!  ordinary differential equations of 1st order                  *      
!                  Y' = F(X,Y).                                  *      
!  The system of N ordinary differential equations of 1st order  *      
!  must be provided by a SUBROUTINE DES.                         *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  X     - initial value for the independent variable X          *      
!  H     - step size                                             *      
!  Y     - vector Y(1:N); value for the solution of the          *      
!          differential equation at X                            *      
!  N     - number of differential equations  ( 1 <= N <= 20 )    *      
!  DES   - right hand side of the differential equation. It has  *      
!          provided by the user as a SUBROUTINE in the form:     *      
!             SUBROUTINE DES (X, Y, N, F)                        *      
!          (starting with: DOUBLE PRECISION Y(N), F(N), etc. ).  *      
!          Here F denotes the right hand side of the differential*      
!          equation at (X,Y). ( In the calling program DES has to*      
!          declared as EXTERNAL)                                 *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  Y2    - vector Y2(1:N); approximate solution using the 2nd    *      
!          order method for solving the differential equation    *      
!          at X+H                                                *      
!  Y3    - vector Y3(1:N); approximate solution using the 3rd    *      
!          order method for solving the differential equation    *      
!          at X+H                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required : none                                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Klaus Niederdrenk                                  *      
!  date     : 11.01.1985                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION Y (N), Y2 (N), Y3 (N), DUMMY (20) 
      DOUBLEPRECISION K1 (20), K2 (20), K3 (20) 
!                                                                       
      CALL DES (X, Y, N, K1) 
      DO 10 I = 1, N 
         DUMMY (I) = Y (I) + H * K1 (I) 
   10 END DO 
      CALL DES (X + H, DUMMY, N, K2) 
      DO 20 I = 1, N 
         DUMMY (I) = Y (I) + 0.25D0 * H * (K1 (I) + K2 (I) ) 
   20 END DO 
      CALL DES (X + 0.5D0 * H, DUMMY, N, K3) 
!                                                                       
      DO 100 I = 1, N 
         Y2 (I) = Y (I) + 0.5D0 * H * (K1 (I) + K2 (I) ) 
         Y3 (I) = Y (I) + H / 6.0D0 * (K1 (I) + K2 (I) + 4.0D0 * K3 (I) &
         )                                                              
  100 END DO 
      RETURN 
      END SUBROUTINE RUKU23                         
!                                                                       
!                                                                       
      SUBROUTINE ENGL45 (X, H, Y, N, DES, Y4, Y5) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Starting from an approximation Y at X, this program uses the  *      
!  ENGLAND embedding formula of 4th and 5th order to find appro- *      
!  ximations Y4 and Y5 at X + H for the solution of the system of*      
!  differential equations of 1st order                           *      
!                  Y' = F(X,Y) .                                 *      
!  The system contains N ordinary differential equations of 1st  *      
!  order, with the right hand side F(X,Y) provided by the user   *      
!  in a SUBROUTINE DES.                                          *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  X     - initial value for the independent variable X          *      
!  H     - step size                                             *      
!  Y     - vector Y(1:N); value for the solution of the          *      
!          differential equation at X                            *      
!  N     - number of differential equations  ( 1 <= N <= 20 )    *      
!  DES   - right hand side of the differential equation. It has  *      
!          provided by the user as a SUBROUTINE in the form:     *      
!             SUBROUTINE DES (X, Y, N, F)                        *      
!          (starting with: DOUBLE PRECISION Y(N), F(N), etc. ).  *      
!          Here F denotes the right hand side of the differential*      
!          equation at (X,Y). ( In the calling program DES has to*      
!          declared as EXTERNAL)                                 *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  Y4    - vector Y4(1:N); approximate solution using the 4th    *      
!          order method for solving the differential equation    *      
!          at X+H                                                *      
!  Y5    - vector Y5(1:N); approximate solution using the 5th    *      
!          order method for solving the differential equation    *      
!          at X+H                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required : none                                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Klaus Niederdrenk                                  *      
!  date     : 11.01.1985                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION Y (N), Y4 (N), Y5 (N) 
      DOUBLEPRECISION DUMMY (20) 
      DOUBLEPRECISION K1 (20), K2 (20), K3 (20), K4 (20), K5 (20),      &
      K6 (20)                                                           
!                                                                       
      CALL DES (X, Y, N, K1) 
      DO 10 I = 1, N 
         DUMMY (I) = Y (I) + 0.5D0 * H * K1 (I) 
   10 END DO 
      CALL DES (X + 0.5D0 * H, DUMMY, N, K2) 
!                                                                       
      DO 20 I = 1, N 
         DUMMY (I) = Y (I) + 0.25D0 * H * (K1 (I) + K2 (I) ) 
   20 END DO 
      CALL DES (X + 0.5D0 * H, DUMMY, N, K3) 
!                                                                       
      DO 30 I = 1, N 
         DUMMY (I) = Y (I) + H * ( - K2 (I) + 2.0D0 * K3 (I) ) 
   30 END DO 
      CALL DES (X + H, DUMMY, N, K4) 
!                                                                       
      DO 40 I = 1, N 
         DUMMY (I) = Y (I) + H / 27.0D0 * (7.0D0 * K1 (I) + 10.0D0 * K2 &
         (I) + K4 (I) )                                                 
   40 END DO 
      CALL DES (X + 2.0D0 / 3.0D0 * H, DUMMY, N, K5) 
!                                                                       
      DO 50 I = 1, N 
         DUMMY (I) = Y (I) + 0.16D-02 * H * (28.0D0 * K1 (I) - 125.0D0 *&
         K2 (I) + 546.0D0 * K3 (I) + 54.0D0 * K4 (I) - 378.0D0 * K5 (I) &
         )                                                              
   50 END DO 
      CALL DES (X + 0.2D0 * H, DUMMY, N, K6) 
!                                                                       
      DO 100 I = 1, N 
         Y4 (I) = Y (I) + H / 6.0D0 * (K1 (I) + 4.0D0 * K3 (I) + K4 (I) &
         )                                                              
         Y5 (I) = Y (I) + H / 336.0D0 * (14.0D0 * K1 (I) + 35.0D0 * K4 (&
         I) + 162.0D0 * K5 (I) + 125.0D0 * K6 (I) )                     
  100 END DO 
      RETURN 
      END SUBROUTINE ENGL45                         
!                                                                       
!                                                                       
      SUBROUTINE PRDO45 (X, H, Y, N, DES, Y4, Y5, ST1, ST2) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Starting from an approximation Y at X, this program uses the  *      
!  Prince-Dormand embedding formula of 4th and 5th order to find *      
!  apprximations Y4 and Y5 at X + H for the solution of the      *      
!  system of differential equations of 1st order                 *      
!                  Y' = F(X,Y) .                                 *      
!  The system contains N ordinary differential equations of 1st  *      
!  order, with the right hand side F(X,Y) provided by the user   *      
!  in a SUBROUTINE DES.                                          *      
!  This program tests for stiffness in two ways.                 *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  X     - initial value for the independent variable X          *      
!  H     - step size                                             *      
!  Y     - vector Y(1:N); value for the solution of the          *      
!          differential equation at X                            *      
!  N     - number of differential equations  ( 1 <= N <= 20 )    *      
!  DES   - right hand side of the differential equation. It has  *      
!          provided by the user as a SUBROUTINE in the form:     *      
!             SUBROUTINE DES (X, Y, N, F)                        *      
!          (starting with: DOUBLE PRECISION Y(N), F(N), etc. ).  *      
!          Here F denotes the right hand side of the differential*      
!          equation at (X,Y). ( In the calling program DES has to*      
!          be declared as EXTERNAL)                              *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  Y4    - vector Y4(1:N); approximate solution using the 4th    *      
!          order method for solving the differential equation    *      
!          at X+H                                                *      
!  Y5    - vector Y5(1:N); approximate solution using the 5th    *      
!          order method for solving the differential equation    *      
!          at X+H                                                *      
!  ST1   - logical variable: .TRUE. , if the test for stiffness  *      
!          (using the approximate dominat eigenvalue) was        *      
!          positive ; else the value of ST! remains unaltered.   *      
!  ST2   - logical variable: .TRUE. , if the test of stiffness   *      
!          due to  Shampine  and  Hiebert  is positive; else the *      
!          value of ST2 is set to  .FALSE. .                     *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required : DVNORM                                 *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Klaus Niederdrenk                                  *      
!  date     : 1.4.1995                                           *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION Y (N), Y4 (N), Y5 (N) 
      DOUBLEPRECISION DUMMY (20) 
      DOUBLEPRECISION K1 (20), K2 (20), K3 (20), K4 (20), K5 (20),      &
      K6 (20), K7 (20), G6 (20), G7 (20)                                
      LOGICAL ST1, ST2 
!                                                                       
      CALL DES (X, Y, N, K1) 
      DO 10 I = 1, N 
         DUMMY (I) = Y (I) + 0.2D0 * H * K1 (I) 
   10 END DO 
      CALL DES (X + 0.2D0 * H, DUMMY, N, K2) 
!                                                                       
      DO 20 I = 1, N 
         DUMMY (I) = Y (I) + 0.075D0 * H * (K1 (I) + 3.0D0 * K2 (I) ) 
   20 END DO 
      CALL DES (X + 0.3D0 * H, DUMMY, N, K3) 
!                                                                       
      DO 30 I = 1, N 
         DUMMY (I) = Y (I) + H / 45.0D0 * (44.0D0 * K1 (I) - 168.0D0 *  &
         K2 (I) + 160.0D0 * K3 (I) )                                    
   30 END DO 
      CALL DES (X + 0.8D0 * H, DUMMY, N, K4) 
!                                                                       
      DO 40 I = 1, N 
         DUMMY (I) = Y (I) + H / 6561.0D0 * (19372.0D0 * K1 (I) -       &
         76080.0D0 * K2 (I) + 64448.0D0 * K3 (I) - 1908.0D0 * K4 (I) )  
   40 END DO 
      CALL DES (X + 8.0D0 / 9.0D0 * H, DUMMY, N, K5) 
!                                                                       
      DO 50 I = 1, N 
         G6 (I) = Y (I) + H / 167904.0D0 * (477901.0D0 * K1 (I) -       &
         1806240.0D0 * K2 (I) + 1495424.0D0 * K3 (I) + 46746.0D0 * K4 ( &
         I) - 45927.0D0 * K5 (I) )                                      
   50 END DO 
      CALL DES (X + H, G6, N, K6) 
!                                                                       
      DO 60 I = 1, N 
         G7 (I) = Y (I) + H / 142464.0D0 * (12985.0D0 * K1 (I) +        &
         64000.0D0 * K3 (I) + 92750.0D0 * K4 (I) - 45927.0D0 * K5 (I)   &
         + 18656.0D0 * K6 (I) )                                         
   60 END DO 
      CALL DES (X + H, G7, N, K7) 
!                                                                       
      DO 100 I = 1, N 
         Y5 (I) = G7 (I) 
         Y4 (I) = Y (I) + H / 21369600.0D0 * (1921409.0D0 * K1 (I)      &
         + 9690880.0D0 * K3 (I) + 13122270.0D0 * K4 (I) - 5802111.0D0 * &
         K5 (I) + 1902912.0D0 * K6 (I) + 534240.0D0 * K7 (I) )          
  100 END DO 
!                                                                       
!**  Test for stiffness via dominant eigenvalue (approximately)         
!                                                                       
      IF (DVNORM (K7, K6, N) .GT.3.3 * DVNORM (G7, G6, N) ) ST1 =       &
      .TRUE.                                                            
!                                                                       
!**  One step of testing stiffness according to Shampine and Hiebert    
!                                                                       
      DO 110 I = 1, N 
         G6 (I) = H * (2.2D0 * K2 (I) + 0.13D0 * K4 (I) + 0.144D0 * K5 (&
         I) )                                                           
         G7 (I) = H * (2.134D0 * K1 (I) + 0.24D0 * K3 (I) + 0.1D0 * K6 (&
         I) )                                                           
  110 END DO 
      IF (DVNORM (G6, G7, N) .LT.DVNORM (Y4, Y5, N) ) THEN 
         ST2 = .TRUE. 
      ELSE 
         ST2 = .FALSE. 
      ENDIF 
      RETURN 
      END SUBROUTINE PRDO45                         
