      SUBROUTINE BVP (A, B, H, YSTART, N, DEQ, BDCD, IPROC, EPSIVP,     &
      EPSBC, IFMAX, ITMAX, ITER, IERR, AUXF, LDF, WORK1, WORK21, WORK22,&
      WORK23)                                                           
!                                                                       
!*****************************************************************      
!                                                                *      
!  This program is designed to solve a general first order       *      
!  boundary value problem given in system form as                *      
!                                                                *      
!    Y' = F(X,Y)  for  A <= X <= B  with   R(Y(A), Y(B)) = 0.    *      
!                                                                *      
!  It uses the shooting method to determine an approximation     *      
!  YSTART for the correct initial value Y(A), which then can be  *      
!  used to find the solution Y using an initial value problem    *      
!  solver via the SUBROUTINE IVP. Here IPROC serves as a label   *      
!  to select the desired initial value problem solver of IVP.    *      
!  The non-linear system of equations that occurs in the         *      
!  shooting method is solved by  NEWTON's method.                *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  A      - left endpoint of the interval of integration         *      
!  B      - right endpoint of the interval; B must be > A        *      
!  H      - appropriate starting step size for the approximate   *      
!           solution of the associated initial value problem for *      
!           the shooting method                                  *      
!  YSTART - starting approximation for the initial value Y(A) of *      
!           the solution Y of the boundary value problem         *      
!  N      - number of differential equations with N > 1 and:     *      
!           for IPROC = 1 :         N <= 20;                     *      
!                     = 2, 4 or 5 : N <= 100;                    *      
!                     = 3 :         N <= 12.                     *      
!  DEQ    - right-hand side of the differential equation, i.e.,  *      
!           the inhomogeneity, that has to be provided as a      *      
!           SUBROUTINE in the form                               *      
!              SUBROUTINE DEQ (X, Y, N, F)                       *      
!           (beginning with: DOUBLE PRECISION Y(N), F(N)).       *      
!           In it F denotes the inhomogeneity of the differential*      
!           equation at (X,Y). In the calling program DEQ has to *      
!           be declared as EXTERNAL.                             *      
!  BDCD   - boundary condition that has to be provided as a      *      
!           SUBROUTINE of the form                               *      
!              SUBROUTINE BDCD (YA, YB, N, R)                    *      
!           (starting with: DOUBLE PRECISION YA(N), YB(N), R(N)).*      
!           Here R denotes the value of R(YA,YB). In the calling *      
!           program BDCD has to be listed as EXTERNAL.           *      
!  IPROC  - label to select the initial value problem solver     *      
!           to be used inside the shooting algorithm:            *      
!           = 1 : Runge-Kutta embedding formula of 4th/5th order *      
!                 (England formula) from the subroutine IVP;     *      
!                 Restriction: 1 < N <= 20 .                     *      
!           = 2 : Predictor-corrector method of 4th order of     *      
!                 Adams-Bashforth-Moulton (subroutine DESABM;    *      
!                 constant IFMAX set to 1000; needs auxiliary    *      
!                 array AUXF ).                                  *      
!           = 3 : Runge-Kutta embedding formula of 7th/8th order *      
!                 (subroutine RKTRB; 1 < N <= 12 ; IFMAX = 10000)*      
!           = 4 : Extrapolation method of Bulirsch-Stoer;        *      
!                 (subroutine DESEXT; IFMAX = 1400; auxiliary    *      
!                 array  AUXF ).                                 *      
!           = 5 : implicit Runge-Kutta-Gauss method;             *      
!                 (subroutine IRKDRV; IFMAX = 200000, auxiliary  *      
!                 arrays  WORK1 , WORK21 , WORK22  and WORK23 ). *      
!  EPSIVP - accuracy bound for the approximate solution of the   *      
!           corresponding initial value problems for the shooting*      
!           method.                                              *      
!  EPSBC  - accuracy bound for the approximation YSTART for Y(A) *      
!           in the boundary condition.                           *      
!  IFMAX  - upper bound for the number of allowed functional     *      
!           evaluations of the right-hand side F while solving   *      
!           an associated initial value problem.                 *      
!           (concerns only IPROC = 1; otherwise set as follows:  *      
!           if IPROC = 2 : IFMAX =   1000;                       *      
!                    = 3 : IFMAX =  10000;                       *      
!                    = 4 : IFMAX =   1400;                       *      
!                    = 5 : IFMAX = 200000. )                     *      
!  ITMAX  - upper bound for the number of NEWTON-iterations while*      
!           solving the non-linear system of equations in the    *      
!           shooting method.                                     *      
!  AUXF   - array of size (1:LDF,1:12) where LDF >= N ; auxiliary*      
!           array for IPROC = 2 or 4                             *      
!  LDF    - leading dimension of AUXF; defined in calling program*      
!  WORK1  - vector of length (1:5*N+78) ; used for IPROC = 5.    *      
!  WORK21 - array of size (1:N,1:N) ; used for IPROC = 5.        *      
!           bei der Wahl IPROC = 5                               *      
!  WORK22 - array of size (1:N+9,1:9) ; for IPROC = 5.           *      
!  WORK23 - array of size (1:2*N+10,1:10) ; for IPROC = 5.       *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  YSTART - approximation for the initial value Y(A) of the      *      
!           solution Y of the boundary value problem             *      
!  ITER   - number of NEWTON-iterations actually performed       *      
!  IERR   - error parameter:                                     *      
!           = 0: everything o.k.                                 *      
!           = 1: at least one of the error parameters EPS... is  *      
!                too small    (relative to machine constant)     *      
!           = 2: B <= A       (within the machine constant)      *      
!           = 3: step width H <= 0                               *      
!                             (relative to machine constant)     *      
!           = 4: N > NMAX or N <= 0 where                        *      
!                    NMAX = 20  , for IPROC = 1                  *      
!                    NMAX = 100 , for IPROC = 2 , 4 or 5         *      
!                    NMAX = 12  , for IPROC = 3 .                *      
!           = 5: IPROC <= 0 or IPROC > 5                         *      
!           = 6: IFMAX function evaluations are not sufficient   *      
!                for approximately solving the associated initial*      
!                value problem with the shooting method.         *      
!           = 7: ITER > ITMAX: the number of allowed NEWTON-     *      
!                iterations is not sufficient to determine an    *      
!                approximate initial value YSTART within the     *      
!                desired accuracy.                               *      
!           = 8: the JACOBI matrix for the NEWTON method is      *      
!                numerically singular; it is impossible to       *      
!                perform a NEWTON-iteration step.                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  required directory (for IPROC = 5): directory number 9        *      
!                                                                *      
!  required subroutines: GAUSS, MACHPD and, depending on IPROC:  *      
!                        one of: IVP, DESABM, RKTRB, DESEXT, or  *      
!                                IRKDRV                          *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Klaus Niederdrenk                                  *      
!  date     : 11.01.1985 / 5.31.1994                             *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION YSTART (N) 
      DIMENSION AMAT (100, 100), YK (100), YK1 (100), YAJ (100),        &
      R (100)                                                           
      DIMENSION RJ (100), D (100), SCAL (100) 
      DIMENSION WOR1RK (4), WOR2RK (16, 16), IWORRK (2), IFLARK (3) 
      DIMENSION AUXF (LDF, 12) 
      DIMENSION WORK1 (5 * N + 78), WORK21 (N, N), WORK22 (N + 9, 9) 
      DIMENSION WORK23 (2 * N + 10, 10), G (100) 
      DIMENSION IPIVOT (100) 
      EXTERNAL DEQ 
!                                                                       
!** FMACHP is the machine constant of the computer being used           
!** EPS1 serves as a step size bound for forming the approximate Jacobi 
!** matrix for the Newton method.                                       
!** EPS2 is used to find quantities sufficiently close to zero.         
!                                                                       
      FMACHP = 1.0D0 
   10 FMACHP = 0.5D0 * FMACHP 
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 10 
      FMACHP = 2.0D0 * FMACHP 
      EPS1 = DSQRT (FMACHP) 
      EPS2 = 100.0D0 * FMACHP 
!                                                                       
!** Check input parameters                                              
!                                                                       
      IERR = 0 
      IF (EPSBC.LT.EPS2) THEN 
         IERR = 1 
      ELSEIF (B.LE.A) THEN 
         IERR = 2 
      ELSEIF (H.LT.EPS2 * DABS (A) ) THEN 
         IERR = 3 
      ELSEIF (IPROC.EQ.1) THEN 
         IF (N.LE.0.OR.N.GT.20) IERR = 4 
      ELSEIF (IPROC.EQ.2.OR.IPROC.EQ.4.OR.IPROC.EQ.5) THEN 
         IF (N.LE.0.OR.N.GT.100) IERR = 4 
      ELSEIF (IPROC.EQ.3) THEN 
         IF (N.LE.0.OR.N.GT.12) IERR = 4 
      ELSEIF (IPROC.LE.0.OR.IPROC.GT.5) THEN 
         IERR = 5 
      ENDIF 
      IF (IERR.NE.0) RETURN 
!                                                                       
!** Preassign local parameters                                          
!                                                                       
      DATA G / 100 * 1.0D0 / 
      ITER = 0 
      EPSABS = 0.5D0 * EPSIVP 
      EPSREL = EPSABS 
      IND = 1 
      IFLARK (1) = 11 
      IFLARK (2) = 0 
      IFLARK (3) = 0 
      IFLRKG = 1 
!                                                                       
!** in case  YSTART  is a sufficiently good approximation               
!** of Y(A) : Return                                                    
!                                                                       
   50 CONTINUE 
      DO 60 I = 1, N 
         YK (I) = YSTART (I) 
   60 END DO 
      XK = A 
      HK = H 
!                                                                       
      IF (IPROC.EQ.1) THEN 
         CALL IVP (XK, HK, YK, N, DEQ, B, EPSABS, EPSREL, 1, IFMAX, IFA,&
         IFEHL)                                                         
         IF (IFEHL.EQ.5) IERR = 6 
      ELSEIF (IPROC.EQ.2) THEN 
         IFEHL = 0 
         CALL DESABM (XK, YK, DEQ, N, B, HK, B - A, EPSABS, EPSREL, IND,&
         IFEHL, AUXF, LDF)                                              
         IF (IFEHL.EQ.2.OR.IFEHL.EQ.3) THEN 
            HK = B - XK 
            IFEHL = 0 
            CALL DESABM (XK, YK, DEQ, N, B, HK, B - XK, EPSABS, EPSREL, &
            IND, IFEHL, AUXF, LDF)                                      
            IF (IFEHL.EQ.2) IERR = 6 
         ENDIF 
      ELSEIF (IPROC.EQ.3) THEN 
         CALL RKTRB (XK, B, N, DEQ, YK, EPSABS, EPSREL, IFLARK, WOR1RK, &
         WOR2RK, IWORRK, IFEHL)                                         
         IFLARK (3) = 1 
         WOR1RK (4) = H 
         IF (IFEHL.EQ. - 2) IERR = 6 
      ELSEIF (IPROC.EQ.4) THEN 
         IFEHL = 0 
         CALL DESEXT (XK, YK, DEQ, N, B, HK, B - A, EPSABS, EPSREL,     &
         IFEHL, AUXF, LDF)                                              
         IF (IFEHL.EQ.1.OR.IFEHL.EQ.2) THEN 
            HK = B - XK 
            IFEHL = 0 
            CALL DESEXT (XK, YK, DEQ, N, B, HK, B - XK, EPSABS, EPSREL, &
            IFEHL, AUXF, LDF)                                           
            IF (IFEHL.EQ.1) IERR = 6 
         ENDIF 
      ELSE 
         DO 65 I = 1, N 
            YK1 (I) = YK (I) 
   65    END DO 
         CALL IRKDRV (DEQ, N, 10, IFLRKG, 9, 0, 0, IFEHL, FMACHP,       &
         EPSIVP, G, XK, B, YK1, YK, WORK1, WORK21, WORK22, WORK23)      
         IFLRKG = 0 
         IF (IFEHL.NE.0) IERR = 6 
      ENDIF 
!                                                                       
      IF (IERR.NE.0) RETURN 
!                                                                       
      CALL BDCD (YSTART, YK, N, R) 
!                                                                       
      RMAX = 0.0D0 
      DO 70 K = 1, N 
         RMAX = DMAX1 (RMAX, DABS (R (K) ) ) 
   70 END DO 
      IF (RMAX.LT.EPSBC) RETURN 
!                                                                       
!** If  ITMAX  NEWTON iterations have been performed                    
!** without finding a sufficiently good approximation                   
!** YSTART for Y(A): Return                                             
!                                                                       
      ITER = ITER + 1 
      IF (ITER.GT.ITMAX) THEN 
         IERR = 7 
         RETURN 
      ENDIF 
!                                                                       
!** Find a better approximation YSTART for                              
!** Y(A) via NEWTON's method. Form the JACOBI                           
!** matrix AMAT approximately from onesided                             
!** difference quotients                                                
!                                                                       
      DO 100 JACOBI = 1, N 
         DO 80 I = 1, N 
            YK (I) = YSTART (I) 
            YAJ (I) = YK (I) 
   80    END DO 
         IF (DABS (YK (JACOBI) ) .LT.EPS2) THEN 
            YK (JACOBI) = YK (JACOBI) + EPS1 
            DELTA = 1.0D0 / EPS1 
         ELSE 
            YK (JACOBI) = YK (JACOBI) * (1.0D0 + EPS1) 
            DELTA = 1.0D0 / (EPS1 * YK (JACOBI) ) 
         ENDIF 
         YAJ (JACOBI) = YK (JACOBI) 
         XK = A 
         HK = H 
!                                                                       
         IF (IPROC.EQ.1) THEN 
            CALL IVP (XK, HK, YK, N, DEQ, B, EPSABS, EPSREL, 1, IFMAX,  &
            IFA, IFEHL)                                                 
            IF (IFEHL.EQ.5) IERR = 6 
         ELSEIF (IPROC.EQ.2) THEN 
            IFEHL = 0 
            CALL DESABM (XK, YK, DEQ, N, B, HK, B - A, EPSABS, EPSREL,  &
            IND, IFEHL, AUXF, LDF)                                      
            IF (IFEHL.EQ.2.OR.IFEHL.EQ.3) THEN 
               HK = B - XK 
               IFEHL = 0 
               CALL DESABM (XK, YK, DEQ, N, B, HK, B - XK, EPSABS,      &
               EPSREL, IND, IFEHL, AUXF, LDF)                           
               IF (IFEHL.EQ.2) IERR = 6 
            ENDIF 
         ELSEIF (IPROC.EQ.3) THEN 
            CALL RKTRB (XK, B, N, DEQ, YK, EPSABS, EPSREL, IFLARK,      &
            WOR1RK, WOR2RK, IWORRK, IFEHL)                              
            WOR1RK (4) = H 
            IF (IFEHL.EQ. - 2) IERR = 6 
         ELSEIF (IPROC.EQ.4) THEN 
            IFEHL = 0 
            CALL DESEXT (XK, YK, DEQ, N, B, HK, B - A, EPSABS, EPSREL,  &
            IFEHL, AUXF, LDF)                                           
            IF (IFEHL.EQ.1.OR.IFEHL.EQ.2) THEN 
               HK = B - XK 
               IFEHL = 0 
               CALL DESEXT (XK, YK, DEQ, N, B, HK, B - XK, EPSABS,      &
               EPSREL, IFEHL, AUXF, LDF)                                
               IF (IFEHL.EQ.1) IERR = 6 
            ENDIF 
         ELSE 
            DO 85 I = 1, N 
               YK1 (I) = YK (I) 
   85       END DO 
            CALL IRKDRV (DEQ, N, 10, IFLRKG, 9, 0, 0, IFEHL, FMACHP,    &
            EPSIVP, G, XK, B, YK1, YK, WORK1, WORK21, WORK22, WORK23)   
            IF (IFEHL.NE.0) IERR = 6 
         ENDIF 
!                                                                       
         IF (IERR.NE.0) RETURN 
!                                                                       
         CALL BDCD (YAJ, YK, N, RJ) 
!                                                                       
         DO 90 K = 1, N 
            AMAT (K, JACOBI) = (RJ (K) - R (K) ) * DELTA 
   90    END DO 
  100 END DO 
!                                                                       
      CALL GAUSS (N, AMAT, 100, R, D, IFLAG, SCAL, IPIVOT) 
!                                                                       
!** Return if the JACOBI matrix is singular                             
!                                                                       
      IF (IFLAG.EQ.0) THEN 
         IERR = 8 
         RETURN 
      ENDIF 
!                                                                       
      DO 110 I = 1, N 
         YSTART (I) = YSTART (I) - D (I) 
  110 END DO 
      GOTO 50 
      END SUBROUTINE BVP                            
