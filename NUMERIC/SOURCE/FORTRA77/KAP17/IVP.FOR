C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 17.3.7.2}
C[  {Adaptive Methods for Initial Value Problems}
C[  {Automatic Step Size Control,
C[   Adaptive Methods for Initial Value Problems}*)
C[LE*)
      SUBROUTINE  IVP (XK, HK, YK, N, DES, XE, EPSABS, 
     1                 EPSREL, INDEX, NMAX, NUSED, IERR)
C[IX{IVP}*)
C
C*****************************************************************
C                                                                *
C  For an approximation YK of the solution Y at XK of a system   *
C  of ordinary differential equations of 1st order               *
C                    Y' = F(X,Y),                                *
C  this program computes an approximation for the solution Y     *
C  at XE.                                                        *
C  We use step size control in such a way that the error of the  *
C  computed approximation falls either absolutely or relatively  *
C  within the given error bounds EPSABS or EPSREL.               *
C  In case of the Prince-Dormand embedding formula we also check *
C  for stiffness of the system.                                  *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C XK     - initial value of the independent variable X           *
C HK     - proposed step size for the next step                  *
C YK     - vector YK(1:N); value of the solution of the          *
C          differential equation at XK                           *
C N      - number of differential equations ( 1 <= N <= 20 )     *
C DES    - right-hand side of the differential equation which    *
C          must be provided as a SUBROUTINE in the form:         *
C             SUBROUTINE DES (X, Y, F, N)                        *
C          (starting with:  DOUBLE PRECISION Y(N), F(N), etc.).  *
C          Here F represents the value of the differential       *
C          equation's right-hand side at (X,Y). (DES has to be   *
C          defined as EXTERNAL in the calling program)           *
C XE     - location where the solution is desired; XE may not    *
C          be chosen smaller than XK.                            *
C EPSABS - error bound for the absolute accuracy of the desired  *
C          solution. EPSABS has to be >= 0; if EPSABS = 0, only  *
C          the relative accuracy is considered.                  *
C EPSREL - error bound for the relative accuracy of the desired  *
C          solution. EPSREL has to be >= 0; if EPSREL = 0,       *
C          only the absolute accuracy is considered.             *
C INDEX  - chooses the embedding formula with step size control: *
C             = 0: RUNGE-KUTTA method  2nd/3rd order             *
C             =-1: Prince-Dormand embedding formula of 4th/5th   *
C                  order (checking for stiffness of the system   *
C                  of the DEs here, see output IERR)             *       
C            else: formula of ENGLAND  4th/5th order             *
C NMAX   - upper limit for the number of function evaluations    *
C          allowed for the right-hand side F.                    *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C XK     - location that was reached during last integration.    *
C          If IERR = 0 normally XK is equal to XE                *
C HK     - local step size used last ( it should remain un-      *
C          changed for the next step )                           *
C YK     - approximate value of the solution Y at the new        *
C          location XK                                           *
C NUSED  - number of function evaluations actually used          *
C IERR   - error parameter:                                      *
C          = 0: everything o.k.                                  *
C               If Prince-Dormand was specified, the sytsem was  *
C               not found to be stiff.                           *
C          = 1: both error bounds  EPS...  are too small         *
C               (relative to the machine constant)               *
C          = 2: XE <= XK       (wrt machine constant)            *
C          = 3: step size  HK <= 0   (wrt machine constant)      *
C          = 4: N > 20  or  N <= 0                               *
C          = 5: NUSED > NMAX:  the number of allowed functional  *
C               evaluations is insufficient to determine an      *
C               adequate approximate solution with the required  *
C               accuracy; on termination, XK and HK contain      *
C               the current values                               *
C          =-1: The computations terminated ok, but the Prince-  *
C               Dormand formula has detected possible stiffness. *
C          =-2: The computations have terminated ok, but the     *
C               Prince-Dormand formula has recognized the system *
C               as stiff using two criteria: we recommend to use *
C               a method suited for stiff DEs instead.           *   
C          =-3: NUSED > NMAX:  The number of allowed functional  *
C               evaluations does not suffice, see IERR = 5.      *
C               Moreover, the Prince-Dormand formula seems to    *
C               indicate that the system is stiff; we recommend  *
C               to use a suitable stiff DE solver.               *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Required subroutines: DVNORM, RUKU23, ENGL45, MACHPD, PRDO45  *
C                                                                *
C*****************************************************************
C                                                                *
C  Author      : Klaus Niederdrenk                               *
C  Date        : 11.01.1985 / 4.1.1995                           *
C  Source code : FORTRAN 77                                      *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION YK(N)
      DOUBLE PRECISION Y(20), YT(20), Y00(20)
      LOGICAL IEND, LEPSI, STEIF1, STEIF2, STEIFA
      EXTERNAL DES
      SAVE FMACHP, EPS1, EPS2, LEPSI
      DATA  LEPSI / .TRUE. / , Y00 / 20 * 0.0D0 /
C
C** FMACHP is the machine constant
C** (i.e., the smallest positive machine number for which  1 + FMACHP > 1)
C** EPS1 bounds the admissable initial step size HK below.
C** EPS2 is used to determine a level for negligable quantities.
C
      IF ( LEPSI ) THEN
        FMACHP = 1.0D0
  10    FMACHP = 0.5D0 * FMACHP
        IF (MACHPD(1.0D0 + FMACHP) .EQ. 1) GOTO 10
        FMACHP = 2.0D0 * FMACHP
        EPS1 = FMACHP ** 0.75D0
        EPS2 = 100.0D0 * FMACHP
        LEPSI = .FALSE.
      ENDIF
C
C** Preassign local variables
C
      SG = DSIGN (1.0D0, XE)
      XEND = (1.0D0 - SG * EPS2) * XE
      IERR = 0
      NUSED = 0
      IEND = .FALSE.
      STEIF1 = .FALSE.
      STEIF2 = .FALSE.
      STEIFA = .FALSE.
      IANZ = 0
C
C** Check input parameters
C
      YMAX = DVNORM(YK, Y00, N)
      IF (EPSABS .LE. EPS2*YMAX .AND. EPSREL .LE. EPS2) THEN
        IERR = 1
      ELSE IF (XEND .LT. XK) THEN
        IERR = 2
      ELSE IF (HK .LT. EPS2 * DABS(XK)) THEN
        IERR = 3
      ELSE IF (N .LE. 0 .OR. N .GT. 20) THEN
        IERR = 4
      ENDIF
      IF (IERR .NE. 0)  R E T U R N
C
C**********  C O N T R O L L I N G   a l g o r i t h  **********
C
       IF (XK+HK .GT. XEND) THEN
         HK = XE - XK
         DUMMY = HK
         IEND = .TRUE.
       ENDIF
C
C** Integrate on the interval *[*XK, XE*]* in suitable steps
C
  50  CONTINUE
C
C** Call the desired one step method
C
      IF (INDEX .EQ. 0) THEN
         CALL RUKU23 (XK, HK, YK, N, DES, Y, YT)
         NUSED=NUSED+3
      ELSE IF (INDEX .EQ. -1) THEN
         CALL PRDO45 (XK, HK, YK, N, DES, Y, YT, STEIF1, STEIFA)
         NUSED = NUSED+7 
         IF (STEIFA) THEN
           IANZ = IANZ+1
           IF (IANZ .GE. 3) STEIF2 = .TRUE.
         ELSE
           IANZ = 0
         ENDIF
      ELSE
         CALL ENGL45 (XK, HK, YK, N, DES, Y, YT)
         NUSED=NUSED+6
      ENDIF
C
      DIFF = DVNORM(Y, YT, N)
C
      IF (DIFF .LT. EPS2) THEN
         S = 2.0D0
      ELSE
         YMAX = DVNORM(YT, Y00, N)
         S = DSQRT(HK * (EPSABS + EPSREL*YMAX) / DIFF)
         IF (INDEX .NE. 0)   S = DSQRT(S)
      ENDIF
C        
      IF (S .GT. 1.0D0) THEN
C
C** accept the performed integration with step size HK
C
        DO 60 I = 1, N
          YK(I) = YT(I)
   60   CONTINUE
C
        XK = XK + HK
C
C** if then endpoint XE has been reached or if more than the 
C** allowable function evaluations were used: go back
C
  70    IF (IEND) THEN
          HK = DUMMY
          IF (INDEX .EQ. -1) THEN 
            IF (STEIF1 .OR.  STEIF2) IERR = -1
            IF (STEIF1 .AND. STEIF2) IERR = -2
          ENDIF
          R E T U R N
        ELSE IF (NUSED .GT. NMAX) THEN
          IERR = 5
          IF (INDEX .EQ. -1 .AND. (STEIF1 .OR. STEIF2)) IERR = -3
          R E T U R N
        ENDIF
C
C** increase the step size for the next step maximally by a factor of two
C
        HK = HK * DMIN1(2.0D0, 0.98D0*S)
C
        IF ((XK+HK) .GE. XEND) THEN
          DUMMY = HK
          HK = XE - XK
          IEND = .TRUE.
C
C** if very close to XE: go back
C
          IF (HK .LT. EPS1 * DABS(XE)) GOTO 70
        ENDIF
      ELSE
C
C** the previous step is unaccaptable, the step size HK must be decreased,
C** at most it must be halved
C
        HK = HK * DMAX1(0.5D0, 0.98D0*S)
        IEND = .FALSE.
      ENDIF
C
      GOTO 50
C
      END
C
C
C[BA*)
C[LE*)
      DOUBLE PRECISION FUNCTION DVNORM (F1, F2, N)
C[IX{DVNORM}*)
C
C*****************************************************************
C                                                                *
C  This program finds the maximum norm of the difference F1 - F2 *
C  of two vectors F1 and F2 of length N.                         *
C                                                                *
C----------------------------------------------------------------*
C[BE*)
C                                                                *
C  required subroutines: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  Author      : Klaus Niederdrenk                               *
C  Date        : 11.01.1985                                      *
C  Source code : FORTRAN 77                                      *
C                                                                *
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION F1(N), F2(N)
      DVNORM=0.0D0
      DO 10 I = 1, N
        DVNORM = DMAX1 (DVNORM, DABS(F1(I)-F2(I)))
   10 CONTINUE
      R E T U R N
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE RUKU23 (X, H, Y, N, DES, Y2, Y3)
C[IX{RUKU23}*)
C
C[BE*)
C*****************************************************************
C                                                                *
C[BA*)
C*****************************************************************
C                                                                *
C  Starting with an approximation Y at X, this program computes  * 
C  approximations Y2 and Y3 at X + H using the RUNGE-KUTTA       *
C  embedding formula of 2nd and 3rd order for a system of        *
C  ordinary differential equations of 1st order                  *
C                  Y' = F(X,Y).                                  *
C  The system of N ordinary differential equations of 1st order  *
C  must be provided by a SUBROUTINE DES.                         *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  X     - initial value for the independent variable X          *
C  H     - step size                                             *
C  Y     - vector Y(1:N); value for the solution of the          *
C          differential equation at X                            *
C  N     - number of differential equations  ( 1 <= N <= 20 )    *
C  DES   - right hand side of the differential equation. It has  *
C          provided by the user as a SUBROUTINE in the form:     *
C             SUBROUTINE DES (X, Y, N, F)                        *
C          (starting with: DOUBLE PRECISION Y(N), F(N), etc. ).  *
C          Here F denotes the right hand side of the differential*
C          equation at (X,Y). ( In the calling program DES has to*
C          declared as EXTERNAL)                                 *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  Y2    - vector Y2(1:N); approximate solution using the 2nd    *
C          order method for solving the differential equation    *
C          at X+H                                                *
C  Y3    - vector Y3(1:N); approximate solution using the 3rd    *
C          order method for solving the differential equation    *
C          at X+H                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required : none                                   *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Klaus Niederdrenk                                  *
C  date     : 11.01.1985                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION Y(N), Y2(N), Y3(N), DUMMY(20)
      DOUBLE PRECISION K1(20), K2(20), K3(20)
C
      CALL DES (X, Y, N, K1)
      DO 10 I = 1, N
        DUMMY(I)=Y(I)+H*K1(I)
  10  CONTINUE
      CALL DES (X+H, DUMMY, N, K2)
      DO 20 I = 1, N
        DUMMY(I)=Y(I)+0.25D0*H*(K1(I)+K2(I))
  20  CONTINUE
      CALL DES (X+0.5D0*H, DUMMY, N, K3)
C
      DO 100 I = 1, N
        Y2(I)=Y(I)+0.5D0*H*(K1(I)+K2(I))
        Y3(I)=Y(I)+H/6.0D0*(K1(I)+K2(I)+4.0D0*K3(I))
 100  CONTINUE
      R E T U R N
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE ENGL45 (X, H, Y, N, DES, Y4, Y5)
C[IX{ENGL45}*)
C
C*****************************************************************
C                                                                *
C  Starting from an approximation Y at X, this program uses the  *
C  ENGLAND embedding formula of 4th and 5th order to find appro- *
C  ximations Y4 and Y5 at X + H for the solution of the system of*
C  differential equations of 1st order                           *
C                  Y' = F(X,Y) .                                 *
C  The system contains N ordinary differential equations of 1st  *
C  order, with the right hand side F(X,Y) provided by the user   *
C  in a SUBROUTINE DES.                                          *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  X     - initial value for the independent variable X          *
C  H     - step size                                             *
C  Y     - vector Y(1:N); value for the solution of the          *
C          differential equation at X                            *
C  N     - number of differential equations  ( 1 <= N <= 20 )    *
C  DES   - right hand side of the differential equation. It has  *
C          provided by the user as a SUBROUTINE in the form:     *
C             SUBROUTINE DES (X, Y, N, F)                        *
C          (starting with: DOUBLE PRECISION Y(N), F(N), etc. ).  *
C          Here F denotes the right hand side of the differential*
C          equation at (X,Y). ( In the calling program DES has to*
C          declared as EXTERNAL)                                 *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  Y4    - vector Y4(1:N); approximate solution using the 4th    *
C          order method for solving the differential equation    *
C          at X+H                                                *
C  Y5    - vector Y5(1:N); approximate solution using the 5th    *
C          order method for solving the differential equation    *
C          at X+H                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required : none                                   *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Klaus Niederdrenk                                  *
C  date     : 11.01.1985                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION Y(N), Y4(N), Y5(N)
      DOUBLE PRECISION DUMMY(20)
      DOUBLE PRECISION K1(20),K2(20),K3(20),K4(20),K5(20),K6(20)
C
      CALL DES (X, Y, N, K1)
      DO 10  I = 1, N
        DUMMY(I)=Y(I)+0.5D0*H*K1(I)
  10  CONTINUE
      CALL DES (X+0.5D0*H, DUMMY, N, K2)
C
      DO 20  I = 1, N
        DUMMY(I)=Y(I)+0.25D0*H*(K1(I)+K2(I))
  20  CONTINUE
      CALL DES (X+0.5D0*H, DUMMY, N, K3)
C
      DO 30  I = 1, N
        DUMMY(I)=Y(I)+H*(-K2(I)+2.0D0*K3(I))
  30  CONTINUE
      CALL DES (X+H, DUMMY, N, K4)
C
      DO 40  I = 1, N
        DUMMY(I)=Y(I)+H/27.0D0*(7.0D0*K1(I)+10.0D0*K2(I)+K4(I))
  40  CONTINUE
      CALL DES (X+2.0D0/3.0D0*H, DUMMY, N, K5)
C
      DO 50  I = 1, N
        DUMMY(I)=Y(I)+0.16D-02*H*(28.0D0*K1(I)-125.0D0*K2(I)+
     +           546.0D0*K3(I)+54.0D0*K4(I)-378.0D0*K5(I))
  50  CONTINUE
      CALL DES (X+0.2D0*H, DUMMY, N, K6)
C
      DO 100  I = 1, N
        Y4(I)=Y(I)+H/6.0D0*(K1(I)+4.0D0*K3(I)+K4(I))
        Y5(I)=Y(I)+H/336.0D0*(14.0D0*K1(I)+35.0D0*K4(I)
     +            +162.0D0*K5(I)+125.0D0*K6(I))
 100  CONTINUE
      R E T U R N
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE PRDO45 (X, H, Y, N, DES, Y4, Y5, ST1, ST2)
C[IX{PRDO45}*)
C
C*****************************************************************
C                                                                *
C  Starting from an approximation Y at X, this program uses the  *
C  Prince-Dormand embedding formula of 4th and 5th order to find *
C  apprximations Y4 and Y5 at X + H for the solution of the      *
C  system of differential equations of 1st order                 *
C                  Y' = F(X,Y) .                                 *
C  The system contains N ordinary differential equations of 1st  *
C  order, with the right hand side F(X,Y) provided by the user   *
C  in a SUBROUTINE DES.                                          *
C  This program tests for stiffness in two ways.                 *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  X     - initial value for the independent variable X          *
C  H     - step size                                             *
C  Y     - vector Y(1:N); value for the solution of the          *
C          differential equation at X                            *
C  N     - number of differential equations  ( 1 <= N <= 20 )    *
C  DES   - right hand side of the differential equation. It has  *
C          provided by the user as a SUBROUTINE in the form:     *
C             SUBROUTINE DES (X, Y, N, F)                        *
C          (starting with: DOUBLE PRECISION Y(N), F(N), etc. ).  *
C          Here F denotes the right hand side of the differential*
C          equation at (X,Y). ( In the calling program DES has to*
C          be declared as EXTERNAL)                              *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  Y4    - vector Y4(1:N); approximate solution using the 4th    *
C          order method for solving the differential equation    *
C          at X+H                                                *
C  Y5    - vector Y5(1:N); approximate solution using the 5th    *
C          order method for solving the differential equation    *
C          at X+H                                                *
C  ST1   - logical variable: .TRUE. , if the test for stiffness  *
C          (using the approximate dominat eigenvalue) was        *
C          positive ; else the value of ST! remains unaltered.   * 
C  ST2   - logical variable: .TRUE. , if the test of stiffness   *
C          due to  Shampine  and  Hiebert  is positive; else the *
C          value of ST2 is set to  .FALSE. .                     *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required : DVNORM                                 *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Klaus Niederdrenk                                  *
C  date     : 1.4.1995                                           *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION Y(N), Y4(N), Y5(N)
      DOUBLE PRECISION DUMMY(20)
      DOUBLE PRECISION K1(20),K2(20),K3(20),K4(20),K5(20),K6(20),
     #                 K7(20),G6(20),G7(20)
      LOGICAL ST1, ST2
C
      CALL DES (X, Y, N, K1)
      DO 10  I = 1, N
        DUMMY(I)=Y(I)+0.2D0*H*K1(I)
  10  CONTINUE
      CALL DES (X+0.2D0*H, DUMMY, N, K2)
C
      DO 20  I = 1, N
        DUMMY(I)=Y(I)+0.075D0*H*(K1(I)+3.0D0*K2(I))
  20  CONTINUE
      CALL DES (X+0.3D0*H, DUMMY, N, K3)
C
      DO 30  I = 1, N
        DUMMY(I)=Y(I)+H/45.0D0*(44.0D0*K1(I)-168.0D0*K2(I)
     +           +160.0D0*K3(I))
  30  CONTINUE
      CALL DES (X+0.8D0*H, DUMMY, N, K4)
C
      DO 40  I = 1, N
        DUMMY(I)=Y(I)+H/6561.0D0*(19372.0D0*K1(I)-76080.0D0*K2(I)
     +           +64448.0D0*K3(I)-1908.0D0*K4(I))
  40  CONTINUE
      CALL DES (X+8.0D0/9.0D0*H, DUMMY, N, K5)
C
      DO 50  I = 1, N
        G6(I)=Y(I)+H/167904.0D0*(477901.0D0*K1(I)-1806240.0D0*K2(I)
     +        +1495424.0D0*K3(I)+46746.0D0*K4(I)-45927.0D0*K5(I))
  50  CONTINUE
      CALL DES (X+H, G6, N, K6)
C
      DO 60 I = 1, N
        G7(I)=Y(I)+H/142464.0D0*(12985.0D0*K1(I)+64000.0D0*K3(I) 
     +        +92750.0D0*K4(I)-45927.0D0*K5(I)+18656.0D0*K6(I))
  60  CONTINUE
      CALL DES (X+H, G7, N, K7)
C
      DO 100  I = 1, N
        Y5(I)=G7(I)
        Y4(I)=Y(I)+H/21369600.0D0*(1921409.0D0*K1(I)+9690880.0D0
     +        *K3(I)+13122270.0D0*K4(I)-5802111.0D0*K5(I)+
     +        1902912.0D0*K6(I)+534240.0D0*K7(I))
 100  CONTINUE
C
C**  Test for stiffness via dominant eigenvalue (approximately)
C
      IF (DVNORM(K7,K6,N) .GT. 3.3*DVNORM(G7,G6,N)) ST1 = .TRUE.
C
C**  One step of testing stiffness according to Shampine and Hiebert
C
      DO 110 I = 1, N
        G6(I) = H*(2.2D0*K2(I)+0.13D0*K4(I)+0.144D0*K5(I))
        G7(I) = H*(2.134D0*K1(I)+0.24D0*K3(I)+0.1D0*K6(I))
 110  CONTINUE
      IF (DVNORM(G6,G7,N) .LT. DVNORM(Y4,Y5,N)) THEN
        ST2 = .TRUE.
      ELSE
        ST2 = .FALSE.
      ENDIF
      R E T U R N
      END
