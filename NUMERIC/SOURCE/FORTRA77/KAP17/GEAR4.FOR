C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 17.7.3}
C[  {Gear's method for integrating stiff systems of DEs}
C[  {Gear's method for integrating stiff systems of DEs}*)
C[LE*)
      SUBROUTINE GEAR4 (XK, HK, YK, N, DES, XE, EPSABS,
     1                 EPSREL, NMAX, NUSED, IERR)
C[IX{GEAR4}*)
C
C*****************************************************************
C                                                                *
C  Starting from an approximation YK for the solution Y of a     *
C  system of ordinary differential equations of first order      *
C                    Y' = F(X,Y)                                 *
C  at XK, this program computes an approximate solution YE at XE.*
C  Here we compute internally with a step size control, that     *
C  ensures the error of the computed solution to be less than the*
C  given absolute or relative error bounds EPSABS and EPSREL.    *
C  These bounds must be specified small enough for good results. *
C  The method used is the multistep method of Gear of fourth     *
C  order which is highly capable of solving stiff DEs. (Stiff DEs*
C  are those DE systems which have solution components of very   *
C  disparate growths.)                                           *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  XK    - starting value for X                                  *
C  HK    - proposed step size for first step                     *
C  YK    - vector YK(1:N); Y value of the solution to the DE     *
C          at XK                                                 *
C  N     - number of DEs  ( 1 <= N <= 20 )                       *
C  DES   - right hand side of the DE, given as a subroutine:     *
C             SUBROUTINE DES (X, Y, N, F)                        *
C          ( starting with: DOUBLE PRECISION Y(N), F(N), X ).    *
C          Here F is the value of the right hand side at (X,Y).  *
C          ( DES must be declared as EXTERNAL in the calling     *
C            program.)                                           *
C  XE    - X value for desired solution; XE > XK.                *
C EPSABS - error bound for absolute error; EPSABS >= 0; if       *
C          EPSABS = 0 the algorithm maintains only the relative  *
C          accuracy.                                             *
C EPSREL - error bound for relative error; EPSREL >= 0; if       *
C          EPSREL = 0 the algorithm maintains only the absolute  *
C          accuracy.                                             *
C  NMAX  - maximal number of evaluations of the right hand side  *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  XK    - final X value of the integration. If IERR = 0,        *
C          usually XK = XE (within machine precision).           *
C  HK    - terminal step size used; should be used for subsequent*
C          integrations                                          *
C  YK    - approximate value for the solution at XK              *
C  NUSED - number of actual evaluations of the right hand side   *
C  IERR  - error parameter:                                      *
C          = 0: all is o.k.                                      *
C          = 1: both error bounds  EPS...  too small             *
C                              (relative to the machine constant)*
C          = 2: XE <= XK       (relative to the machine constant)*
C          = 3: step size  HK <= 0  (rel. to machine precision)  *
C          = 4: N > 20   or   N <= 0                             *
C          = 5: NUSED > NMAX: Number of allowed function         *
C               evaluations was exceeded; try to restart with    *
C               XK, YK and HK.                                   *
C          = 6: The Jacobi matrix is singular; XK, YK, HK contain*
C               the values reached.                              *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Subroutines used:  IVP, GAUSSP, GAUSSS, DVNORM, MACHPD        *
C                                                                *
C*****************************************************************
C                                                                *
C  Author      : Klaus Niederdrenk                               *
C  Date        : 1.22.1996                                       *
C  Source code : FORTRAN 77                                      *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION YK(1:N)
      PARAMETER ( NDGL = 20 )
      DIMENSION ZJ(0:4,1:NDGL), ZJP1(0:4,1:NDGL), F(1:NDGL)
      DIMENSION FS(1:NDGL,1:NDGL), HELP(1:NDGL), Y0(1:NDGL)
      DIMENSION YKP1(1:NDGL), CON(1:NDGL)
      DIMENSION D(1:NDGL), IPIVOT(1:NDGL), FSG(1:NDGL,1:NDGL)
      LOGICAL IEND, LEPSI
      EXTERNAL DES
      SAVE EPS1, EPS2, LEPSI, Y0, HS
      DATA LEPSI / .TRUE. / , Y0 / NDGL * 0.0D0 /
C
C** Using the machine constant FMACHP, we determine EPS1 in order
C** to avoid excessively small final steps near XE, EPS2 to check
C** for zero and HS as the optimal step size for approximating the
C** Jacobi matrix. (This is done only once at the start.)
C
      IF ( LEPSI ) THEN
        FMACHP = 1.0D0
  10    FMACHP = 0.5 * FMACHP
        IF ( MACHPD(1.0 + FMACHP) .EQ. 1) GOTO 10
        FMACHP = 2.0D0 * FMACHP
        EPS1   = FMACHP ** 0.75D0
        EPS2   = 100.0D0 * FMACHP
        HS     = 10.0D0 * SQRT(FMACHP)
        LEPSI  = .FALSE.
      ENDIF
C
C** Initialize
C
      SG    = DSIGN(1.0D0, XE)
      XEND  = (1.0D0 - SG*EPS2) * XE
      IERR  = 0
      NUSED = 0
      IEND  = .FALSE.
C
C** Check input parameters
C
      YMAX = DVNORM(YK, Y0, N)
      IF (EPSABS .LE. EPS2*YMAX .AND. EPSREL .LE. EPS2) THEN
        IERR = 1
      ELSE IF (XEND .LT. XK) THEN
        IERR = 2
      ELSE IF (HK .LT. EPS2*DABS(XK)) THEN
        IERR = 3
      ELSE IF ( N .LE. 0 .OR. N .GT. NDGL ) THEN
        IERR = 4
      ENDIF
      IF (IERR .NE. 0)  R E T U R N
C
C****  compute first integration   ****
C
      IF (XK+HK .GT. XEND) THEN
        HK = XE - XK
        DUMMY = HK
        IEND = .TRUE.
      ENDIF
      DO 20 I = 1, N
        HELP(I) = YK(I)
  20  CONTINUE
      XKA = XK
      XKE = XKA
      HKA = 0.25*HK
      HK1 = HKA
      DO 40 K = 1, 4
        XKE = XKE + HKA
        CALL IVP (XKA, HK1, HELP, N, DES, XKE, EPSABS, EPSREL,
     #            1, NMAX-NUSED, NANL, IERR)
        NUSED = NUSED + NANL
        IF( IERR .NE. 0 ) R E T U R N
        DO 30 I = 1, N
          ZJP1(K,I) = HELP(I)
  30    CONTINUE
  40  CONTINUE
      CALL DES (XK, YK, N, F)
      NUSED = NUSED + 1
C
C** Determine first Gear-Nordsieck approximation
C
      DO 50 I = 1, N
        ZJ(0,I) = YK(I)
        ZJ(1,I) = HK*F(I)
        ZJ(2,I) = 1.0D0/24.0D0*(35.0D0*YK(I) - 104.0D0*ZJP1(1,I)
     #            + 114.0D0*ZJP1(2,I) - 56.0D0*ZJP1(3,I)
     #            +  11.0D0*ZJP1(4,I))
        ZJ(3,I) = 1.0D0/12.0D0*(-5.0D0*YK(I) + 18.0D0*ZJP1(1,I)
     #            - 24.0D0*ZJP1(2,I) + 14.0D0*ZJP1(3,I)
     #            -  3.0D0*ZJP1(4,I))
        ZJ(4,I) = 1.0D0/24.0D0*(YK(I) - 4.0D0*ZJP1(1,I)
     #            + 6.0D0*ZJP1(2,I) - 4.0D0*ZJP1(3,I) + ZJP1(4,I))
  50  CONTINUE
C
C
C****  S t e p  S i z e  A l g o r i t h m   ****
C
C
  75  CONTINUE
C
C** Compute implicit approximation using Newton method
C
      DO 90 I = 1, N
        YKP1(I) = ZJ(0,I)+ZJ(1,I)+ZJ(2,I)+ZJ(3,I)+ZJ(4,I)
  90  CONTINUE
      CALL DES (XK+HK, YKP1, N, F)
      DO 120 K = 1, N
        DO 100 I = 1, N
          HELP(I) = YKP1(I)
 100    CONTINUE
        HELP(K) = HELP(K) - HS
        CALL DES (XK+HK, HELP, N, FS(1,K))
        DO 110 I = 1, N
          FS(I,K) = -HK * 0.48D0 * (F(I) - FS(I,K)) / HS
 110    CONTINUE
        FS(K,K) = FS(K,K) + 1.0D0
 120  CONTINUE
      NUSED = NUSED + N + 1
      DO 190 I = 1, N
        CON(I) = YKP1(I) - 0.48D0 * ( ZJ(1,I) + 2.0D0*ZJ(2,I)
     #           + 3.0D0*ZJ(3,I) + 4.0D0*ZJ(4,I) )
        DO 180 K = 1, N
          FSG(K,I) = FS(K,I)
 180    CONTINUE
 190  CONTINUE
      CALL GAUSSP (N, FSG, NDGL, IPIVOT, MARK, D)
      IF ( MARK .EQ. 0 ) THEN
        IERR = 6
        R E T U R N
      ENDIF
      DO 220 ITER = 1, 3
        DO 210 I = 1, N
          HELP(I) = -YKP1(I)
          DO 200 K = 1, N
            HELP(I) = HELP(I) + FS(I,K)*YKP1(K)
 200      CONTINUE
          HELP(I) = HK*0.48D0*F(I) + HELP(I) + CON(I)
 210    CONTINUE
        CALL GAUSSS (N, FSG, NDGL, IPIVOT, HELP, YKP1)
        CALL DES (XK+HK, YKP1, N, F)
 220  CONTINUE
      NUSED = NUSED + 3
C
C** Determine corresponding Gear-Nordsieck approximation
C
      DO 230 I = 1, N
        HELP(I) = HK*F(I) - ZJ(1,I) - 2.0D0*ZJ(2,I)
     #              - 3.0D0*ZJ(3,I) - 4.0D0*ZJ(4,I)
 230  CONTINUE
      DO 250 I = 1, N
        ZJP1(0,I) = YKP1(I)
        ZJP1(1,I) = HK*F(I)
        ZJP1(2,I) = ZJ(2,I) + 3.0D0*ZJ(3,I) + 6.0D0*ZJ(4,I)
     #              + 0.7D0*HELP(I)
        ZJP1(3,I) = ZJ(3,I) + 4.0D0*ZJ(4,I) + 0.2D0*HELP(I)
        ZJP1(4,I) = ZJ(4,I) + 0.02D0*HELP(I)
 250  CONTINUE
C
C** Determine whether the last step should be accepted
C
      DO 260 I = 1, N
        HELP(I) = ZJP1(4,I)
        CON(I)  = ZJ(4,I)
 260  CONTINUE
      DIFF = DVNORM(HELP, CON, N)
      YMAX = DVNORM(YKP1, Y0, N)
      EPS = (EPSABS + EPSREL*YMAX) / 6.0D0
      Q = DSQRT(DSQRT(EPS/DIFF))/1.2
      IF ( DIFF .LT. EPS ) THEN
C
C** Accept last step; prepare for next integration step
C
        XK = XK + HK
        DO 270 I = 1, N
          YK(I) = YKP1(I)
 270    CONTINUE
C
C** Jump back if the interval endpoint XE has been reached or
C** if the right hand side has been called too often.
C
 275    IF ( IEND ) THEN
          HK = DUMMY
          R E T U R N
        ELSE IF ( NUSED .GT. NMAX ) THEN
          IERR = 5
          R E T U R N
        ENDIF
C
C** adapt step size for next step
C
        HALT = HK
        HK = DMIN1(Q, 2.0D0) * HK
        IF ( XK + HK .GE. XEND ) THEN
          DUMMY = HK
          HK = XE - XK
          IEND = .TRUE.
C
C** jump back if sufficiently close to XE
C
          IF ( HK .LT. EPS1*DABS(XE) ) GOTO 275
        ENDIF
C
C** Set up the Gera-Nordsieck approximation for the next
C** integration
C
        QUOT1 = HK / HALT
        QUOT2 = QUOT1 * QUOT1
        QUOT3 = QUOT2 * QUOT1
        QUOT4 = QUOT3 * QUOT1
        DO 280 I = 1, N
          ZJ(0,I) = ZJP1(0,I)
          ZJ(1,I) = QUOT1 * ZJP1(1,I)
          ZJ(2,I) = QUOT2 * ZJP1(2,I)
          ZJ(3,I) = QUOT3 * ZJP1(3,I)
          ZJ(4,I) = QUOT4 * ZJP1(4,I)
 280    CONTINUE
      ELSE
C
C** Repeat last step for a smaller step size
C** and modify the Gear-Nordsieck approximation accordingly
C
        HALT = HK
        HK = DMAX1(0.5D0, Q) * HK
        QUOT1 = HK / HALT
        QUOT2 = QUOT1 * QUOT1
        QUOT3 = QUOT2 * QUOT1
        QUOT4 = QUOT3 * QUOT1
        DO 290 I = 1, N
          ZJ(1,I) = QUOT1 * ZJ(1,I)
          ZJ(2,I) = QUOT2 * ZJ(2,I)
          ZJ(3,I) = QUOT3 * ZJ(3,I)
          ZJ(4,I) = QUOT4 * ZJ(4,I)
 290    CONTINUE
        IEND = .FALSE.
      ENDIF
C
      GOTO 75
C
      END
