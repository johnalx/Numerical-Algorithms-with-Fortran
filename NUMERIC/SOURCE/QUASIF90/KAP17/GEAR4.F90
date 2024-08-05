      SUBROUTINE GEAR4 (XK, HK, YK, N, DES, XE, EPSABS, EPSREL, NMAX,   &
      NUSED, IERR)                                                      
!                                                                       
!*****************************************************************      
!                                                                *      
!  Starting from an approximation YK for the solution Y of a     *      
!  system of ordinary differential equations of first order      *      
!                    Y' = F(X,Y)                                 *      
!  at XK, this program computes an approximate solution YE at XE.*      
!  Here we compute internally with a step size control, that     *      
!  ensures the error of the computed solution to be less than the*      
!  given absolute or relative error bounds EPSABS and EPSREL.    *      
!  These bounds must be specified small enough for good results. *      
!  The method used is the multistep method of Gear of fourth     *      
!  order which is highly capable of solving stiff DEs. (Stiff DEs*      
!  are those DE systems which have solution components of very   *      
!  disparate growths.)                                           *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  XK    - starting value for X                                  *      
!  HK    - proposed step size for first step                     *      
!  YK    - vector YK(1:N); Y value of the solution to the DE     *      
!          at XK                                                 *      
!  N     - number of DEs  ( 1 <= N <= 20 )                       *      
!  DES   - right hand side of the DE, given as a subroutine:     *      
!             SUBROUTINE DES (X, Y, N, F)                        *      
!          ( starting with: DOUBLE PRECISION Y(N), F(N), X ).    *      
!          Here F is the value of the right hand side at (X,Y).  *      
!          ( DES must be declared as EXTERNAL in the calling     *      
!            program.)                                           *      
!  XE    - X value for desired solution; XE > XK.                *      
! EPSABS - error bound for absolute error; EPSABS >= 0; if       *      
!          EPSABS = 0 the algorithm maintains only the relative  *      
!          accuracy.                                             *      
! EPSREL - error bound for relative error; EPSREL >= 0; if       *      
!          EPSREL = 0 the algorithm maintains only the absolute  *      
!          accuracy.                                             *      
!  NMAX  - maximal number of evaluations of the right hand side  *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  XK    - final X value of the integration. If IERR = 0,        *      
!          usually XK = XE (within machine precision).           *      
!  HK    - terminal step size used; should be used for subsequent*      
!          integrations                                          *      
!  YK    - approximate value for the solution at XK              *      
!  NUSED - number of actual evaluations of the right hand side   *      
!  IERR  - error parameter:                                      *      
!          = 0: all is o.k.                                      *      
!          = 1: both error bounds  EPS...  too small             *      
!                              (relative to the machine constant)*      
!          = 2: XE <= XK       (relative to the machine constant)*      
!          = 3: step size  HK <= 0  (rel. to machine precision)  *      
!          = 4: N > 20   or   N <= 0                             *      
!          = 5: NUSED > NMAX: Number of allowed function         *      
!               evaluations was exceeded; try to restart with    *      
!               XK, YK and HK.                                   *      
!          = 6: The Jacobi matrix is singular; XK, YK, HK contain*      
!               the values reached.                              *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Subroutines used:  IVP, GAUSSP, GAUSSS, DVNORM, MACHPD        *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author      : Klaus Niederdrenk                               *      
!  Date        : 1.22.1996                                       *      
!  Source code : FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION YK (1:N) 
      PARAMETER (NDGL = 20) 
      DIMENSION ZJ (0:4, 1:NDGL), ZJP1 (0:4, 1:NDGL), F (1:NDGL) 
      DIMENSION FS (1:NDGL, 1:NDGL), HELP (1:NDGL), Y0 (1:NDGL) 
      DIMENSION YKP1 (1:NDGL), CON (1:NDGL) 
      DIMENSION D (1:NDGL), IPIVOT (1:NDGL), FSG (1:NDGL, 1:NDGL) 
      LOGICAL IEND, LEPSI 
      EXTERNAL DES 
      SAVE EPS1, EPS2, LEPSI, Y0, HS 
      DATA LEPSI / .TRUE. /, Y0 / NDGL * 0.0D0 / 
!                                                                       
!** Using the machine constant FMACHP, we determine EPS1 in order       
!** to avoid excessively small final steps near XE, EPS2 to check       
!** for zero and HS as the optimal step size for approximating the      
!** Jacobi matrix. (This is done only once at the start.)               
!                                                                       
      IF (LEPSI) THEN 
         FMACHP = 1.0D0 
   10    FMACHP = 0.5 * FMACHP 
         IF (MACHPD (1.0 + FMACHP) .EQ.1) GOTO 10 
         FMACHP = 2.0D0 * FMACHP 
         EPS1 = FMACHP**0.75D0 
         EPS2 = 100.0D0 * FMACHP 
         HS = 10.0D0 * SQRT (FMACHP) 
         LEPSI = .FALSE. 
      ENDIF 
!                                                                       
!** Initialize                                                          
!                                                                       
      SG = DSIGN (1.0D0, XE) 
      XEND = (1.0D0 - SG * EPS2) * XE 
      IERR = 0 
      NUSED = 0 
      IEND = .FALSE. 
!                                                                       
!** Check input parameters                                              
!                                                                       
      YMAX = DVNORM (YK, Y0, N) 
      IF (EPSABS.LE.EPS2 * YMAX.AND.EPSREL.LE.EPS2) THEN 
         IERR = 1 
      ELSEIF (XEND.LT.XK) THEN 
         IERR = 2 
      ELSEIF (HK.LT.EPS2 * DABS (XK) ) THEN 
         IERR = 3 
      ELSEIF (N.LE.0.OR.N.GT.NDGL) THEN 
         IERR = 4 
      ENDIF 
      IF (IERR.NE.0) RETURN 
!                                                                       
!****  compute first integration   ****                                 
!                                                                       
      IF (XK + HK.GT.XEND) THEN 
         HK = XE-XK 
         DUMMY = HK 
         IEND = .TRUE. 
      ENDIF 
      DO 20 I = 1, N 
         HELP (I) = YK (I) 
   20 END DO 
      XKA = XK 
      XKE = XKA 
      HKA = 0.25 * HK 
      HK1 = HKA 
      DO 40 K = 1, 4 
         XKE = XKE+HKA 
         CALL IVP (XKA, HK1, HELP, N, DES, XKE, EPSABS, EPSREL, 1, NMAX &
         - NUSED, NANL, IERR)                                           
         NUSED = NUSED+NANL 
         IF (IERR.NE.0) RETURN 
         DO 30 I = 1, N 
            ZJP1 (K, I) = HELP (I) 
   30    END DO 
   40 END DO 
      CALL DES (XK, YK, N, F) 
      NUSED = NUSED+1 
!                                                                       
!** Determine first Gear-Nordsieck approximation                        
!                                                                       
      DO 50 I = 1, N 
         ZJ (0, I) = YK (I) 
         ZJ (1, I) = HK * F (I) 
         ZJ (2, I) = 1.0D0 / 24.0D0 * (35.0D0 * YK (I) - 104.0D0 * ZJP1 &
         (1, I) + 114.0D0 * ZJP1 (2, I) - 56.0D0 * ZJP1 (3, I) + 11.0D0 &
         * ZJP1 (4, I) )                                                
         ZJ (3, I) = 1.0D0 / 12.0D0 * ( - 5.0D0 * YK (I) + 18.0D0 *     &
         ZJP1 (1, I) - 24.0D0 * ZJP1 (2, I) + 14.0D0 * ZJP1 (3, I)      &
         - 3.0D0 * ZJP1 (4, I) )                                        
         ZJ (4, I) = 1.0D0 / 24.0D0 * (YK (I) - 4.0D0 * ZJP1 (1, I)     &
         + 6.0D0 * ZJP1 (2, I) - 4.0D0 * ZJP1 (3, I) + ZJP1 (4, I) )    
   50 END DO 
!                                                                       
!                                                                       
!****  S t e p  S i z e  A l g o r i t h m   ****                       
!                                                                       
!                                                                       
   75 CONTINUE 
!                                                                       
!** Compute implicit approximation using Newton method                  
!                                                                       
      DO 90 I = 1, N 
         YKP1 (I) = ZJ (0, I) + ZJ (1, I) + ZJ (2, I) + ZJ (3, I)       &
         + ZJ (4, I)                                                    
   90 END DO 
      CALL DES (XK + HK, YKP1, N, F) 
      DO 120 K = 1, N 
         DO 100 I = 1, N 
            HELP (I) = YKP1 (I) 
  100    END DO 
         HELP (K) = HELP (K) - HS 
         CALL DES (XK + HK, HELP, N, FS (1, K) ) 
         DO 110 I = 1, N 
            FS (I, K) = - HK * 0.48D0 * (F (I) - FS (I, K) ) / HS 
  110    END DO 
         FS (K, K) = FS (K, K) + 1.0D0 
  120 END DO 
      NUSED = NUSED+N + 1 
      DO 190 I = 1, N 
         CON (I) = YKP1 (I) - 0.48D0 * (ZJ (1, I) + 2.0D0 * ZJ (2, I)   &
         + 3.0D0 * ZJ (3, I) + 4.0D0 * ZJ (4, I) )                      
         DO 180 K = 1, N 
            FSG (K, I) = FS (K, I) 
  180    END DO 
  190 END DO 
      CALL GAUSSP (N, FSG, NDGL, IPIVOT, MARK, D) 
      IF (MARK.EQ.0) THEN 
         IERR = 6 
         RETURN 
      ENDIF 
      DO 220 ITER = 1, 3 
         DO 210 I = 1, N 
            HELP (I) = - YKP1 (I) 
            DO 200 K = 1, N 
               HELP (I) = HELP (I) + FS (I, K) * YKP1 (K) 
  200       END DO 
            HELP (I) = HK * 0.48D0 * F (I) + HELP (I) + CON (I) 
  210    END DO 
         CALL GAUSSS (N, FSG, NDGL, IPIVOT, HELP, YKP1) 
         CALL DES (XK + HK, YKP1, N, F) 
  220 END DO 
      NUSED = NUSED+3 
!                                                                       
!** Determine corresponding Gear-Nordsieck approximation                
!                                                                       
      DO 230 I = 1, N 
         HELP (I) = HK * F (I) - ZJ (1, I) - 2.0D0 * ZJ (2, I) - 3.0D0 *&
         ZJ (3, I) - 4.0D0 * ZJ (4, I)                                  
  230 END DO 
      DO 250 I = 1, N 
         ZJP1 (0, I) = YKP1 (I) 
         ZJP1 (1, I) = HK * F (I) 
         ZJP1 (2, I) = ZJ (2, I) + 3.0D0 * ZJ (3, I) + 6.0D0 * ZJ (4, I)&
         + 0.7D0 * HELP (I)                                             
         ZJP1 (3, I) = ZJ (3, I) + 4.0D0 * ZJ (4, I) + 0.2D0 * HELP (I) 
         ZJP1 (4, I) = ZJ (4, I) + 0.02D0 * HELP (I) 
  250 END DO 
!                                                                       
!** Determine whether the last step should be accepted                  
!                                                                       
      DO 260 I = 1, N 
         HELP (I) = ZJP1 (4, I) 
         CON (I) = ZJ (4, I) 
  260 END DO 
      DIFF = DVNORM (HELP, CON, N) 
      YMAX = DVNORM (YKP1, Y0, N) 
      EPS = (EPSABS + EPSREL * YMAX) / 6.0D0 
      Q = DSQRT (DSQRT (EPS / DIFF) ) / 1.2 
      IF (DIFF.LT.EPS) THEN 
!                                                                       
!** Accept last step; prepare for next integration step                 
!                                                                       
         XK = XK + HK 
         DO 270 I = 1, N 
            YK (I) = YKP1 (I) 
  270    END DO 
!                                                                       
!** Jump back if the interval endpoint XE has been reached or           
!** if the right hand side has been called too often.                   
!                                                                       
  275    IF (IEND) THEN 
            HK = DUMMY 
            RETURN 
         ELSEIF (NUSED.GT.NMAX) THEN 
            IERR = 5 
            RETURN 
         ENDIF 
!                                                                       
!** adapt step size for next step                                       
!                                                                       
         HALT = HK 
         HK = DMIN1 (Q, 2.0D0) * HK 
         IF (XK + HK.GE.XEND) THEN 
            DUMMY = HK 
            HK = XE-XK 
            IEND = .TRUE. 
!                                                                       
!** jump back if sufficiently close to XE                               
!                                                                       
            IF (HK.LT.EPS1 * DABS (XE) ) GOTO 275 
         ENDIF 
!                                                                       
!** Set up the Gera-Nordsieck approximation for the next                
!** integration                                                         
!                                                                       
         QUOT1 = HK / HALT 
         QUOT2 = QUOT1 * QUOT1 
         QUOT3 = QUOT2 * QUOT1 
         QUOT4 = QUOT3 * QUOT1 
         DO 280 I = 1, N 
            ZJ (0, I) = ZJP1 (0, I) 
            ZJ (1, I) = QUOT1 * ZJP1 (1, I) 
            ZJ (2, I) = QUOT2 * ZJP1 (2, I) 
            ZJ (3, I) = QUOT3 * ZJP1 (3, I) 
            ZJ (4, I) = QUOT4 * ZJP1 (4, I) 
  280    END DO 
      ELSE 
!                                                                       
!** Repeat last step for a smaller step size                            
!** and modify the Gear-Nordsieck approximation accordingly             
!                                                                       
         HALT = HK 
         HK = DMAX1 (0.5D0, Q) * HK 
         QUOT1 = HK / HALT 
         QUOT2 = QUOT1 * QUOT1 
         QUOT3 = QUOT2 * QUOT1 
         QUOT4 = QUOT3 * QUOT1 
         DO 290 I = 1, N 
            ZJ (1, I) = QUOT1 * ZJ (1, I) 
            ZJ (2, I) = QUOT2 * ZJ (2, I) 
            ZJ (3, I) = QUOT3 * ZJ (3, I) 
            ZJ (4, I) = QUOT4 * ZJ (4, I) 
  290    END DO 
         IEND = .FALSE. 
      ENDIF 
!                                                                       
      GOTO 75 
!                                                                       
      END SUBROUTINE GEAR4                          
