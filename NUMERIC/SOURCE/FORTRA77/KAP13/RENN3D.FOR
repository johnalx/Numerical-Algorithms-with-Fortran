C[BA*)
C[LE*)
      SUBROUTINE RENN3D (N,XN,YN,ZN,NK,BETA,B,C,D,T,HELP,IMARK,IERR)
C[IX{RENN3D}*)
C
C*****************************************************************
C                                                                *
C  The programm RENN3D computes the coefficient vectors B, C and *
C  D, as well as the lengths T(I) of the parameter intervals of  *
C  a closed or open parametric 3-dimensional Renner spline.      *
C[BE*)
C  The subsplines are represented in vector form:                *
C                                                                *
C  [X(T)]=[PIX(T)]=[XN(I)+B(I,1)*T+C(I,1)*T**2+D(I,1)*T**3]      *
C  [Y(T)]=[PIY(T)]=[YN(I)+B(I,2)*T+C(I,2)*T**2+D(I,2)*T**3]      *
C  [Z(T)]=[PIZ(T)]=[ZN(I)+B(I,3)*T+C(I,3)*T**2+D(I,3)*T**3]      *
C                                                                *
C  for I=0, ..., N-1 and T a point in the interval [0, T(I)].    *
C                                                                *
C                                                                *
C  ASSUMPTIONS:                                                  *
C  ============    1. N >= 4 or NK >= 4                          *
C                  2. If 0.0 < BETA < 1.0: we must have          *
C                     NK >= N + INT((N+1)/2); otherwise NK = N   *
C                  2. The node (XN(I), YN(I), ZN(I)) cannot be   *
C                     equal to  (XN(I+1),YN(I+1), ZN(I+1))       *
C                     for any I = 0, ..., N-1                    *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  N       : Index of the last node                              *
C  XN      : DOUBLE PRECISION (NK+1)-vector XN(0:NK), containing *
C            the x-coordinates XN(I) of the nodes for I=0,..., N *
C  YN      : DOUBLE PRECISION (NK+1)-vector YN(0:NK), containing *
C            the y-coordinates YN(I) of the nodes for I=0,..., N *
C  ZN      : DOUBLE PRECISION (NK+1)-vector ZN(0:NK), containing *
C            the y-coordinates ZN(I) of the nodes for I=0,..., N *
C  NK      : NK = N + INT((N+1)/2), the maximal number of nodes  *
C            when using rounded corners,                         *
C            i.e. for 0.0 < BETA < 1.0.                          *
C            Without rounding of corners we must have: NK = N.   *
C  BETA    : for 0.0 < BETA < 1.0 corners are rounded; otherwise *
C            corners are kept as corners                         *
C                                                                *
C                                                                *
C  AUXILIARY PARAMETERS:                                         *
C  =====================                                         *
C  HELP    : DOUBLE PRECISION array HELP(-2:NK+1,1:14)           *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  N       : Index of the last node. If 0.0 < BETA < 1.0 , then  *
C            the output value for N can differ from its input.   *
C            When corners are rounded, the node list can be en-  *
C            larged by at most INT((N+1)/2) points.              *
C  T       : DOUBLE PRECISION NK-vector T(0:NK-1),] nodes        *
C            containing the parameter values      ] and          *
C  XN      : DOUBLE PRECISION vector XN(0:NK)     ] coefficients *
C  YN      : DOUBLE PRECISION vector YN(0:NK)     ] of           *
C  ZN      : DOUBLE PRECISION vector ZN(0:NK)     ] the          *
C  B       : DOUBLE PRECISION array B(0:NK-1,1:3) ] subsplines   *
C  C       : DOUBLE PRECISION array C(0:NK-1,1:3) ] for          *
C  D       : DOUBLE PRECISION array D(0:NK-1,1:3) ] I=0, ...,N-1.*
C            If 0.0 < BETA < 1.0, then the nodes XN(I) and FN(I) *
C            can be different from their input values for I=0,...*
C            ..., N                                              *
C                                                                *
C  IMARK   : Pointer                                             *
C            = 0, Subspline is not closed                        *
C            = 1, Subspline is closed                            *
C                                                                *
C  IERR    : Error parameter                                     *
C            =  0, all is ok                                     *
C            = -1, N < 4  or  NK < 4                             *
C            = -2, the assumption                                *
C                  (XN(I), YN(I), ZN(I)) not equal to            *
C                  (XN(I+1), YN(I+1), ZN(I+1)) is violated for   *
C                  some I=0, ..., N-1                            *
C            = -3, NK < N+INT((N+1)/2) while 0.0 < BETA < 1.0    *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Required subroutines: MACHPD                                  *
C                                                                *
C*****************************************************************
C                                                                *
C  Author     : Gisela Engeln-Muellges                           *
C  Date       : 4.7.1993                                         *
C  Source code: FORTRAN 77                                       *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C     Declarations
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XN(0:NK), YN(0:NK), ZN(0:NK), B(0:NK-1, 1:3),
     +                 C(0:NK-1, 1:3), D(0:NK-1, 1:3), T(0:NK-1),
     +                 HELP(-2:NK+1, 1:14)
      double precision sklprd
C
C     Check input
C
      IERR = 0
      IF (N .LT. 4 .OR. NK .LT. 4) THEN
         IERR = -1
         RETURN
      ENDIF
C
C     If  0.0 < BETA < 1.0  corners are rounded
C
      IF (BETA .GT. 0.0D0             .AND.
     F    BETA .LT. 1.0D0             .AND.
     F    NK   .LT. N+INT(0.5D0*(N+1))     ) THEN
            IERR = -3
            RETURN
      ENDIF
C
C     Compute the machine constant
C
      FMACHP = 1.0D0
    5 FMACHP = 0.5D0 * FMACHP
      IF (MACHPD(1.0D0+FMACHP) .EQ. 1) GOTO 5
      FMACHP = 2.0D0 * FMACHP
      EPS = 4.0D0 * FMACHP
C
C     Calculate the chordal vectors; store their three coordinates
C     in the first three columns of the auxiliary array HELP,
C     the fourth column of HELP contains the length of these chordal vectors.
C
      DO 10, I = 0, N-1
         HELP(I, 1) = XN(I+1) - XN(I)
         HELP(I, 2) = YN(I+1) - YN(I)
         HELP(I, 3) = ZN(I+1) - ZN(I)
         HELP(I, 4) = DSQRT(HELP(I, 1) * HELP(I, 1) +
     F                      HELP(I, 2) * HELP(I, 2) +
     F                      HELP(I, 3) * HELP(I, 3)   )
C
C     we check that all chordal vectors are nonzero, thus
C     verifying that consecutive nodes are distinct.
C
         IF(HELP(I,4) .LE. EPS) THEN
            IERR=-2
            RETURN
         ELSE
C
C    Columns 5 to 7 of HELP contain the three components
C    of the the chordal unit vectors
C
            HELP(I,5)=HELP(I,1)/HELP(I,4)
            HELP(I,6)=HELP(I,2)/HELP(I,4)
            HELP(I,7)=HELP(I,3)/HELP(I,4)
         ENDIF
   10 CONTINUE
C
C     We compute the area of the parallelograms generated by
C     two consecutive chordal vectors and store the results in
C     column 8 of the auxiliary array HELP
C
      DO 20, I = 0, N-2
         sklprd = HELP(I,5)*HELP(I+1,5) +
     F            HELP(I,6)*HELP(I+1,6) +
     F            HELP(I,7)*HELP(I+1,7)
         HELP(I,8)=DSQRT(1.0D0 - sklprd * sklprd)
   20 CONTINUE
C
C     Determine whether the curve is closed or not
C
      IMARK = 0
      IF (DABS(XN(0) - XN(N)) .LE. EPS .AND.
     F    DABS(YN(0) - YN(N)) .LE. EPS .AND.
     F    DABS(ZN(0) - ZN(N)) .LE. EPS      ) THEN
         IMARK = 1
      ENDIF
C
      IF(BETA.GT.0.0D0 .AND. BETA.LT.1.0D0 )THEN
C
C     If 0.0 < BETA < 1.0, round existing corners
C
         IF (IMARK .EQ. 1) THEN
C
C     Curve is closed
C
            XL=HELP(N-2,8)+HELP(0,8)
            sklprd = HELP(N-1,5)*HELP(0,5) +
     F               HELP(N-1,6)*HELP(0,6) +
     F               HELP(N-1,7)*HELP(0,7)
            XR=DSQRT(1.0D0 - sklprd * sklprd)
C
C     Condition for a corner
C
            IF(XL.LE.EPS .AND. XR.GT.EPS) THEN
C
C     Relabel points
C
               DO 30, I = 0, N-1
                  XN(I) = XN(I+1)
                  YN(I) = YN(I+1)
                  ZN(I) = ZN(I+1)
   30          CONTINUE
               XN(N) = XN(0)
               YN(N) = YN(0)
               ZN(N) = ZN(0)
C
C     Relabel chordal data
C
               DO 40 K = 1, 7
                  HELP(N,K) = HELP(0,K)
   40          CONTINUE
               DO 50, I = 0, N-1
                  DO 60, K = 1, 7
                     HELP(I,K) = HELP(I+1,K)
   60             CONTINUE
   50          CONTINUE
C
C     Relabel areas
C
               DO 70, I = 0, N-3
                  HELP(I,8) = HELP(I+1,8)
   70          CONTINUE
               sklprd = HELP(N-2,5)*HELP(N-1,5) +
     F                  HELP(N-2,6)*HELP(N-1,6) +
     F                  HELP(N-2,7)*HELP(N-1,7)
               HELP(N-2,8) = DSQRT(1.0D0 - sklprd * sklprd)
            ENDIF
C
C     Prepare loop
C
            I = 1
            IMAX = N-1
            sklprd = HELP(N-1,5)*HELP(0,5) +
     F               HELP(N-1,6)*HELP(0,6) +
     F               HELP(N-1,7)*HELP(0,7)
            HELP(-1,8) = DSQRT( 1.0D0 - sklprd * sklprd)
            HELP(N-1,8) = HELP(-1,8)
         ELSE
C
C     Curve not closed
C
C     Prepare loop
C
            I = 2
            IMAX = N-2
         ENDIF

   75    XL = HELP(I-2,8)+HELP(I,8)
         XR = HELP(I-1,8)
C
C     Round existing corners
C
         IF(XL.LE.EPS .AND. XR.GT.EPS)THEN
C
C     Relabel points  I to N
C
            DO 80, J = N, I, -1
               XN(J+1) = XN(J)
               YN(J+1) = YN(J)
               ZN(J+1) = ZN(J)
   80       CONTINUE
C
C    Relabel chordal data for I to N-1
C
            DO 90, J = N-1, I, -1
               DO 100, K = 1, 7
                  HELP(J+1,K) = HELP(J,K)
  100          CONTINUE
   90       CONTINUE
C
C     Relabel areas I to IMAX
C
            DO 110, J = IMAX, I, -1
               HELP(J+1,8) = HELP(J,8)
  110       CONTINUE
C
C     Create two new nodes  I and I+1
C
            XL = HELP(I-1,4)
            XR = HELP(I+1,4)
            XB = BETA*DMIN1(XL,XR)
            XLAMDA = XB/XL
            XMUE   = XB/XR

            XN(I) = XN(I)-XLAMDA*HELP(I-1,1)
            YN(I) = YN(I)-XLAMDA*HELP(I-1,2)
            ZN(I) = ZN(I)-XLAMDA*HELP(I-1,3)
            XN(I+1) = XN(I+1)+XMUE*HELP(I+1,1)
            YN(I+1) = YN(I+1)+XMUE*HELP(I+1,2)
            ZN(I+1) = ZN(I+1)+XMUE*HELP(I+1,3)
C
C     Compute new chordal data at I-1, I und I+1 and store
C
            DO 120, J = I-1, I+1
               HELP(J, 1) = XN(J+1) - XN(J)
               HELP(J, 2) = YN(J+1) - YN(J)
               HELP(J, 3) = ZN(J+1) - ZN(J)
               HELP(J, 4) = DSQRT(HELP(J, 1) * HELP(J, 1) +
     F                            HELP(J, 2) * HELP(J, 2) +
     F                            HELP(J, 3) * HELP(J, 3)   )
               HELP(J, 5) = HELP(J,1)/HELP(J,4)
               HELP(J, 6) = HELP(J,2)/HELP(J,4)
               HELP(J, 7) = HELP(J,3)/HELP(J,4)
  120       CONTINUE
C
C     Compute new areas at  I-2, I-1, I and I+1, and store
C
            DO 130, J = I-2, I+1
               sklprd = HELP(J,5)*HELP(J+1,5) +
     F                  HELP(J,6)*HELP(J+1,6) +
     F                  HELP(J,7)*HELP(J+1,7)
               HELP(J,8) = DSQRT(1.0D0 - sklprd * sklprd)
  130       CONTINUE
C
C     Increase point count
C
            N = N + 1
            IMAX = IMAX + 1
         ENDIF
C
C     Set index for next point
C
         I = I + 1
         IF(I.LE.IMAX) GOTO 75
C
      ENDIF
C
      IF(IMARK.EQ.1)THEN
C
C     Curve is closed
C
C     Compute additional chordal data and store
C
         DO 140, I = N-2, N-1
            DO 150, K = 1, 7
               HELP(I-N,K) = HELP(I,K)
  150       CONTINUE
  140    CONTINUE
         DO 160, I = 0, 1
            DO 170, K = 1, 7
               HELP(I+N,K) = HELP(I,K)
  170       CONTINUE
  160    CONTINUE
      ELSE
C
C     Curve not closed
C
         DO 180, K = 1, 3
            HELP(-2,K) = 3.D0*HELP(0,K)-2.D0*HELP(1,K)
            HELP(-1,K) = 2.D0*HELP(0,K)-HELP(1,K)
  180    CONTINUE
C
C     Compute additional chordal data at -2, -1, N and N+1
C
         DO 190, I = -2, -1
            HELP(I,4) = DSQRT( HELP(I,1)*HELP(I,1)
     F                        +HELP(I,2)*HELP(I,2)
     F                        +HELP(I,3)*HELP(I,3))
            IF(HELP(I,4).GT.EPS)THEN
               HELP(I, 5) = HELP(I,1)/HELP(I,4)
               HELP(I, 6) = HELP(I,2)/HELP(I,4)
               HELP(I, 7) = HELP(I,3)/HELP(I,4)
            ELSE
               HELP(I, 5) = 0.0D0
               HELP(I, 6) = 0.0D0
               HELP(I, 7) = 0.0D0
            ENDIF
  190    CONTINUE
         DO 200, K = 1, 3
            HELP(N,K)   = 2.D0*HELP(N-1,K)-HELP(N-2,K)
            HELP(N+1,K) = 3.D0*HELP(N-1,K)-2.D0*HELP(N-2,K)
  200    CONTINUE
         DO 210, I = N, N+1
            HELP(I,4) = DSQRT( HELP(I,1)*HELP(I,1)
     F                        +HELP(I,2)*HELP(I,2)
     F                        +HELP(I,3)*HELP(I,3))
            IF(HELP(I,4).GT.EPS)THEN
               HELP(I, 5) = HELP(I,1)/HELP(I,4)
               HELP(I, 6) = HELP(I,2)/HELP(I,4)
               HELP(I, 7) = HELP(I,3)/HELP(I,4)
            ELSE
               HELP(I, 5) = 0.0D0
               HELP(I, 6) = 0.0D0
               HELP(I, 7) = 0.0D0
            ENDIF
  210    CONTINUE
      ENDIF
C
C     Compute new areas for  -2, -1, N-1 and N
C
      DO 220, I = -2, -1
         sklprd = HELP(I,5)*HELP(I+1,5) +
     F            HELP(I,6)*HELP(I+1,6) +
     F            HELP(I,7)*HELP(I+1,7)
         HELP(I,8) = DSQRT(1.0D0 - sklprd * sklprd)
  220 CONTINUE
      DO 230, I = N-1, N
         sklprd = HELP(I,5)*HELP(I+1,5) +
     F            HELP(I,6)*HELP(I+1,6) +
     F            HELP(I,7)*HELP(I+1,7)
         HELP(I,8) = DSQRT(1.0D0 - sklprd * sklprd)
  230 CONTINUE
C
C     Compute the left and right hand unit tangent vectors and store in
C     columns 9 to 11 and 12 to 14 of HELP, respectively
C
      DO 240, I = 0, N
         XL = HELP(I-2,8)
         XR = HELP(I,8)
         IF(XL+XR.GT.EPS)THEN
            ALPHA = XL/(XL+XR)
            DO 250, K = 1,3
               HELP(I,K+8) = HELP(I-1,K)+ALPHA*(HELP(I,K)-HELP(I-1,K))
  250       CONTINUE
            HELPT = DSQRT(  HELP(I,9)*HELP(I,9)
     F                    + HELP(I,10)*HELP(I,10)
     F                    + HELP(I,11)*HELP(I,11) )
            HELP(I,9) = HELP(I,9)/HELPT
            HELP(I,10) = HELP(I,10)/HELPT
            HELP(I,11) = HELP(I,11)/HELPT
            HELP(I,12) = HELP(I,9)
            HELP(I,13) = HELP(I,10)
            HELP(I,14) = HELP(I,11)
         ELSE
            DO 260, K = 5, 7
               HELP(I,K+4) = HELP(I-1,K)
               HELP(I,K+7) = HELP(I,K)
  260       CONTINUE
         ENDIF
  240 CONTINUE
C
C     Compute the lengths of the parameter intervals
C
      DO 270, I = 0, N-1
         TS1 = HELP(I,12)+HELP(I+1,9)
         TS2 = HELP(I,13)+HELP(I+1,10)
         TS3 = HELP(I,14)+HELP(I+1,11)
         A1 = 16.0D0-(TS1*TS1+TS2*TS2+TS3*TS3)
         B1 = 6.0D0*(HELP(I,1)*TS1+HELP(I,2)*TS2+HELP(I,3)*TS3)
         C1 = 36.0D0*HELP(I,4)*HELP(I,4)
         T(I) = (-B1+DSQRT(B1*B1+A1*C1))/A1
         HELPT = 1.0D0/T(I)
C
C     Compute the spline coefficient vectors
C
         DO 280, K = 1, 3
            B(I,K) = HELP(I,K+11)
            C(I,K) = (3.0D0*HELPT*HELP(I,K)-2.0D0*B(I,K)
     F                                     -HELP(I+1,K+8))*HELPT
            D(I,K) = (B(I,K)+HELP(I+1,K+8)-2.0D0*HELPT*HELP(I,K))
     F                                              *HELPT*HELPT
  280    CONTINUE
  270 CONTINUE
C
      RETURN
      END
