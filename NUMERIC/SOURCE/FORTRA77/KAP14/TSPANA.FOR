C[BA*)
C[LE*)
      SUBROUTINE TSPANA (PHI,N,PHIN,A,B,C,D,PHIR,PX,PY,S,S1,S2,S3,
C[IX{TSPANA}*)
     +                   XK,YK,C1,CCR,IERR1,IERR2)
C
C*****************************************************************
C                                                                *
C  Evaluation program for transformed parametric cubic spline    *
C  functions                                                     *
C[BE*)
C  formatted as follows                                          *
C                                                                *
C  S(PHI) = A(I) + B(I)(PHI-PHIN(I)) + C(I)(PHI-PHIN(I))**2 +    *
C                                    + D(I)(PHI-PHIN(I))**3      *
C                                                                *
C  for PHI in the interval [PHIN(I),PHIN(I+1)], I=0,1,..,N-1.    *
C                                                                *
C  This program determines the functional value and that of the  *
C  1st, 2nd and 3rd derivative of a spline function S(PHI), as   *
C  well as the values XK=F(S(PHI)), YK=F(S(PHI)), and the 1st    *
C  derivative and the curvature of the curve K at a given PHI.   *
C                                                                *
C                                                                *
C  NOTE:  This evaluation program is not well suited for making  *
C  =====  a table of values for S(PHI) or the curve determined   *
C         by the points XK, YK.                                  *
C         If one only wants to evaluate the spline function      *
C         S(PHI), we recommend the subroutine SPVAL instead.     *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  PHI  :  Position (in radians) where we want to evaluate the   *
C          spline                                                *
C  N    :  Number of the final node PHIN(N)                      *
C  PHIN :  vector PHIN(0:N); the nodes PHIN(I), I=0,1,...,N      *
C  A    :  ] N+1-vectors ..(0:N);                                *
C  B    :  ] the elements in positions 0 to N-1 describe the     *
C  C    :  ] coefficients of the spline function S(PHI)          *
C  D    :  ]                                                     *
C                                                                *
C  PHIR :  ] the rotation angle PHIR and the translation vector  *
C  PX   :  ] (PX,PY) are outputs of the subroutine ISPLTR for    *
C  PY   :  ] interpolating splines and of CFSPTR for fitting     *
C          ] splines.                                            *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  S  :  Function value of the spline function at PHI            *
C  S1 :  1st derivative "   "     "      "     "   "             *
C  S2 :  2nd derivative "   "     "      "     "   "             *
C  S3 :  3rd derivative "   "     "      "     "   "             *
C  XK : ] Coordinates of curve K at PHI                          *
C  YK : ]                                                        *
C  C1 :  1st derivative of curve K at PHI                        *
C        It is determined by the equation:                       *
C        C1 = (S1*SIN(RHO)+S*COS(RHO))/(S1*COS(RHO)-S*SIN(RHO)), *
C              with RHO = PHI + PHIR.                            *
C  CCR:  Curvature of curve K at location PHI.                   *
C        It is determined by the equation:                       *
C        CCR = (2*S1**2 - S*S2 + S**2)/((S1**2 + S**2)**1.5).    *
C                                                                *
C  IERR1 :  Error parameter for determining C1                   *
C           = 0 : Everything o.k.                                *
C           = 1 : The denominator of the equation for C1 is zero;*
C                 C1 was not determined                          *
C           = 2 : The magnitude of the denominator in the        *
C                 equation for C1 is not zero; but it is less    *
C                 than four times the machine constant.          *
C                 A value for C1 cannot be accurately determined.*
C  IERR2 :  Error parameter for determining CCR                  *
C           = 0 : Everything o.k.                                *
C           = 1 : The denominator of the equation for CCR is zero*
C                 CCR was not determined                         *
C           = 2 : The magnitude of the denominator in the        *
C                 equation for CCR is not zero, but it is less   *
C                 than four times the machine constant.          *
C                 A value for CCR could not be accurately        *
C                 determined.                                    *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Subroutines required: SPLFVD, MACHPD                          *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Guenter Palm                                       *
C  date     : 04.15.1988                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C-----declarations------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION PHIN(0:N), A(0:N), B(0:N), C(0:N), D(0:N)
      LOGICAL FLAG
      SAVE FMACHP,FLAG
C
C-----initializing------------------------------------------------
C
      DATA FLAG /.TRUE./
      TWOPI = 8.0D0*DATAN(1.0D0)
      IERR1 = 0
      IERR2 = 0
C
C-----determine the machine constant (only on 1st call)----------
C
      IF (FLAG) THEN
        FMACHP = 1.0D0
   10   FMACHP = 0.5D0*FMACHP
        IF (MACHPD(1.0D0+FMACHP) .EQ. 1) GOTO 10
        FMACHP = 8.0D0*FMACHP
        FLAG = .FALSE.
      ENDIF
C
C-----assign PHI to the auxiliary variable PHIX, -----------------
C     if necessary convert PHIX so that PHIX lies
C     in the interval [0,2*PI]
C
      IF (PHI .LT. 0.0D0) THEN
        L = ABS(INT(PHI/TWOPI)) + 1
        PHIX = L*TWOPI - PHI
      ELSEIF (PHI .GT. TWOPI) THEN
        L = INT(PHI/TWOPI)
        PHIX = PHI - L*TWOPI
      ELSE
        PHIX = PHI
      ENDIF
C
C-----determine the functional value S and that of the derivatives-
C     at PHIX in SUBROUTINE SPLFVD
C
      CALL SPLFVD (PHIX,N,PHIN,A,B,C,D,S,S1,S2,S3)
C
C-----Determine the coordinates XK, YK of the curve, ------------
C     as well as the 1st derivative and the curvature
C
      RHO  = PHIX + PHIR
      COSA = DCOS(RHO)
      SINA = DSIN(RHO)
      XK   = S*COSA + PX
      YK   = S*SINA + PY
      HZ   = S1*SINA + S*COSA
      HN   = S1*COSA - S*SINA
      IF (HN .EQ. 0.0D0) THEN
        IERR1 = 1
      ELSE
        IF (DABS(HN) .LE. FMACHP) IERR1 = 2
        C1 = HZ/HN
      ENDIF
      HZ = 2.0D0*S1*S1 - S*S2 + S*S
      HN = (S1*S1 + S*S)**1.5D0
      IF (HN .EQ. 0.0D0) THEN
        IERR2 = 1
      ELSE
        IF (DABS(HN). LE. FMACHP) IERR2 = 2
        CCR = HZ/HN
      ENDIF
      RETURN
      END
