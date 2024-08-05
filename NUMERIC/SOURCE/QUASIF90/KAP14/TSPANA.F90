      SUBROUTINE TSPANA (PHI, N, PHIN, A, B, C, D, PHIR, PX, PY, S, S1, &
      S2, S3, XK, YK, C1, CCR, IERR1, IERR2)                            
!                                                                       
!*****************************************************************      
!                                                                *      
!  Evaluation program for transformed parametric cubic spline    *      
!  functions formatted as follows                                *      
!                                                                *      
!  S(PHI) = A(I) + B(I)(PHI-PHIN(I)) + C(I)(PHI-PHIN(I))**2 +    *      
!                                    + D(I)(PHI-PHIN(I))**3      *      
!                                                                *      
!  for PHI in the interval [PHIN(I),PHIN(I+1)], I=0,1,..,N-1.    *      
!                                                                *      
!  This program determines the functional value and that of the  *      
!  1st, 2nd and 3rd derivative of a spline function S(PHI), as   *      
!  well as the values XK=F(S(PHI)), YK=F(S(PHI)), and the 1st    *      
!  derivative and the curvature of the curve K at a given PHI.   *      
!                                                                *      
!                                                                *      
!  NOTE:  This evaluation program is not well suited for making  *      
!  =====  a table of values for S(PHI) or the curve determined   *      
!         by the points XK, YK.                                  *      
!         If one only wants to evaluate the spline function      *      
!         S(PHI), we recommend the subroutine SPVAL instead.     *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  PHI  :  Position (in radians) where we want to evaluate the   *      
!          spline                                                *      
!  N    :  Number of the final node PHIN(N)                      *      
!  PHIN :  vector PHIN(0:N); the nodes PHIN(I), I=0,1,...,N      *      
!  A    :  ] N+1-vectors ..(0:N);                                *      
!  B    :  ] the elements in positions 0 to N-1 describe the     *      
!  C    :  ] coefficients of the spline function S(PHI)          *      
!  D    :  ]                                                     *      
!                                                                *      
!  PHIR :  ] the rotation angle PHIR and the translation vector  *      
!  PX   :  ] (PX,PY) are outputs of the subroutine ISPLTR for    *      
!  PY   :  ] interpolating splines and of CFSPTR for fitting     *      
!          ] splines.                                            *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  S  :  Function value of the spline function at PHI            *      
!  S1 :  1st derivative "   "     "      "     "   "             *      
!  S2 :  2nd derivative "   "     "      "     "   "             *      
!  S3 :  3rd derivative "   "     "      "     "   "             *      
!  XK : ] Coordinates of curve K at PHI                          *      
!  YK : ]                                                        *      
!  C1 :  1st derivative of curve K at PHI                        *      
!        It is determined by the equation:                       *      
!        C1 = (S1*SIN(RHO)+S*COS(RHO))/(S1*COS(RHO)-S*SIN(RHO)), *      
!              with RHO = PHI + PHIR.                            *      
!  CCR:  Curvature of curve K at location PHI.                   *      
!        It is determined by the equation:                       *      
!        CCR = (2*S1**2 - S*S2 + S**2)/((S1**2 + S**2)**1.5).    *      
!                                                                *      
!  IERR1 :  Error parameter for determining C1                   *      
!           = 0 : Everything o.k.                                *      
!           = 1 : The denominator of the equation for C1 is zero;*      
!                 C1 was not determined                          *      
!           = 2 : The magnitude of the denominator in the        *      
!                 equation for C1 is not zero; but it is less    *      
!                 than four times the machine constant.          *      
!                 A value for C1 cannot be accurately determined.*      
!  IERR2 :  Error parameter for determining CCR                  *      
!           = 0 : Everything o.k.                                *      
!           = 1 : The denominator of the equation for CCR is zero*      
!                 CCR was not determined                         *      
!           = 2 : The magnitude of the denominator in the        *      
!                 equation for CCR is not zero, but it is less   *      
!                 than four times the machine constant.          *      
!                 A value for CCR could not be accurately        *      
!                 determined.                                    *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Subroutines required: SPLFVD, MACHPD                          *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Guenter Palm                                       *      
!  date     : 04.15.1988                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
!-----declarations------------------------------------------------      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION PHIN (0:N), A (0:N), B (0:N), C (0:N), D (0:N) 
      LOGICAL FLAG 
      SAVE FMACHP, FLAG 
!                                                                       
!-----initializing------------------------------------------------      
!                                                                       
      DATA FLAG / .TRUE. / 
      TWOPI = 8.0D0 * DATAN (1.0D0) 
      IERR1 = 0 
      IERR2 = 0 
!                                                                       
!-----determine the machine constant (only on 1st call)----------       
!                                                                       
      IF (FLAG) THEN 
         FMACHP = 1.0D0 
   10    FMACHP = 0.5D0 * FMACHP 
         IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 10 
         FMACHP = 8.0D0 * FMACHP 
         FLAG = .FALSE. 
      ENDIF 
!                                                                       
!-----assign PHI to the auxiliary variable PHIX, -----------------      
!     if necessary convert PHIX so that PHIX lies                       
!     in the interval [0,2*PI]                                          
!                                                                       
      IF (PHI.LT.0.0D0) THEN 
         L = ABS (INT (PHI / TWOPI) ) + 1 
         PHIX = L * TWOPI - PHI 
      ELSEIF (PHI.GT.TWOPI) THEN 
         L = INT (PHI / TWOPI) 
         PHIX = PHI - L * TWOPI 
      ELSE 
         PHIX = PHI 
      ENDIF 
!                                                                       
!-----determine the functional value S and that of the derivatives-     
!     at PHIX in SUBROUTINE SPLFVD                                      
!                                                                       
      CALL SPLFVD (PHIX, N, PHIN, A, B, C, D, S, S1, S2, S3) 
!                                                                       
!-----Determine the coordinates XK, YK of the curve, ------------       
!     as well as the 1st derivative and the curvature                   
!                                                                       
      RHO = PHIX + PHIR 
      COSA = DCOS (RHO) 
      SINA = DSIN (RHO) 
      XK = S * COSA + PX 
      YK = S * SINA + PY 
      HZ = S1 * SINA + S * COSA 
      HN = S1 * COSA - S * SINA 
      IF (HN.EQ.0.0D0) THEN 
         IERR1 = 1 
      ELSE 
         IF (DABS (HN) .LE.FMACHP) IERR1 = 2 
         C1 = HZ / HN 
      ENDIF 
      HZ = 2.0D0 * S1 * S1 - S * S2 + S * S 
      HN = (S1 * S1 + S * S) **1.5D0 
      IF (HN.EQ.0.0D0) THEN 
         IERR2 = 1 
      ELSE 
         IF (DABS (HN) .LE.FMACHP) IERR2 = 2 
         CCR = HZ / HN 
      ENDIF 
      RETURN 
      END SUBROUTINE TSPANA                         
