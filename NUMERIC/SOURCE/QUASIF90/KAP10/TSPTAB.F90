      SUBROUTINE TSPTAB (N, NL, PBEG, PEND, PHIN, A, B, C, D, PHIR, PX, &
      PY, NT, XTAB, YTAB, IERR)                                         
!                                                                       
!*****************************************************************      
!                                                                *      
!  This program creates a table of values for a transformed      *      
!  parametric cubic spline function, given in the form:          *      
!                                                                *      
!  S(PHI) = A(I) + B(I)(PHI-PHIN(I)) + C(I)(PHI-PHIN(I))**2 +    *      
!                                    + D(I)(PHI-PHIN(I))**3      *      
!                                                                *      
!  for PHI in the interval [PHIN(I),PHIN(I+1)], I=0,1,...,N-1.   *      
!                                                                *      
!  The program creates a table containing                        *      
!      XTAB := XTAB(PHI) = S(PHI)*COS(PHI+PHIR) + PX, and        *      
!      YTAB := YTAB(PHI) = S(PHI)*SIN(PHI+PHIR) + PY.            *      
!  Here PHI lies in [PBEG,PEND] and the following conventions    *      
!  are used:                                                     *      
!   - if PBEG < PHIN(0), the end polynomial P(0) is evaluated    *      
!     for all values XTAB < PHIN(0)                              *      
!   - if PEND > PHIN(N) the end polynomial P(N-1) is evaluated   *      
!     for all values XTAB > PHIN(N)                              *      
!   - in every table the interval end points PBEG and PEND and   *      
!     all nodes PHIN(I) in between will be used in the table     *      
!   - in each subinterval [PHIN(I),PHIN(I+1)] the table is       *      
!     created for equidistant steps of size H. Thus H will always*      
!     depend on the length of the given interval and on the      *      
!     length NL of the table.                                    *      
!   - the input parameter NL presets an approximate table length;*      
!     the actual table length is NT+1 (NT denotes the final index*      
!     of the table of values). We must have  0 < NT < NL+N+3.    *      
!                                                                *      
!                                                                *      
!  ASSUMPTIONS:   PBEG <  PEND                                   *      
!  ============   NL   >= 0                                      *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N    :  Index of the final node PHIN(N)                       *      
!  NL   :  Table length given for dimensioning of the vectors    *      
!          XTAB and YTAB                                         *      
!  XBEG :  Starting table value                                  *      
!  XEND :  Final table value                                     *      
!  PHIN :  vector PHIN(0:N); the nodes PHIN(I), I=0,1,..,N       *      
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
!  NT   :  Final index for the table; equal to the actual table  *      
!          length - 1                                            *      
!  XTAB :  vector XTAB(0:NL+N+2) ] The elements in positions 0   *      
!  YTAB :  vector YTAB(0:NL+N+2) ] to NT form the table of values*      
!  IERR :  Error parameter                                       *      
!          = 0 : Everything o.k.                                 *      
!          = 1 : Stop because PBEG >= PEND                       *      
!          = 2 : Stop because NL < 0                             *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Guenter Palm                                       *      
!  date     : 03.28.1989                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
!-----declarations------------------------------------------------      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DOUBLEPRECISION PHIN (0:N), A (0:N), B (0:N), C (0:N), D (0:N),   &
      XTAB (0:NL + N + 2), YTAB (0:NL + N + 2)                          
!                                                                       
!-----checking the input data-------------------------------------      
!                                                                       
      IF (PEND.LE.PBEG) THEN 
         IERR = 1 
         RETURN 
      ELSEIF (NL.LT.0) THEN 
         IERR = 2 
         RETURN 
      ENDIF 
      IERR = 0 
!                                                                       
!-----determine the interval [PHIN(I),PHIN(I+1)] ----------------       
!     which includes PBEG; label it IBEG                                
!                                                                       
      I = 0 
      K = N 
   10 M = (I + K) / 2 
      IF (PBEG.LT.PHIN (M) ) THEN 
         K = M 
      ELSE 
         I = M 
      ENDIF 
      IF (K.GT.I + 1) GOTO 10 
      IBEG = I 
!                                                                       
!-----determine the interval [PHIN(I),PHIN(I+1)] ----------             
!     which includes PEND; label it IEND                                
!                                                                       
      K = N 
   20 M = (I + K) / 2 
      IF (PEND.LT.PHIN (M) ) THEN 
         K = M 
      ELSE 
         I = M 
      ENDIF 
      IF (K.GT.I + 1) GOTO 20 
      IEND = I 
!                                                                       
!-----determine the values XTAB(I), YTAB(I), I=0,1,...,NT ------        
!                                                                       
!     initialize                                                        
!                                                                       
      HP = PEND-PBEG 
      FC = NL / HP 
      NT = 0 
      PW = PBEG 
!                                                                       
      IF (IBEG.NE.IEND) THEN 
!                                                                       
         IF (PBEG.LT.PHIN (0) ) THEN 
            IP = 0 
         ELSE 
            IP = 1 
         ENDIF 
!                                                                       
         IF (PEND.GT.PHIN (N) ) THEN 
            IM = 0 
         ELSE 
            IM = 1 
         ENDIF 
!                                                                       
!        determine the table values from PBEG to                        
!        PHIN(IBEG+IP)                                                  
!                                                                       
         I = IBEG 
         PD = PW - PHIN (I) 
         S = ( (D (I) * PD+C (I) ) * PD+B (I) ) * PD+A (I) 
         RHO = PW + PHIR 
         XTAB (NT) = S * DCOS (RHO) + PX 
         YTAB (NT) = S * DSIN (RHO) + PY 
         DIF = PHIN (IBEG + IP) - PBEG 
         TIV = DIF * FC 
         ITV = INT (TIV) 
         IF ( (TIV - ITV) .GT.0.0D0) ITV = ITV + 1 
         IF (ITV.GT.0) H = DIF / ITV 
         DO 30 J = 1, ITV - 1, 1 
            NT = NT + 1 
            PW = PW + H 
            PD = PW - PHIN (I) 
            S = ( (D (I) * PD+C (I) ) * PD+B (I) ) * PD+A (I) 
            RHO = PW + PHIR 
            XTAB (NT) = S * DCOS (RHO) + PX 
            YTAB (NT) = S * DSIN (RHO) + PY 
   30    END DO 
         NT = NT + 1 
         IF ( (IEND-IBEG) .NE.1) THEN 
!                                                                       
!           determine the table values from PHIN(IBEG+IP)               
!           to PHIN(IEND-IM+1)                                          
!                                                                       
            IBP = IBEG + IP 
            IEM = IEND-IM 
            DO 40 I = IBP, IEM, 1 
               PW = PHIN (I) 
               RHO = PW + PHIR 
               XTAB (NT) = A (I) * DCOS (RHO) + PX 
               YTAB (NT) = A (I) * DSIN (RHO) + PY 
               DIF = PHIN (I + 1) - PHIN (I) 
               TIV = DIF * FC 
               ITV = INT (TIV) 
               IF ( (TIV - ITV) .GT.0.0D0) ITV = ITV + 1 
               IF (ITV.GT.0) H = DIF / ITV 
               DO 50 J = 1, ITV - 1, 1 
                  NT = NT + 1 
                  PW = PW + H 
                  PD = PW - PHIN (I) 
                  S = ( (D (I) * PD+C (I) ) * PD+B (I) ) * PD+A (I) 
                  RHO = PW + PHIR 
                  XTAB (NT) = S * DCOS (RHO) + PX 
                  YTAB (NT) = S * DSIN (RHO) + PY 
   50          END DO 
               NT = NT + 1 
   40       END DO 
         ENDIF 
         PW = PHIN (IEND-IM + 1) 
      ENDIF 
!                                                                       
!     determine the table values from the location which                
!     was tabulated last to PEND                                        
!                                                                       
      PD = PW - PHIN (IEND) 
      S = ( (D (IEND) * PD+C (IEND) ) * PD+B (IEND) ) * PD+A (IEND) 
      RHO = PW + PHIR 
      XTAB (NT) = S * DCOS (RHO) + PX 
      YTAB (NT) = S * DSIN (RHO) + PY 
      DIF = PEND-PW 
      TIV = DIF * FC 
      ITV = INT (TIV) 
      IF ( (TIV - ITV) .GT.0.0D0) ITV = ITV + 1 
      IF (ITV.GT.0) H = DIF / ITV 
      DO 60 J = 1, ITV - 1, 1 
         NT = NT + 1 
         PW = PW + H 
         PD = PW - PHIN (IEND) 
         S = ( (D (IEND) * PD+C (IEND) ) * PD+B (IEND) ) * PD+A (IEND) 
         RHO = PW + PHIR 
         XTAB (NT) = S * DCOS (RHO) + PX 
         YTAB (NT) = S * DSIN (RHO) + PY 
   60 END DO 
      NT = NT + 1 
      PD = PEND-PHIN (IEND) 
      S = ( (D (IEND) * PD+C (IEND) ) * PD+B (IEND) ) * PD+A (IEND) 
      RHO = PEND+PHIR 
      XTAB (NT) = S * DCOS (RHO) + PX 
      YTAB (NT) = S * SIN (RHO) + PY 
      RETURN 
      END SUBROUTINE TSPTAB                         
