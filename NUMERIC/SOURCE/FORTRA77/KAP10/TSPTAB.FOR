C[BA*)
C[LE*)
      SUBROUTINE TSPTAB (N,NL,PBEG,PEND,PHIN,A,B,C,D,PHIR,PX,PY,
C[IX{TSPTAB}*)
     +                   NT,XTAB,YTAB,IERR)
C
C*****************************************************************
C                                                                *
C  This program creates a table of values for a transformed      *
C  parametric cubic spline function                              *
C[BE*)
C  given in the form:                                            *
C                                                                *
C  S(PHI) = A(I) + B(I)(PHI-PHIN(I)) + C(I)(PHI-PHIN(I))**2 +    *
C                                    + D(I)(PHI-PHIN(I))**3      *
C                                                                *
C  for PHI in the interval [PHIN(I),PHIN(I+1)], I=0,1,...,N-1.   *
C                                                                *
C  The program creates a table containing                        *
C      XTAB := XTAB(PHI) = S(PHI)*COS(PHI+PHIR) + PX, and        *
C      YTAB := YTAB(PHI) = S(PHI)*SIN(PHI+PHIR) + PY.            *
C  Here PHI lies in [PBEG,PEND] and the following conventions    *
C  are used:                                                     *
C   - if PBEG < PHIN(0), the end polynomial P(0) is evaluated    *
C     for all values XTAB < PHIN(0)                              *
C   - if PEND > PHIN(N) the end polynomial P(N-1) is evaluated   *
C     for all values XTAB > PHIN(N)                              *
C   - in every table the interval end points PBEG and PEND and   *
C     all nodes PHIN(I) in between will be used in the table     *
C   - in each subinterval [PHIN(I),PHIN(I+1)] the table is       *
C     created for equidistant steps of size H. Thus H will always*
C     depend on the length of the given interval and on the      *
C     length NL of the table.                                    *
C   - the input parameter NL presets an approximate table length;*
C     the actual table length is NT+1 (NT denotes the final index*
C     of the table of values). We must have  0 < NT < NL+N+3.    *
C                                                                *
C                                                                *
C  ASSUMPTIONS:   PBEG <  PEND                                   *
C  ============   NL   >= 0                                      *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  N    :  Index of the final node PHIN(N)                       *
C  NL   :  Table length given for dimensioning of the vectors    *
C          XTAB and YTAB                                         *
C  XBEG :  Starting table value                                  *
C  XEND :  Final table value                                     *
C  PHIN :  vector PHIN(0:N); the nodes PHIN(I), I=0,1,..,N       *
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
C  NT   :  Final index for the table; equal to the actual table  *
C          length - 1                                            *
C  XTAB :  vector XTAB(0:NL+N+2) ] The elements in positions 0   *
C  YTAB :  vector YTAB(0:NL+N+2) ] to NT form the table of values*
C  IERR :  Error parameter                                       *
C          = 0 : Everything o.k.                                 *
C          = 1 : Stop because PBEG >= PEND                       *
C          = 2 : Stop because NL < 0                             *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Guenter Palm                                       *
C  date     : 03.28.1989                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C-----declarations------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION PHIN(0:N), A(0:N), B(0:N), C(0:N), D(0:N),
     +                 XTAB(0:NL+N+2), YTAB(0:NL+N+2)
C
C-----checking the input data-------------------------------------
C
      IF (PEND .LE. PBEG) THEN
        IERR = 1
        RETURN
      ELSEIF (NL .LT. 0) THEN
        IERR = 2
        RETURN
      ENDIF
      IERR = 0
C
C-----determine the interval [PHIN(I),PHIN(I+1)] ----------------
C     which includes PBEG; label it IBEG
C
      I = 0
      K = N
   10 M = (I+K)/2
      IF (PBEG .LT. PHIN(M)) THEN
        K = M
      ELSE
        I = M
      ENDIF
      IF (K .GT. I+1) GOTO 10
      IBEG = I
C
C-----determine the interval [PHIN(I),PHIN(I+1)] ----------
C     which includes PEND; label it IEND
C
      K = N
   20 M = (I+K)/2
      IF (PEND .LT. PHIN(M)) THEN
        K = M
      ELSE
        I = M
      ENDIF
      IF (K .GT .I+1) GOTO 20
      IEND = I
C
C-----determine the values XTAB(I), YTAB(I), I=0,1,...,NT ------
C
C     initialize
C
      HP = PEND - PBEG
      FC = NL/HP
      NT = 0
      PW = PBEG
C
      IF (IBEG .NE. IEND) THEN
C
        IF (PBEG .LT. PHIN(0)) THEN
          IP = 0
        ELSE
          IP = 1
        ENDIF
C
        IF (PEND .GT. PHIN(N)) THEN
          IM = 0
        ELSE
          IM = 1
        ENDIF
C
C        determine the table values from PBEG to
C        PHIN(IBEG+IP)
C
        I = IBEG
        PD = PW - PHIN(I)
        S = ((D(I)*PD + C(I))*PD + B(I))*PD + A(I)
        RHO = PW + PHIR
        XTAB(NT) = S*DCOS(RHO) + PX
        YTAB(NT) = S*DSIN(RHO) + PY
        DIF = PHIN(IBEG+IP) - PBEG
        TIV = DIF*FC
        ITV = INT(TIV)
        IF ((TIV-ITV) .GT. 0.0D0) ITV = ITV+1
        IF (ITV .GT. 0) H = DIF/ITV
        DO 30 J = 1,ITV-1,1
          NT = NT + 1
          PW = PW + H
          PD = PW - PHIN(I)
          S = ((D(I)*PD + C(I))*PD + B(I))*PD + A(I)
          RHO = PW + PHIR
          XTAB(NT) = S*DCOS(RHO) + PX
          YTAB(NT) = S*DSIN(RHO) + PY
   30   CONTINUE
        NT = NT + 1
        IF ((IEND-IBEG) .NE. 1) THEN
C
C           determine the table values from PHIN(IBEG+IP)
C           to PHIN(IEND-IM+1)
C
          IBP = IBEG + IP
          IEM = IEND - IM
          DO 40 I = IBP,IEM,1
            PW = PHIN(I)
            RHO = PW + PHIR
            XTAB(NT) = A(I)*DCOS(RHO) + PX
            YTAB(NT) = A(I)*DSIN(RHO) + PY
            DIF = PHIN(I+1) - PHIN(I)
            TIV = DIF*FC
            ITV = INT(TIV)
            IF ((TIV-ITV) .GT. 0.0D0) ITV = ITV+1
            IF (ITV .GT. 0) H = DIF/ITV
            DO 50 J = 1,ITV-1,1
              NT = NT + 1
              PW = PW + H
              PD = PW - PHIN(I)
              S = ((D(I)*PD + C(I))*PD + B(I))*PD + A(I)
              RHO = PW + PHIR
              XTAB(NT) = S*DCOS(RHO) + PX
              YTAB(NT) = S*DSIN(RHO) + PY
   50       CONTINUE
            NT = NT + 1
   40     CONTINUE
        ENDIF
        PW = PHIN(IEND-IM+1)
      ENDIF
C
C     determine the table values from the location which
C     was tabulated last to PEND
C
      PD = PW - PHIN(IEND)
      S = ((D(IEND)*PD + C(IEND))*PD + B(IEND))*PD + A(IEND)
      RHO = PW + PHIR
      XTAB(NT) = S*DCOS(RHO) + PX
      YTAB(NT) = S*DSIN(RHO) + PY
      DIF = PEND - PW
      TIV = DIF*FC
      ITV = INT(TIV)
      IF ((TIV-ITV) .GT. 0.0D0) ITV = ITV+1
      IF (ITV .GT. 0) H = DIF/ITV
      DO 60 J = 1,ITV-1,1
        NT = NT + 1
        PW = PW + H
        PD = PW - PHIN(IEND)
        S = ((D(IEND)*PD + C(IEND))*PD + B(IEND))*PD + A(IEND)
        RHO = PW + PHIR
        XTAB(NT) = S*DCOS(RHO) + PX
        YTAB(NT) = S*DSIN(RHO) + PY
   60 CONTINUE
      NT = NT + 1
      PD = PEND - PHIN(IEND)
      S = ((D(IEND)*PD + C(IEND))*PD + B(IEND))*PD + A(IEND)
      RHO = PEND + PHIR
      XTAB(NT) = S*DCOS(RHO) + PX
      YTAB(NT) = S*SIN(RHO) + PY
      RETURN
      END
