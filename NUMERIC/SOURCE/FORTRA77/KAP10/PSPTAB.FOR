C[BA*)
C[LE*)
      SUBROUTINE PSPTAB (N,NL,TBEG,TEND,T,AX,BX,CX,DX,
C[IX{PSPTAB}*)
     +                   AY,BY,CY,DY,NT,XTAB,YTAB,IERR)
C
C*****************************************************************
C                                                                *
C  Program to create a table of values for parametric cubic      *
C  splines                                                       *
C[BE*)
C  with component functions SX(T), SY(T), given in the           *
C  following form:                                               *
C                                                                *
C  SX := SX(T) = AX(I) + BX(I)(T-T(I)) + CX(I)(T-T(I))**2 +      *
C                                      + DX(I)(T-T(I))**3        *
C                                                                *
C  SY := SY(T) = AY(I) + BY(I)(T-T(I)) + CY(I)(T-T(I))**2 +      *
C                                      + DY(I)(T-T(I))**3        *
C                                                                *
C  for T in the interval [T(I),T(I+1)], I=0,1,...,N-1.           *
C                                                                *
C                                                                *
C  This program creates a table of function values XTAB = SX(TW) *
C  and YTAB = SY(TW) where TW lies in [TBEG,TEND]. We use the    *
C  following conventions:                                        *
C   - if TBEG < T(0), the end polynomial P(0) will be evaluated  *
C     for all values XTAB < T(0)                                 *
C   - if TEND > T(N), the end polynomial P(N-1) will be evaluated*
C     for all values XTAB > T(N)                                 *
C   - in every table the interval end points TBEG and TEND and   *
C     all nodes T(I) in between will occur in the table          *
C   - in each subinterval [T(I),T(I+1)] the table is created for *
C     equidistant steps of size H. Thus H will always depend on  *
C     the length of the given interval and on the length NL of   *
C     the table.                                                 *
C   - the input parameter NL presets an approximate table length;*
C     the actual table length is NT+1 (NT denotes the final index*
C     of the table of values). We must have  0 < NT < NL+N+3.    *
C                                                                *
C                                                                *
C  ASSUMPTIONS:   TBEG <  TEND                                   *
C  ============   NL   >= 0                                      *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  N    :  Index of the final knot T(N)                          *
C  NL   :  Given table length for dimensioning of vectors XTAB   *
C          and YTAB                                              *
C  TBEG :  starting value of the table                           *
C  TEND :  ending value of the table                             *
C  T    :  N+1-vector T(0:N); the nodes T(I), I=0,1,...,N        *
C  AX :  ] N+1-vectors ..(0:N);                                  *
C  BX :  ] the components in positions 0 to N-1 contain the      *
C  CX :  ] spline coefficients for the component function SX(T)  *
C  DX :  ]                                                       *
C                                                                *
C  AY :  ] N+1-vectors ..(0:N);                                  *
C  BY :  ] the components in positions 0 to N-1 contain the      *
C  CY :  ] spline coefficients for the component function SY(T)  *
C  DY :  ]                                                       *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  NT   :  final table index; equals actual table length - 1     *
C  XTAB :  vector XTAB(0:NL+N+2) ] the elements in positions     *
C  YTAB :  vector YTAB(0:NL+N+2) ] 0 to NT contain the table     *
C  IERR :  Error parameter                                       *
C          = 0 : Everything o.k.                                 *
C          = 1 : Stop caused by TBEG >= TEND                     *
C          = 2 : Stop caused by NL < 0                           *
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
      DOUBLE PRECISION XTAB(0:NL+N+2), YTAB(0:NL+N+2), T(0:N),
     +                 AX(0:N), BX(0:N), CX(0:N), DX(0:N),
     +                 AY(0:N), BY(0:N), CY(0:N), DY(0:N)
C
C-----checking the input------------------------------------------
C
      IF (TEND .LE. TBEG) THEN
        IERR = 1
        RETURN
      ELSEIF (NL .LT. 0) THEN
        IERR = 2
        RETURN
      ENDIF
      IERR = 0
C
C-----determine the interval [T(I),T(I+1)]-------------------------
C     that contains TBEG; index it by IBEG
C
      I = 0
      K = N
   10 M = (I+K)/2
      IF (TBEG .LT. T(M)) THEN
        K = M
      ELSE
        I = M
      ENDIF
      IF (K .GT. I+1) GOTO 10
      IBEG = I
C
C-----determine the interval [T(I),T(I+1)]------------------------
C     that contains TEND; index it by IEND
C
      K = N
   20 M = (I+K)/2
      IF (TEND .LT. T(M)) THEN
        K = M
      ELSE
        I = M
      ENDIF
      IF (K .GT. I+1) GOTO 20
      IEND = I
C
C-----determine the table values XTAB(I), YTAB(I), I=0,1,...,NT---
C
C     Initialize
C
      HP = TEND - TBEG
      FC = NL/HP
      NT = 0
      TW = TBEG
C
      IF (IBEG .NE. IEND) THEN
C
        IF (TBEG .LT. T(0)) THEN
          IP = 0
        ELSE
          IP = 1
        ENDIF
C
        IF (TEND .GT. T(N)) THEN
          IM = 0
        ELSE
          IM = 1
        ENDIF
C
C     determine the values for the table at TBEG to T(IBEG+IP)
C
        I = IBEG
        TD = TW - T(I)
        XTAB(NT) = ((DX(I)*TD + CX(I))*TD + BX(I))*TD + AX(I)
        YTAB(NT) = ((DY(I)*TD + CY(I))*TD + BY(I))*TD + AY(I)
        DIF = T(I+IP) - TBEG
        TIV = DIF*FC
        ITV = INT(TIV)
        IF ((TIV-ITV) .GT. 0.0D0) ITV = ITV+1
        IF (ITV .GT. 0) H = DIF/ITV
        DO 30 J = 1,ITV-1,1
          NT = NT + 1
          TW = TW + H
          TD = TW - T(I)
          XTAB(NT) = ((DX(I)*TD + CX(I))*TD + BX(I))*TD + AX(I)
          YTAB(NT) = ((DY(I)*TD + CY(I))*TD + BY(I))*TD + AY(I)
   30   CONTINUE
        NT = NT + 1
        IF ((IEND-IBEG) .NE. 1) THEN
C
C           determine the table values at T(IBEG+IP) to
C           T(IEND-IM+1)
C
          IBP = IBEG + IP
          IEM = IEND - IM
          DO 40 I = IBP,IEM,1
            TW = T(I)
            XTAB(NT) = AX(I)
            YTAB(NT) = AY(I)
            DIF = T(I+1) - T(I)
            TIV = DIF*FC
            ITV = INT(TIV)
            IF ((TIV-ITV) .GT. 0.0D0) ITV = ITV+1
            IF (ITV .GT. 0) H = DIF/ITV
            DO 50 J = 1,ITV-1,1
              NT = NT + 1
              TW = TW + H
              TD = TW - T(I)
              XTAB(NT) = ((DX(I)*TD+CX(I))*TD+BX(I))*TD+AX(I)
              YTAB(NT) = ((DY(I)*TD+CY(I))*TD+BY(I))*TD+AY(I)
   50       CONTINUE
            NT = NT + 1
   40     CONTINUE
        ENDIF
        TW = T(IEND-IM+1)
      ENDIF
C
C     determine the table values from the location
C     tabulated last to TEND
C
      I = IEND
      TD = TW - T(I)
      XTAB(NT) = ((DX(I)*TD + CX(I))*TD + BX(I))*TD + AX(I)
      YTAB(NT) = ((DY(I)*TD + CY(I))*TD + BY(I))*TD + AY(I)
      DIF = TEND - TW
      TIV = DIF*FC
      ITV = INT(TIV)
      IF ((TIV-ITV) .GT. 0.0D0) ITV = ITV+1
      IF (ITV .GT. 0) H = DIF/ITV
      DO 60 J = 1,ITV-1,1
        NT = NT + 1
        TW = TW + H
        TD = TW - T(I)
        XTAB(NT) = ((DX(I)*TD + CX(I))*TD + BX(I))*TD + AX(I)
        YTAB(NT) = ((DY(I)*TD + CY(I))*TD + BY(I))*TD + AY(I)
   60 CONTINUE
      NT = NT + 1
      TD = TEND - T(I)
      XTAB(NT) = ((DX(I)*TD + CX(I))*TD + BX(I))*TD + AX(I)
      YTAB(NT) = ((DY(I)*TD + CY(I))*TD + BY(I))*TD + AY(I)
      RETURN
      END
