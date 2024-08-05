C[BA*)
C[LE*)
      SUBROUTINE SPTAB (N,NL,XBEG,XEND,XN,A,B,C,D,NT,XTAB,YTAB,
C[IX{SPTAB}*)
     +                  IERR)
C
C*****************************************************************
C                                                                *
C  This program produces a table of values for non-parametric    *
C  cubic spline functions                                        *
C[BE*)
C  of the form:                                                  *
C                                                                *
C  S := S(X) = A(I) + B(I)(X-XN(I)) + C(I)(X-XN(I))**2 +         *
C                                   + D(I)(X-XN(I))**3           *
C                                                                *
C  for X in the interval [XN(I),XN(I+1)] for I=0,1,...,N-1.      *
C                                                                *
C  This program creates a table of values XTAB and YTAB=S(XTAB), *
C  where XTAB lies in [XBEG,XEND]. The following conventions are *
C  used:                                                         *
C   - if XBEG < XN(0), the end polynomial P(0) is evaluated for  *
C     all values XTAB < XN(0)                                    *
C   - if XEND > XN(N), the end polynomial P(N-1) is evaluated for*
C     all values XTAB > XN(N)                                    *
C   - in every table the interval end points XBEG and XEND and   *
C     all nodes XN(I) in between will be used in the table       *
C   - in each subinterval [XN(I),XN(I+1)] the table is created   *
C     for equidistant steps of size H. Thus H will always depend *
C     on the length of the given interval and on the length NL of*
C     the table.                                                 *
C   - the input parameter NL presets an approximate table length;*
C     the actual table length is NT+1 (NT denotes the final index*
C     of the table of values). We must have  0 < NT < NL+N+3.    *
C                                                                *
C                                                                *
C  ASSUMPTIONS:   XBEG <  XEND                                   *
C  ============   NL   >= 0                                      *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  N    :  Index of the final node XN(N)                         *
C  NL   :  Table length given for dimensioning of the vectors    *
C          XTAB and YTAB                                         *
C  XBEG :  Starting table value                                  *
C  XEND :  Final table value                                     *
C  XN   :  vector XN(0:N); the nodes XN(I), I=0,1,..,N           *
C  A    :  ] N+1-vectors A(0:N);                                 *
C  B    :  ] the elements in positions 0 to N-1 describe the     *
C  C    :  ] coefficients of the spline function S(X)            *
C  D    :  ]                                                     *
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
C          = 1 : Stop because XBEG >= XEND                       *
C          = 2 : Stop because NL < 0                             *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Subroutins required: none                                     *
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
      DOUBLE PRECISION XN(0:N), A(0:N), B(0:N), C(0:N), D(0:N),
     +                 XTAB(0:NL+N+2), YTAB(0:NL+N+2)
C
C-----checking the input data-------------------------------------
C
      IF (XEND .LE. XBEG) THEN
        IERR = 1
        RETURN
      ELSEIF (NL .LT. 0) THEN
        IERR = 2
        RETURN
      ENDIF
      IERR = 0
C
C-----determine the interval [XN(I),XN(I+1)] ---------------------
C     that contains XBEG; label it with IBEG
C
      I = 0
      K = N
   10 M = (I+K)/2
      IF (XBEG .LT. XN(M)) THEN
        K = M
      ELSE
        I = M
      ENDIF
      IF (K .GT. I+1) GOTO 10
      IBEG = I
C
C-----determine the interval [XN(I),XN(I+1)] ---------------------
C     that contains XEND; label it with IEND
C
      K = N
   20 M = (I+K)/2
      IF (XEND .LT. XN(M)) THEN
        K = M
      ELSE
        I = M
      ENDIF
      IF (K .GT. I+1) GOTO 20
      IEND = I
C
C-----determine the table values XTAB(I), YTAB(I), I=0,1,...,NT---
C
C     initialize
C
      HP = XEND - XBEG
      FC = NL/HP
      NT = 0
      XTAB(NT) = XBEG
C
      IF (IBEG .NE. IEND) THEN
C
        IF (XBEG .LT. XN(0)) THEN
          IP = 0
        ELSE
          IP = 1
        ENDIF
C
        IF (XEND .GT. XN(N)) THEN
          IM = 0
        ELSE
          IM = 1
        ENDIF
C
C        determine the table values for XBEG to
C        XN(IBEG+IP)
C
        XD = XTAB(NT) - XN(IBEG)
        YTAB(NT) = ((D(IBEG)*XD+C(IBEG))*XD+B(IBEG))*XD+A(IBEG)
        DIF = XN(IBEG+IP) - XBEG
        TIV = DIF*FC
        ITV = INT(TIV)
        IF ((TIV-ITV) .GT. 0.0D0) ITV = ITV+1
        IF (ITV .GT. 0) H = DIF/ITV
        DO 30 J = 1,ITV-1,1
          NT = NT + 1
          XTAB(NT) = XTAB(NT-1) + H
          XD = XTAB(NT) - XN(IBEG)
          YTAB(NT) = ((D(IBEG)*XD + C(IBEG))*XD +
     +                              B(IBEG))*XD + A(IBEG)
   30   CONTINUE
        NT = NT + 1
        IF ((IEND-IBEG) .NE. 1) THEN
C
C           determine the table values for XN(IBEG+IP)
C           to XN(IEND-IM+1)
C
          IBP = IBEG + IP
          IEM = IEND - IM
          DO 40 I = IBP,IEM,1
            XTAB(NT) = XN(I)
            YTAB(NT) = A(I)
            DIF = XN(I+1) - XN(I)
            TIV = DIF*FC
            ITV = INT(TIV)
            IF ((TIV-ITV) .GT. 0.0D0) ITV = ITV+1
            IF (ITV. GT. 0) H = DIF/ITV
            DO 50 J = 1,ITV-1,1
              NT = NT + 1
              XTAB(NT) = XTAB(NT-1) + H
              XD = XTAB(NT) - XN(I)
              YTAB(NT) = ((D(I)*XD + C(I))*XD + B(I))*XD + A(I)
   50       CONTINUE
            NT = NT + 1
   40     CONTINUE
        ENDIF
        XTAB(NT) = XN(IEND-IM+1)
      ENDIF
C
C     determine the table values from the last tabulated
C     location to XEND
C
      XD = XTAB(NT) - XN(IEND)
      YTAB(NT) = ((D(IEND)*XD+C(IEND))*XD+ B(IEND))*XD+A(IEND)
      DIF = XEND - XTAB(NT)
      TIV = DIF*FC
      ITV = INT(TIV)
      IF ((TIV-ITV) .GT. 0.0D0) ITV = ITV+1
      IF (ITV .GT. 0) H = DIF/ITV
      DO 60 J = 1,ITV-1,1
        NT = NT + 1
        XTAB(NT) = XTAB(NT-1) + H
        XD = XTAB(NT) - XN(IEND)
        YTAB(NT) = ((D(IEND)*XD + C(IEND))*XD +
     +                            B(IEND))*XD + A(IEND)
   60 CONTINUE
      NT = NT + 1
      XD = XEND - XN(IEND)
      XTAB(NT) = XEND
      YTAB(NT) = ((D(IEND)*XD+C(IEND))*XD+B(IEND))*XD+A(IEND)
      RETURN
      END
