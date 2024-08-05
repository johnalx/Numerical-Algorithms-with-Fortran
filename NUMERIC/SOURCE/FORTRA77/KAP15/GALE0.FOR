C[BA*)
C[LE*)
      SUBROUTINE GALE0 (IDGR, FLAG, ALPHA, ZWGH)
C[IX{GALE0}*)
C
C*****************************************************************
C                                                                *
C  The SUBROUTINE GALE0 determines the nodes and weights of the  *
C  GAUSS-quadrature formula of degree IDGR.                      *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  IDGR  : degree of the quadrature formula                      *
C  FLAG  : = .TRUE.  if the the weights shall be determined      *
C          = .FALSE. only the nodes are to be determined         *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS                                             *
C  ==================                                            *
C  ALPHA : vector ALPHA(1:21); contains the nodes of the Gauss   *
C          quadrature formula of degree IDGR                     *
C  ZWGH  : vector ZWGH(1:21); contains the weights for the nodes *
C          ALPHA                                                 *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required : GXPOLY, GXPEGA, MACHPD                 *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Hermann-Josef Rheinbach                            *
C  editor   : Norbert Vogt                                       *
C  date     : 04.10.1989                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION POLD(0:100), PMID(0:100), C(0:100)
      DIMENSION ALPHA(1:IDGR), ZWGH(1:IDGR)
      LOGICAL   FLAG
C
C     determine the coefficients of the LEGENDRE-polynomial
C
      POLD(0) = 1.0D0
      PMID(0) = 0.0D0
      PMID(1) = 1.0D0
      XK      = 0.0D0
      ZPDGR   = DBLE(IDGR)
      KPLUS   = IDGR
C
      DO 10 K = 1, IDGR-1
         XK     = XK + 1.0D0
         XKPLUS = XK + 1.0D0
         XKINV  = 1.0D0 / XKPLUS
         XFA    = (XK + XKPLUS) * XKINV
         DO 20 I = 0, K
            C(I+1) = PMID(I) * XFA
 20      CONTINUE
         C(0) = 0.0D0
         XFA  = XK * XKINV
         DO 30 I = 0, K-1
            C(I)    = C(I) - POLD(I) * XFA
            POLD(I) = PMID(I)
 30      CONTINUE
         POLD(K) = PMID(K)
         DO 40 I = 0, K+1
            PMID(I) = C(I)
 40      CONTINUE
C
   10 CONTINUE
C
C     determine the zeros if lying symmetrically to the
C     origin with corresponding weights
C
      BORDRA = 1.0D0
      ZW     = 3.141592654D0 / (ZPDGR - 0.5D0)
      KDIV2  = INT(ZPDGR*0.5D0)
      DO 50 J = 1, KDIV2
         ZJ     = DBLE(J)
         BORDRB = 0.5D0 * DBLE((DCOS((ZJ-0.5D0)*ZW) + DCOS(ZJ*ZW)))
         CALL GXPEGA (BORDRA, BORDRB, C, KPLUS, XSI)
         ALPHA(J) = XSI
         BORDRA   = BORDRB
 50   CONTINUE
      DO 55 I = 1, KDIV2
         ALPHA(I) = -ALPHA(I)
 55   CONTINUE
      IPOS = KDIV2
      IF ((MOD(IDGR,2) .EQ. 1)) THEN
         IPOS        = KDIV2 + 1
         ALPHA(IPOS) = 0.0D0
         DO 60 I = 1, KDIV2
            ALPHA(IPOS+I) = -ALPHA(IPOS-I)
 60      CONTINUE
      ELSE
         DO 70 I = 1, KDIV2
            ALPHA(IPOS+I) = -ALPHA(IPOS+1-I)
 70      CONTINUE
      ENDIF
C
C     determine the weights for the nodes
C
      IF (FLAG) THEN
         DO 80 I = 1, IDGR
            XNULL = ALPHA(I)
            CALL GXPOLY (XNULL, F, IDGR-1, POLD)
            XW      = XKPLUS*XKPLUS * F*F
            ZWGH(I) = 2.0D0 * ( 1.0D0 - XNULL*XNULL ) / XW
 80      CONTINUE
      ENDIF
      R E T U R N
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE GXPOLY (X, F, N, C)
C[IX{GXPOLY}*)
C
C*****************************************************************
C                                                                *
C  This subroutine evaluates a polynomial of degree N at X       *
C  using the HORNER scheme. The value appears in F.              *
C                                                                *
C----------------------------------------------------------------*
C[BE*)
C                                                                *
C  subroutines required : none                                   *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Hermann-Josef Rheinbach                            *
C  editor   : Norbert Vogt                                       *
C  date     : 04.10.1989                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION C(0:N)
C
      NDGR = N-1
      F    = C(N)
      DO 10 K = 0, NDGR
         F = F * X + C(NDGR-K)
 10   CONTINUE
      R E T U R N
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE GXPEGA (A, B, C, N, XSI)
C[IX{GXPEGA}*)
C
C*****************************************************************
C                                                                *
C  This subroutine determines a zero of a polynomial of degree N,*
C  whose coefficients are given in the vector C(N) which lies    *
C  between A and B.                                              *
C  The method used is a modified version of the PEGASUS method.  *
C                                                                *
C----------------------------------------------------------------*
C[BE*)
C                                                                *
C  subroutines required: GXPOLY, MACHPD                          *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Hermann-Josef Rheinbach                            *
C  editor   : Norbert Vogt                                       *
C  date     : 04.10.1989                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION  C(0:N)
      LOGICAL    FDELTA
C
C     determining the machine constant EPS, and using  EPS * 100
C
      SAVE DELTA, FDELTA
      DATA FDELTA /.TRUE./
C
C     FDELTA = TRUE
C
      IF (FDELTA) THEN
         DELTA = 1.0D0
C
C        Repeat until machine constant found
C
  5      DELTA = 0.5D0 * DELTA
         IF (MACHPD(1.0D0+DELTA) .EQ. 1) GOTO 5
C
         DELTA  = 200.0D0 * DELTA
         FDELTA = .FALSE.
      ENDIF
C
C     initialize
C
      X1 = A
      X2 = B
C
      CALL GXPOLY (X1, F1, N, C)
      CALL GXPOLY (X2, F2, N, C)
      XDIFF = X2-X1
C
      DO 10 I = 1, 50
         S12 = XDIFF / (F2-F1)
         X3  = X2 - F2 * S12
         CALL GXPOLY (X3, F3, N, C)
         IF ((F2*F3) .LT. 0.0D0) THEN
            X1 = X2
            F1 = F2
         ELSE
            F1 = F1 * F2 / (F2 + F3)
         ENDIF
         X2 = X3
         F2 = F3
         IF (DABS(F2) .LT. DELTA) THEN
            XSI = X2
            IF (DABS(F1) .LT. DABS(F2)) XSI = X1
            R E T U R N
         ENDIF
         XDIFF = X2 - X1
         IF (DABS(XDIFF) .LT. DELTA) THEN
            XSI = X2
            IF (DABS(F1) .LT. DABS(F2)) XSI = X1
            R E T U R N
         ENDIF
   10 CONTINUE
      R E T U R N
      END
