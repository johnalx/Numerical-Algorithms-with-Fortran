      PROGRAM TESFKN
C
C*****************************************************************
C                                                                *
C  This program tests the subroutine  FFAKON  with values of the *
C  functions                                                     *
C       SIN(X)       for  X in [0, PI]                           *
C        0           otherwise                                   *
C  and                                                           *
C        0.5         for  X = 0                                  *
C       (PI-X)/PI    for  X in (0, PI]                           *
C        0           otherwise,                                  *
C  whose (nonperiodic) convolution is given by                   *
C       1 - X/PI - COS(X) + SIN(X)/PI   for  X  in [0, PI]       *
C       2 - X/PI + SIN(X)/PI            for  X  in (PI, 2*PI]    *
C       0                               otherwise,               *
C  while their (nonperiodic) correlation is                      *
C       1 + X/PI + SIN(X)/PI            for  X  in [-PI, 0]      *
C       X/PI + COS(X) + SIN(X)/PI       for  X  in (0, PI]       *
C       0                               otherwise.               *
C  We use equidistant nodes. An increase in the number of nodes  *  
C  must increase the accuracy of the approximation.              *
C                                                                *
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER M(8), INDEX, L(8), IFEHL
      DIMENSION F1RE(0:1023),F1IM(0:1023),F2RE(0:1023),F2IM(0:1023)
      DATA M / 10, 10, 41, 41, 128, 128, 500, 500 /
      DATA L / 32, 32, 128, 128, 512, 512, 1024, 1024 /
      INDEX = 1
C
C     Assign function values
C
      PI = 4.0D0*ATAN(1.0D0)
      DO 50 I = 1, 8
        DX = PI/DBLE(M(I))
        DO 10 J = 0, M(I)
          XJ = DBLE(J)*DX
          F1RE(J) = SIN(XJ)
          F1IM(J) = 0.0D0
          F2RE(J) = (PI-XJ)/PI
          F2IM(J) = 0.0D0
  10    CONTINUE
        F2RE(0) = 0.5D0
C
C     Compute convolution (INDEX=0) or correlation (INDEX=1)
C
      INDEX = 1 - INDEX
      CALL FFAKON(M(I),F1RE,F1IM,M(I),F2RE,F2IM,L(I),DX,INDEX,IFEHL)
C
C     Compare with theoretical output
C
        WRITE(4,*)'Error parameter :',IFEHL
        FFALT  = 0.0D0
        FFALTI = 0.0D0
        FKORR  = 0.0D0
        FKORRI = 0.0D0
        IF ( INDEX .EQ. 0) THEN
          DO 20 J = 0, 2*M(I)
            XJ    = DBLE(J)*DX
            IF(J .LE. M(I)) THEN
              FUFAL = 1.0D0-XJ/PI-COS(XJ)+SIN(XJ)/PI
            ELSE
              FUFAL = 2.0D0-XJ/PI+SIN(XJ)/PI
            ENDIF
            FFALT = MAX(FFALT, ABS(F1RE(J)-FUFAL))
            FFALTI= MAX(FFALTI, ABS(F1IM(J)))
  20      CONTINUE
          WRITE(4,900) M(I),'Convolution    ',FFALT,FFALTI
        ELSE
          DO 30 J = 0, 2*M(I)
            XJ = -PI+DBLE(J)*DX
            IF(J .LE. M(I)) THEN
              FUKOR = 1.0D0+XJ/PI+SIN(XJ)/PI
            ELSE
              FUKOR = XJ/PI+COS(XJ)+SIN(XJ)/PI
            ENDIF
            FKORR = MAX(FKORR, ABS(F1RE(J)-FUKOR))
            FKORRI= MAX(FKORRI, ABS(F1IM(J)))
  30      CONTINUE
          WRITE(4,900) M(I),'Correlation',FKORR,FKORRI
        ENDIF
  50  CONTINUE
      S T O P
 900  FORMAT(1X,'M = N = ',I3,'  Error at ',A,' : ',2(1PD10.3,2X))
      END
