      PROGRAM TEST
C
C*****************************************************************
C                                                                *
C  Test program for cubature over rectangular region using the   *
C  Gaussian formulas:                                            *
C                                                                *
C  We test the subroutine K4GAUV                                 *
C                                                                *
C  Our test example produces the results:                        *
C                                                                *
C[                                                              ]*
C[  EXACT SOLUTION: .38682227139506E+00                         ]*
C[                                                              ]*
C[                                                              ]*
C[  METHOD: 0  IERR: 0  APPROXIMATE VALUE:  .38732633996768E+00 ]*
C[                                                              ]*
C[  ESTIMATED ERROR:-.5056E-03  ACTUAL ERROR:  .5041E-03        ]*
C[                                                              ]*
C[  NUMBER OF FUNCTION CALLS:     80                            ]*
C[                                                              ]*
C[                                                              ]*
C[  METHOD: 1  IERR: 0  APPROXIMATE VALUE:  .38682222765159E+00 ]*
C[                                                              ]*
C[  ESTIMATED ERROR: .2191E-06  ACTUAL ERROR: -.4374E-07        ]*
C[                                                              ]*
C[  NUMBER OF FUNCTION CALLS:    320                            ]*
C[                                                              ]*
C[                                                              ]*
C[  METHOD: 2  IERR: 0  APPROXIMATE VALUE:  .38682227139652E+00 ]*
C[                                                              ]*
C[  ESTIMATED ERROR:-.3081E-10  ACTUAL ERROR:  .1465E-11        ]*
C[                                                              ]*
C[  NUMBER OF FUNCTION CALLS:    720                            ]*
C[                                                              ]*
C[                                                              ]*
C[  METHOD: 3  IERR: 0  APPROXIMATE VALUE:  .38682227139506E+00 ]*
C[                                                              ]*
C[  ESTIMATED ERROR: .2239E-14  ACTUAL ERROR:  .5551E-16        ]*
C[                                                              ]*
C[  NUMBER OF FUNCTION CALLS:   1280                            ]*
C[                                                              ]*
C[                                                              ]*
C[  METHOD: 4  IERR: 0  APPROXIMATE VALUE:  .38682227139506E+00 ]*
C[                                                              ]*
C[  ESTIMATED ERROR:-.1110E-15  ACTUAL ERROR: -.4996E-15        ]*
C[                                                              ]*
C[  NUMBER OF FUNCTION CALLS:   2000                            ]*
C[                                                              ]*
C[                                                              ]*
C[  METHOD: 5  IERR: 0  APPROXIMATE VALUE:  .38682227139506E+00 ]*
C[                                                              ]*
C[  ESTIMATED ERROR: .2035E-15  ACTUAL ERROR:  .4996E-15        ]*
C[                                                              ]*
C[  NUMBER OF FUNCTION CALLS:   2880                            ]*
C[                                                              ]*
C[                                                              ]*
C[  METHOD: 6  IERR: 0  APPROXIMATE VALUE:  .38682227139506E+00 ]*
C[                                                              ]*
C[  ESTIMATED ERROR:-.3331E-15  ACTUAL ERROR: -.3886E-15        ]*
C[                                                              ]*
C[  NUMBER OF FUNCTION CALLS:   3920                            ]*
C[                                                              ]*
C[                                                              ]*
C[  METHOD: 7  IERR: 0  APPROXIMATE VALUE:  .38682227139506E+00 ]*
C[                                                              ]*
C[  ESTIMATED ERROR: .2776E-15  ACTUAL ERROR: -.1110E-15        ]*
C[                                                              ]*
C[  NUMBER OF FUNCTION CALLS:   5120                            ]*
C[                                                              ]*
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: SXCY, K4GAUV, K4GINI                    *
C                                                                *
C*****************************************************************
C                                                                *
C  Author      : Volker KrÅger                                   *
C  Date        : 6.12.1991                                       *
C  Source code : FORTRAN 77                                      *
C                                                                *
C*****************************************************************
C
C declarations
C
C
C   number of rectangles in X nd Y directions
C
      PARAMETER (NX=4,NY=4)
      EXTERNAL SXCY
      DOUBLE PRECISION WORK(2,0:7),X(0:NX),Y(0:NY)
      DOUBLE PRECISION A,B,C,D,DX,DY,CREC,DIVIAT,EXACT
      LOGICAL ESTDIV
C
C outer corners
C
      A=0.0D+00
      B=1.0D+00
      C=0.0D+00
      D=1.0D+00
C
C define X vector
C
      X(0)=A
      X(NX)=B
      DX=(B-A)/DBLE(NX)
      DO 10 I=1,NX-1
         X(I)=X(I-1)+DX
10    CONTINUE
C
C define Y vector
C
      Y(0)=C
      Y(NY)=D
      DY=(D-C)/DBLE(NY)
      DO 20 I=1,NY-1
         Y(I)=Y(I-1)+DY
20    CONTINUE
C
C initialize MOLD for first call of K4GAUV 
C
      MOLD=-1
C
C swith error estimation
C
      ESTDIV=.TRUE.
C
C exact solution
C
      EXACT=(COS(B)-COS(A))*(SIN(C)-SIN(D))
      WRITE(*,1000) EXACT
C
C compute approximate solutions and error estimates
C
      DO 30 METHOD= 0,7
         CALL K4GAUV(SXCY,X,NX,Y,NY,METHOD,MOLD,CREC,ESTDIV,
     +               DIVIAT,WORK,IERR,IUFCLL)
         WRITE(*,1100) METHOD,IERR,CREC
         WRITE(*,1200) DIVIAT,CREC-EXACT
         WRITE(*,1300) IUFCLL
30    CONTINUE
      STOP
C
C Formatangaben
C
1000  FORMAT (1X,'C[',T66,']*',/,
     +        1X,'C[',2X,'EXACT SOLUTION:',E20.14,T66,']*',/,
     +        1X,'C[',T66,']*')
1100  FORMAT (1X,'C[',T66,']*',/,
     +        1X,'C[',2X,'METHOD: ',I1,2X,'IERR: ',I1,2X,
     +        'APPROXIMATE VALUE: ',E20.14,T66,']*')
1200  FORMAT (1X,'C[',T66,']*',/,
     +        1X,'C[',2X,'ESTIMATED ERROR:',E10.4,2X,'ACTUAL ERROR: ',
     +        E10.4,T66,']*')
1300  FORMAT (1X,'C[',T66,']*',/,
     +        1X,'C[',2X,'NUMBER OF FUNCTION CALLS:',I7,T66,']*',/,
     +        1X,'C[',T66,']*')
      END
C
C
      DOUBLE PRECISION FUNCTION SXCY (X,Y)
C
C*****************************************************************
C                                                                *
C  Testfunktion zur Kubatur                                      *
C                                                                *
C     SXCY = SIN(X)*COS(Y)                                       *
C                                                                *
C*****************************************************************
C
      DOUBLE PRECISION X,Y
C
C Z-coordinate
C
      SXCY=SIN(X)*COS(Y)
      RETURN
      END
