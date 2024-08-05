      PROGRAM TEST
C
C*****************************************************************
C                                                                *
C  Test program for cubature over rectangular region using       *
C  Newton-Cotes formulas:                                        *
C                                                                *
C  We test subroutine  K4NECN                                    *
C                                                                *
C  The test example produces the output:                         *
C                                                                *
C[                                                              ]*
C[  EXACT SOLUTION: .38682227139506E+00                         ]*
C[                                                              ]*
C[                                                              ]*
C[  METHOD: 1  IERR: 0  APPROXIMATE VALUE:  .38581531547004E+00 ]*
C[                                                              ]*
C[  ESTIMATED ERROR: .1005E-02  ACTUAL ERROR: -.1007E-02        ]*
C[                                                              ]*
C[  NUMBER OF FUNCTION CALLS:    106                            ]*
C[                                                              ]*
C[                                                              ]*
C[  METHOD: 2  IERR: 0  APPROXIMATE VALUE:  .38682233700823E+00 ]*
C[                                                              ]*
C[  ESTIMATED ERROR:-.6571E-07  ACTUAL ERROR:  .6561E-07        ]*
C[                                                              ]*
C[  NUMBER OF FUNCTION CALLS:    370                            ]*
C[                                                              ]*
C[                                                              ]*
C[  METHOD: 3  IERR: 0  APPROXIMATE VALUE:  .38682230055496E+00 ]*
C[                                                              ]*
C[  ESTIMATED ERROR:-.2920E-07  ACTUAL ERROR:  .2916E-07        ]*
C[                                                              ]*
C[  NUMBER OF FUNCTION CALLS:    794                            ]*
C[                                                              ]*
C[                                                              ]*
C[  METHOD: 4  IERR: 0  APPROXIMATE VALUE:  .38682227139353E+00 ]*
C[                                                              ]*
C[  ESTIMATED ERROR: .1528E-11  ACTUAL ERROR: -.1525E-11        ]*
C[                                                              ]*
C[  NUMBER OF FUNCTION CALLS:   1378                            ]*
C[                                                              ]*
C[                                                              ]*
C[  METHOD: 5  IERR: 0  APPROXIMATE VALUE:  .38682227139420E+00 ]*
C[                                                              ]*
C[  ESTIMATED ERROR: .8606E-12  ACTUAL ERROR: -.8588E-12        ]*
C[                                                              ]*
C[  NUMBER OF FUNCTION CALLS:   2122                            ]*
C[                                                              ]*
C[                                                              ]*
C[  METHOD: 6  IERR: 0  APPROXIMATE VALUE:  .38682227139506E+00 ]*
C[                                                              ]*
C[  ESTIMATED ERROR:-.3026E-16  ACTUAL ERROR: -.3331E-15        ]*
C[                                                              ]*
C[  NUMBER OF FUNCTION CALLS:   3026                            ]*
C[                                                              ]*
C[                                                              ]*
C[  METHOD: 7  IERR: 0  APPROXIMATE VALUE:  .38682227139506E+00 ]*
C[                                                              ]*
C[  ESTIMATED ERROR:-.1611E-16  ACTUAL ERROR:  .3331E-15        ]*
C[                                                              ]*
C[  NUMBER OF FUNCTION CALLS:   4090                            ]*
C[                                                              ]*
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  required subroutines : SXCY, K4NECN, K4INIT, GRIDOT           *
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
      EXTERNAL SXCY
      DOUBLE PRECISION WORK(0:9)
      DOUBLE PRECISION A,B,C,D,CREC,DIVIAT,EXACT
      LOGICAL ESTDIV
C
C number of rectangles in X and Y direction
C
      IP = 4
      IQ = 4
C
C outer corners of rectangle
C
      A=0.0D+00
      B=1.0D+00
      C=0.0D+00
      D=1.0D+00
C
C set MOLD for first call of  K4NECN
C
      MOLD=8
C
C with error estimation
C
      ESTDIV=.TRUE.
C
C exact solution
C
      EXACT=(COS(B)-COS(A))*(SIN(C)-SIN(D))
      WRITE(*,1000) EXACT
C
C compute approximate solutions for all methods and 
C estimate the error
C
      DO 10 METHOD= 1,7
         CALL K4NECN(SXCY,A,B,IP,C,D,IQ,METHOD,MOLD,CREC,ESTDIV,
     +               DIVIAT,WORK,IERR,IUFCLL)
         WRITE(*,1100) METHOD,IERR,CREC
         WRITE(*,1200) DIVIAT,CREC-EXACT
         WRITE(*,1300) IUFCLL
10    CONTINUE
      STOP
C
C Format statements
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
