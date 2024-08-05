      PROGRAM TEST
C
C*****************************************************************
C                                                                *
C  Test program for cubature over rectangular region using the   *
C  Gaussian formulas:                                            *
C                                                                *
C  We test the subroutine K4GAUE                                 *
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
C  required subroutines : SXCY, K4GAUE, K4GINI                   *
C                                                                *
C*****************************************************************
C                                                                *
C  Author      : Volker KrÅger                                   *
C  Date        : 6.12.1991                                       *
C  SOURCE code : FORTRAN 77                                      *
C                                                                *
C*****************************************************************
C
C declarations
C
      EXTERNAL SXCY
      DOUBLE PRECISION WORK(2,0:7)
      DOUBLE PRECISION A,B,C,D,CREC,DIVIAT,EXACT
      LOGICAL ESTDIV
C
C number of rectangles in X and Y direction
C
      IP = 4
      IQ = 4
C
C outer corners
C
      A=0.0D+00
      B=1.0D+00
      C=0.0D+00
      D=1.0D+00
C
C initialize MOLD for first call of K4GAUE
C
      MOLD=-1
C
C set flag to estimate the error
C
      ESTDIV=.TRUE.
C
C exact solution
C
      EXACT=(COS(B)-COS(A))*(SIN(C)-SIN(D))
      WRITE(*,1000) EXACT
C
C compute approximations for all methods and find error estimates
C
      DO 10 METHOD= 0,7
         CALL K4GAUE(SXCY,A,B,IP,C,D,IQ,METHOD,MOLD,CREC,ESTDIV,
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
C  Test function for cubature                                    *
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
