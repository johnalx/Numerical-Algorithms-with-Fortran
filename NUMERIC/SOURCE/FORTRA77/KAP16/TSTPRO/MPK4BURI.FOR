      PROGRAM TEST
C
C*****************************************************************
C                                                                *
C  Test program for cubature over rectangular regions using the  *
C  Bulirsch-Richardson method.                                   *
C                                                                *
C  We test the subroutine K4BURI.                                *
C                                                                *
C  With our test example we compute the following results:       *
C                                                                *
C[                                                              ]*
C[  EXACT SOLUTION: .38682227139506E+00                         ]*
C[                                                              ]*
C[                                                              ]*
C[  IERR: 0  APPROXIMATE VALUE:  .38682171161932E+00            ]*
C[                                                              ]*
C[  ESTIMATED ERROR: .4797E-09  ACTUAL ERROR:  .5598E-06        ]*
C[                                                              ]*
C[  NUMBER OF FUNCTION CALLS:    164                            ]*
C[                                                              ]*
C                                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  required subroutines: SXCY, K4BURI, K4BUST, BURIEX, DENOM     *
C                                                                *
C*****************************************************************
C                                                                *
C  Author      : Volker KrÅger                                   *
C  Date        : 6.12.1991                                       *
C  Source code : FORTRAN 77                                      *
C                                                                *
C*****************************************************************
C
C Declarations
C
C    N = number of trapezoidal rule cubatures
C
      PARAMETER (N=4)
      EXTERNAL SXCY
      DOUBLE PRECISION WORK(0:N-1,2)
      DOUBLE PRECISION A,B,C,D,CREC,DIVIAT,EXACT
C
C corners of rectangle
C
      A=0.0D+00
      B=1.0D+00
      C=0.0D+00
      D=1.0D+00
C
C number of rectangles in X and Y direction used for first
C trapezoidal cubature
C
      IP=2
      IQ=2
C
C exact solution
C
      EXACT=(COS(B)-COS(A))*(SIN(C)-SIN(D))
      WRITE(*,1000) EXACT
C
C computed solution and error estimate
C
         CALL K4BURI(SXCY,A,B,IP,C,D,IQ,N,CREC,DIVIAT,WORK,
     +               IERR,IUFCLL)
         WRITE(*,1100) IERR,CREC
         WRITE(*,1200) DIVIAT,ABS(CREC-EXACT)
         WRITE(*,1300) IUFCLL
      STOP
C
C Format statements
C
1000  FORMAT (1X,'C[',T66,']*',/,
     +        1X,'C[',2X,'EXACT SOLUTION:',E20.14,T66,']*',/,
     +        1X,'C[',T66,']*')
1100  FORMAT (1X,'C[',T66,']*',/,
     +        1X,'C[',2X,'IERR: ',I1,2X,'APPROXIMATE VALUE: ',
     +        E20.14,T66,']*')
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
