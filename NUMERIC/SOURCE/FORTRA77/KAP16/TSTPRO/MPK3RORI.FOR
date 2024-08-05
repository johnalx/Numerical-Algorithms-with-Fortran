      PROGRAM TEST
C
C*****************************************************************
C                                                                *
C  Testprogramm zur Kubatur Åber Dreieck-Gebiete nach Newton-    *
C  Cotes-Formeln und 
C  Test program for cubature over triangular regions using the   *
C  Newton-Cotes formulas and Romberg-Richardson extrapolation:   *
C                                                                *             
C  We test the subroutine K3RORI .                               *             
C                                                                *             
C  For F(X,Y)=SIN(X)*COS(Y) we approximate the double integral   *
C  over a rectangular region by dividing the rectangle into two  *
C  triangles. The sum of the two triangular integrals gives the  *
C  value of the rectangular integration.                         *
C                                                                *
C  The test example gives rise to the following results:         *
C                                                                *
C[                                                              ]*
C[  EXACT SOLUTION: .38682227139506E+00                         ]*
C[                                                              ]*
C[                                                              ]*
C[  IERR: 0  APPROXIMATE VALUE:  .15058433969224E+00            ]*
C[           ESTIMATE ERROR: .6282E-08                          ]*
C[                                                              ]*
C[  IERR: 0  APPROXIMATE VALUE:  .23623793189023E+00            ]*
C[           ESTIMATE ERROR: .1155E-08                          ]*
C[                                                              ]*
C[  APPROXIMATE VALUE:  .38682227158247E+00  ERROR:  .1874E-09  ]*
C[                                                              ]*
C[  NUMBER OF FUNCTION CALLS:    300                            ]*
C[                                                              ]*
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines used:   SXCY, K3RORI, RORIEX, K3NEC3              *
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
C   number of cubatures
C
      PARAMETER (N=4)
      EXTERNAL SXCY
      DOUBLE PRECISION A,B,C,D,PX,PY,QX,QY,RX,RY,CTRI,CHILF,
     +                 EXACT,DIVIAT
      DOUBLE PRECISION WORK(0:N-1,2)
C
C corners of rectangle
C
      A = 0.0D+00
      B = 1.0D+00
      C = 0.0D+00
      D = 1.0D+00
C
C exact solution
C
      EXACT=(COS(B)-COS(A))*(SIN(C)-SIN(D))
      WRITE(*,1000) EXACT
C
C coordinates of the outer corners of first triangle
C
      PX = A
      PY = C
      QX = B
      QY = C
      RX = A
      RY = D
C
C compute triangular integral
C
      CALL K3RORI(SXCY,PX,PY,QX,QY,RX,RY,N,WORK,CHILF,DIVIAT,
     +            IERR,IUFCLL)
      WRITE(*,1100) IERR,CHILF,DIVIAT
      CTRI=CHILF
C
C coordinates of the outer corners of second triangle
C
      PX = RX
      PY = RY
      RX = B
      RY = D
C
C compute second triangular integral
C
      CALL K3RORI(SXCY,PX,PY,QX,QY,RX,RY,N,WORK,CHILF,DIVIAT,
     +            IERR,IHELP)
      WRITE(*,1100) IERR,CHILF,DIVIAT
C
C add for rectangular integral approximation
C
      CTRI=CTRI+CHILF
      IUFCLL=IHELP+IUFCLL
      WRITE(*,1200) CTRI,CTRI-EXACT
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
     +        E20.14,T66,']*',/,
     +        1X,'C[',11X,'ESTIMATE ERROR:',E10.4,T66,']*')
1200  FORMAT (1X,'C[',T66,']*',/,
     +        1X,'C[',2X,'APPROXIMATE VALUE: ',E20.14,2X,'ERROR: ',
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
C
      DOUBLE PRECISION X,Y
C
C Z-coordinate
C
      SXCY=SIN(X)*COS(Y)
      RETURN
      END
