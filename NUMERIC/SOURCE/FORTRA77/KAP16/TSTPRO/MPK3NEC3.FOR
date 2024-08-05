      PROGRAM TEST
C
C*****************************************************************
C                                                                *
C  Test program for cubature over triangular regions using the   *
C  NEWTON-COTES formulas:                                        *
C                                                                *             
C  We test the subroutine K3NEC3 .                               *             
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
C[  IERR: 0  APPROXIMATE VALUE:  .15058428623848E+00            ]*
C[                                                              ]*
C[  IERR: 0  APPROXIMATE VALUE:  .23623794484260E+00            ]*
C[                                                              ]*
C[  APPROXIMATE VALUE:  .38682223108107E+00  ERROR: -.4031E-07  ]*
C[                                                              ]*
C[  NUMBER OF FUNCTION CALLS:    330                            ]*
C[                                                              ]*
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  required subroutines :    SXCY, K3NEC3                        *
C                                                                *
C*****************************************************************
C                                                                *
C  Author      : Volker KrÅger                                   *
C  Date        : 6.12.1991                                       *
C  Source code : FORTRAN 77                                      *
C                                                                *
C*****************************************************************
C
C DEclarations
C
      EXTERNAL SXCY
      DOUBLE PRECISION A,B,C,D,PX,PY,QX,QY,RX,RY,CTRI,CHELP,EXACT
C
C corners of rectangle
C
      A = 0.0D+00
      B = 1.0D+00
      C = 0.0D+00
      D = 1.0D+00
C
C number of triangles along one edge
C
      N = 10
C
C find and print exact solution
C
      EXACT=(COS(B)-COS(A))*(SIN(C)-SIN(D))
      WRITE(*,1000) EXACT
C
C set up outer corners of first triangle
C
      PX = A
      PY = C
      QX = B
      QY = C
      RX = A
      RY = D
C
C compute integral 
C
      CALL K3NEC3(SXCY,PX,PY,QX,QY,RX,RY,N,CHELP,IERR,IUFCLL)
      WRITE(*,1100) IERR,CHELP
      CTRI=CHELP
C
C set up outer corners of second triangle
C
      PX = RX
      PY = RY
      RX = B
      RY = D
C
C compute integral
C
      CALL K3NEC3(SXCY,PX,PY,QX,QY,RX,RY,N,CHELP,IERR,IHELP)
      WRITE(*,1100) IERR,CHELP
C
C add results for approximation of rectangular integral
C
      CTRI=CTRI+CHELP
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
     +        1X,'C[',2X,'IERR: ',I1,2X,
     +        'APPROXIMATE VALUE: ',E20.14,T66,']*')
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
