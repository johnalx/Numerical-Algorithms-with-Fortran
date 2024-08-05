      PROGRAM TEST
C                                                                              
C*****************************************************************             
C                                                                *             
C  Test program for cubature over triangular regions using the   *
C  N point GAUSS formula:                                        *             
C                                                                *             
C  We test the subroutine K3GAUN .                               *             
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
C[  METHOD: 1  IERR: 0  APPROXIMATE VALUE:  .15062615253663E+00 ]*
C[                                                              ]*
C[  METHOD: 1  IERR: 0  APPROXIMATE VALUE:  .23630357146945E+00 ]*
C[                                                              ]*
C[  APPROXIMATE VALUE:  .38692972400608E+00  ERROR:  .1075E-03  ]*
C[                                                              ]*
C[  NUMBER OF FUNCTION CALLS:    200                            ]*
C[                                                              ]*
C[                                                              ]*
C[  METHOD: 2  IERR: 0  APPROXIMATE VALUE:  .15062615253663E+00 ]*
C[                                                              ]*
C[  METHOD: 2  IERR: 0  APPROXIMATE VALUE:  .23627075228868E+00 ]*
C[                                                              ]*
C[  APPROXIMATE VALUE:  .38689690482531E+00  ERROR:  .7463E-04  ]*
C[                                                              ]*
C[  NUMBER OF FUNCTION CALLS:    400                            ]*
C[                                                              ]*
C[                                                              ]*
C[  METHOD: 3  IERR: 0  APPROXIMATE VALUE:  .15058434344993E+00 ]*
C[                                                              ]*
C[  METHOD: 3  IERR: 0  APPROXIMATE VALUE:  .23623792745072E+00 ]*
C[                                                              ]*
C[  APPROXIMATE VALUE:  .38682227090064E+00  ERROR: -.4944E-09  ]*
C[                                                              ]*
C[  NUMBER OF FUNCTION CALLS:    600                            ]*
C[                                                              ]*
C[                                                              ]*
C[  METHOD: 7  IERR: 0  APPROXIMATE VALUE:  .15058433946991E+00 ]*
C[                                                              ]*
C[  METHOD: 7  IERR: 0  APPROXIMATE VALUE:  .23623793192523E+00 ]*
C[                                                              ]*
C[  APPROXIMATE VALUE:  .38682227139514E+00  ERROR:  .8815E-13  ]*
C[                                                              ]*
C[  NUMBER OF FUNCTION CALLS:   1400                            ]*
C[                                                              ]*
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  required subroutines :    SXCY, K3GAUN                        *
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
      DOUBLE PRECISION A,B,C,D,PX,PY,QX,QY,RX,RY,CTRI,CHELP,EXACT,
     +                 WORK(3,0:6)
C                                                                              
C corners of rectangle
C                                                                              
      A = 0.0D+00
      B = 1.0D+00
      C = 0.0D+00
      D = 1.0D+00
C                                                                              
C initialize MOLD for first call
C                                                                              
      MOLD=-1
C                                                                              
C number of triangles along one rectangular edge
C                                                                              
      N = 10
C                                                                              
C compute exact solution and put out
C                                                                              
      EXACT=(COS(B)-COS(A))*(SIN(C)-SIN(D))
      WRITE(*,1000) EXACT
C                                                                              
C Compute approximations
C                                                                              
      DO 10 METHOD=1,7
         IF (METHOD.LE.3.OR.METHOD.EQ.7) THEN
C                                                                              
C set up coordinates of outer corner of first triangle
C                                                                              
            PX = A
            PY = C
            QX = B
            QY = C
            RX = A
            RY = D
C                                                                              
C compute approximations and print out
C                                                                              
            CALL K3GAUN(SXCY,PX,PY,QX,QY,RX,RY,N,METHOD,MOLD,
     +                  CHELP,WORK,IERR,IUFCLL)
            WRITE(*,1100) METHOD,IERR,CHELP
            CTRI=CHELP
C                                                     
C set up coordinates of outer corner of second triangle                         
C                                                                              
            PX = RX
            PY = RY
            RX = B
            RY = D
C                                                                              
C compute approximate value and put out
C                                                                              
            CALL K3GAUN(SXCY,PX,PY,QX,QY,RX,RY,N,METHOD,MOLD,
     +                  CHELP,WORK,IERR,IHELP)
            WRITE(*,1100) METHOD,IERR,CHELP
C                                                                              
C compute sum of two triangular integrations and put out
C                                                                              
            CTRI=CTRI+CHELP
            IUFCLL=IHELP+IUFCLL
            WRITE(*,1200) CTRI,CTRI-EXACT
            WRITE(*,1300) IUFCLL
         ENDIF
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
C  Test function for cubature
C                                                                *             
C     SXCY = SIN(X)*COS(Y)                                       *             
C                                                                *             
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
