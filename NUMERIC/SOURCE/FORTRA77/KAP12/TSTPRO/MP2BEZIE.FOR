      PROGRAM TEST
C*****************************************************************
C                                                                *
C  Test program for subroutines on the modified bicubic Bezier   *
C  method.                                                       *
C                                                                *
C*****************************************************************
C                                                                *
C  The nodes of the file  "KUGEL" lie on a sphere.               *
C  This test programm interpolates points on a line on a circum- *
C  ference of the sphere.                                        *
C  Results using a MICROSOFT FORTRAN 5.0 compiler:               *
C[                                                              ]*
C[          X                Y                Z                 ]*
C[  ==================================================          ]*
C[                                                              ]*
C[  -.999929933E+00   .353267931E-04  -.294328553E-04           ]*
C[  -.941279558E+00  -.762395581E-06   .324621498E+00           ]*
C[  -.791644770E+00  -.670370507E-06   .614148729E+00           ]*
C[  -.546057400E+00  -.681622303E-06   .837132256E+00           ]*
C[  -.245958139E+00  -.835384915E-06   .969357971E+00           ]*
C[   .827085438E-01  -.662574349E-06   .996555085E+00           ]*
C[   .401896292E+00  -.773519701E-06   .915795150E+00           ]*
C[   .677288911E+00  -.729535461E-06   .735713789E+00           ]*
C[   .879444182E+00  -.659037560E-06   .475847681E+00           ]*
C[   .986314287E+00  -.796146775E-07   .164532612E+00           ]*
C[   .986316568E+00   .757708536E-06  -.164531160E+00           ]*
C[   .879445010E+00   .444048337E-06  -.475850784E+00           ]*
C[   .677286706E+00  -.872998415E-06  -.735710899E+00           ]*
C[   .401901270E+00  -.426938409E-06  -.915796924E+00           ]*
C[   .827049204E-01   .134294660E-05  -.996556971E+00           ]*
C[  -.245956906E+00  -.655975118E-06  -.969356672E+00           ]*
C[  -.546054460E+00  -.528105080E-06  -.837134062E+00           ]*
C[  -.791647792E+00   .886784928E-06  -.614151466E+00           ]*
C[  -.941281043E+00  -.326553655E-06  -.324619041E+00           ]*
C[  -.999985989E+00   .248030080E-04  -.323925255E-05           ]*
C[                                                              ]*
C[  PROGRAM FINISHED                                            ]*
C                                                                *
C*****************************************************************
C                                                                *
C  REQUIRED SUBROUTINES    : BEZIER, INTPOL, BEZPNT, BEZBRD      *
C                            CALCWP, CALCP                       *
C                                                                *
C*****************************************************************
C                                                                *
C  Author      : Hartmut Turowski                                *
C  Date        : 2.10.1991                                       *
C  Source code : FORTRAN 77                                      *
C                                                                *
C*****************************************************************
C..
C..   number of patches in each direction
C..
      PARAMETER (M=16,N=16)
C..
C..   declarations
C..
      DOUBLE PRECISION B(3,0:3*M,0:3*N),D(3,0:M,0:N),
     F                 WERT(3,0:3*M,0:3*N),PUNKTE(3,M*N),
     F                 EPS,WP
C..
C..   open input file
C..      
      OPEN (9,FILE='KUGEL',STATUS='OLD',ERR=1050)
      REWIND (9)
C..
C..   read number of patches
C..
      READ (9,*,END=1100,ERR=1100)IM,IN
C..
C..   too many patches?
C..
      IF ((IM .GT. M) .OR. (IN .GT. N)) THEN
         WRITE (*,100) M,IM,N,IN
100      FORMAT (//1X,'TOO MANY PATCHES'//
     +          1X,'MAXIMAL NUMBER IN M DIRECTION: ',I3,
     +          ' READ IN: ',I3/
     +          1X,'                 AND IN N DIRECTION: ',I3,
     +          ' READ IN: ',I3)
         STOP
      ENDIF
C..
C..   READ IN  Bezier matrix
C..      
      DO 10 I=1,3
         DO 10 K=0,IN
10          READ (9,*,END=1000,ERR=1000)(D(I,J,K),J=0,IM)
C..
C..   close input file
C..
      CLOSE (9)
C..
C..   compute coordinates of the bezier points
C..   ITYP = 1      : modified method
C..   EPS  = 1.0D-3 : error bound
C..
      ITYP = 1
      EPS  = 1.0D-3
      CALL BEZIER(B,D,WERT,ITYP,IM,IN,EPS)
C..
C..   Interpolation of points on the equator
C..   WP = 0.5: Parameter curve on the equator
C..   NZ = 20 : num,ber of points on parameter line
C..
      WP = 0.5D+00
      NZ=20
      CALL CALCVP (B,IM,IN,WP,NZ,PUNKTE)
C..
C..   output
C..
      WRITE (*,2000) '        X                Y                Z '
      DO 20 I=1,NZ
         WRITE (*,3000) (PUNKTE(J,I),J=1,3)
20    CONTINUE
      WRITE (*,4000) 'PROGRAM FINISHED'
      STOP
C..
C..   output of error; stop
C..      
1000  WRITE (*,4000) 'ERROR IN INPUT FILE '
      STOP
1050  WRITE (*,4000) 'ERROR IN OPENING INPUT FILE'
      STOP
1100  WRITE (*,4000) 'EROR WITH CONTENTS OF INPUT FILE KUGEL'
      STOP
2000  FORMAT (1X, 'C[  ', T66, ']*', /,
     +        1X, 'C[  ', A, T66, ']*', /,
     +        1X, 'C[  ', 50('='), T66, ']*', /,
     +        1X, 'C[  ', T66, ']*')
3000  FORMAT (1X, 'C[  ', 3(E15.9,2X), T66, ']*')
4000  FORMAT (1X, 'C[  ', T66, ']*', /,
     +        1X, 'C[  ', A, T66, ']*')
      END
