      PROGRAM TEST 
!*****************************************************************      
!                                                                *      
!  Test program for subroutines on the modified bicubic Bezier   *      
!  method.                                                       *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  The nodes of the file  "KUGEL" lie on a sphere.               *      
!  This test programm interpolates points on a line on a circum- *      
!  ference of the sphere.                                        *      
!  Results using a MICROSOFT FORTRAN 5.0 compiler:               *      
![                                                              ]*      
![          X                Y                Z                 ]*      
![  ==================================================          ]*      
![                                                              ]*      
![  -.999929933E+00   .353267931E-04  -.294328553E-04           ]*      
![  -.941279558E+00  -.762395581E-06   .324621498E+00           ]*      
![  -.791644770E+00  -.670370507E-06   .614148729E+00           ]*      
![  -.546057400E+00  -.681622303E-06   .837132256E+00           ]*      
![  -.245958139E+00  -.835384915E-06   .969357971E+00           ]*      
![   .827085438E-01  -.662574349E-06   .996555085E+00           ]*      
![   .401896292E+00  -.773519701E-06   .915795150E+00           ]*      
![   .677288911E+00  -.729535461E-06   .735713789E+00           ]*      
![   .879444182E+00  -.659037560E-06   .475847681E+00           ]*      
![   .986314287E+00  -.796146775E-07   .164532612E+00           ]*      
![   .986316568E+00   .757708536E-06  -.164531160E+00           ]*      
![   .879445010E+00   .444048337E-06  -.475850784E+00           ]*      
![   .677286706E+00  -.872998415E-06  -.735710899E+00           ]*      
![   .401901270E+00  -.426938409E-06  -.915796924E+00           ]*      
![   .827049204E-01   .134294660E-05  -.996556971E+00           ]*      
![  -.245956906E+00  -.655975118E-06  -.969356672E+00           ]*      
![  -.546054460E+00  -.528105080E-06  -.837134062E+00           ]*      
![  -.791647792E+00   .886784928E-06  -.614151466E+00           ]*      
![  -.941281043E+00  -.326553655E-06  -.324619041E+00           ]*      
![  -.999985989E+00   .248030080E-04  -.323925255E-05           ]*      
![                                                              ]*      
![  PROGRAM FINISHED                                            ]*      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  REQUIRED SUBROUTINES    : BEZIER, INTPOL, BEZPNT, BEZBRD      *      
!                            CALCWP, CALCP                       *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author      : Hartmut Turowski                                *      
!  Date        : 2.10.1991                                       *      
!  Source code : FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!..                                                                     
!..   number of patches in each direction                               
!..                                                                     
      PARAMETER (M = 16, N = 16) 
!..                                                                     
!..   declarations                                                      
!..                                                                     
      DOUBLEPRECISION B (3, 0:3 * M, 0:3 * N), D (3, 0:M, 0:N), WERT (3,&
      0:3 * M, 0:3 * N), PUNKTE (3, M * N), EPS, WP                     
!..                                                                     
!..   open input file                                                   
!..                                                                     
      OPEN (9, FILE = 'KUGEL', STATUS = 'OLD', ERR = 1050) 
      REWIND (9) 
!..                                                                     
!..   read number of patches                                            
!..                                                                     
      READ (9, *, END = 1100, ERR = 1100) IM, IN 
!..                                                                     
!..   too many patches?                                                 
!..                                                                     
      IF ( (IM.GT.M) .OR. (IN.GT.N) ) THEN 
         WRITE ( *, 100) M, IM, N, IN 
  100 FORMAT    (//1X,'TOO MANY PATCHES'//                              &
     &          1X,'MAXIMAL NUMBER IN M DIRECTION: ',I3,                &
     &          ' READ IN: ',I3/                                        &
     &          1X,'                 AND IN N DIRECTION: ',I3,          &
     &          ' READ IN: ',I3)                                        
         STOP 
      ENDIF 
!..                                                                     
!..   READ IN  Bezier matrix                                            
!..                                                                     
      DO 10 I = 1, 3 
         DO 10 K = 0, IN 
   10 READ (9, *, END = 1000, ERR = 1000) (D (I, J, K), J = 0, IM) 
!..                                                                     
!..   close input file                                                  
!..                                                                     
      CLOSE (9) 
!..                                                                     
!..   compute coordinates of the bezier points                          
!..   ITYP = 1      : modified method                                   
!..   EPS  = 1.0D-3 : error bound                                       
!..                                                                     
      ITYP = 1 
      EPS = 1.0D-3 
      CALL BEZIER (B, D, WERT, ITYP, IM, IN, EPS) 
!..                                                                     
!..   Interpolation of points on the equator                            
!..   WP = 0.5: Parameter curve on the equator                          
!..   NZ = 20 : num,ber of points on parameter line                     
!..                                                                     
      WP = 0.5D+00 
      NZ = 20 
      CALL CALCVP (B, IM, IN, WP, NZ, PUNKTE) 
!..                                                                     
!..   output                                                            
!..                                                                     
      WRITE ( * , 2000) '        X                Y                Z ' 
      DO 20 I = 1, NZ 
         WRITE ( *, 3000) (PUNKTE (J, I), J = 1, 3) 
   20 END DO 
      WRITE ( * , 4000) 'PROGRAM FINISHED' 
      STOP 
!..                                                                     
!..   output of error; stop                                             
!..                                                                     
 1000 WRITE ( * , 4000) 'ERROR IN INPUT FILE ' 
      STOP 
 1050 WRITE ( * , 4000) 'ERROR IN OPENING INPUT FILE' 
      STOP 
 1100 WRITE ( * , 4000) 'EROR WITH CONTENTS OF INPUT FILE KUGEL' 
      STOP 
 2000 FORMAT (1X, 'C[  ', T66, ']*', /,                                 &
     &        1X, 'C[  ', A, T66, ']*', /,                              &
     &        1X, 'C[  ', 50('='), T66, ']*', /,                        &
     &        1X, 'C[  ', T66, ']*')                                    
 3000 FORMAT (1X, 'C[  ', 3(E15.9,2X), T66, ']*') 
 4000 FORMAT (1X, 'C[  ', T66, ']*', /,                                 &
     &        1X, 'C[  ', A, T66, ']*')                                 
      END PROGRAM TEST                              
