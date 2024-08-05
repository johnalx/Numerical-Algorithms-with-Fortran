      PROGRAM TEST 
!                                                                       
!*****************************************************************      
!----------------------------------------------------------------*      
!                 DOUBLE PRECISION Version                       *      
!----------------------------------------------------------------*      
!                                                                *      
!     Test program for the subroutines  BEZIER, CALCWP,          *      
!                                       CALCVP and CALCP         *      
!----------------------------------------------------------------*      
!  required subroutines     : BEZIER, CALCWP, CALCVP, CALCP      *      
!----------------------------------------------------------------*      
!                                                                *      
!  We read in points on the sphere with radius 1 and center at   *      
!  (0,0,0).                                                      *      
!  We use the bicubic Bezier method and the modified bicubic     *      
!  Bezier method to generate the bezier points for the given     *      
!  nodes. Subsequently we interpolate along the parameter lines  *      
!  VP=0.5 and WP=0.5.                                            *      
!                                                                *      
!  Results with the MICROSOFT FORTRAN 5.0 compiler:              *      
!                                                                *      
![                                                              ]*      
![                                                              ]*      
![ RESULTS FOR THE MODIFIED BEZIER METHOD                       ]*      
![                                                              ]*      
![  POINTS WITH VP = 0.5                                        ]*      
![  ====================                                        ]*      
![                                                              ]*      
![      .00021295       .99986308      -.00068808               ]*      
![      .00006434       .68105293       .48127421               ]*      
![      .00000997       .36764110       .85182067               ]*      
![      .00000212       .00000608       .99999418               ]*      
![      .00000116      -.44353105       .85184865               ]*      
![     -.00000020      -.83287884       .48148197               ]*      
![     -.00000114     -1.00000323       .00000307               ]*      
![     -.00000110      -.83288142      -.48147852               ]*      
![     -.00000043      -.44353562      -.85185068               ]*      
![      .00000026       .00000072     -1.00000067               ]*      
![      .00000308       .36763341      -.85186954               ]*      
![      .00002102       .68102530      -.48161338               ]*      
![      .00006977       .99977855      -.00044204               ]*      
![                                                              ]*      
![  POINTS WITH WP = 0.5                                        ]*      
![  ====================                                        ]*      
![                                                              ]*      
![     -.99991353       .00090358      -.00016300               ]*      
![     -.93297203      -.25803785      -.00005425               ]*      
![     -.84299278      -.49959550      -.00001620               ]*      
![     -.70699570      -.70698796      -.00001085               ]*      
![     -.50966382      -.86505024      -.00000711               ]*      
![     -.26633273      -.96508299      -.00000138               ]*      
![     -.00000114     -1.00000323       .00000307               ]*      
![      .26633187      -.96508878       .00000384               ]*      
![      .50966585      -.86506049       .00000212               ]*      
![      .70700000      -.70699999      -.00000001               ]*      
![      .84300234      -.49964594       .00000030               ]*      
![      .93301618      -.25836517       .00000910               ]*      
![     1.00005365      -.00018156       .00003320               ]*      
![                                                              ]*      
![                                                              ]*      
![ RESULTS FOR THE BEZIER METHOD                                ]*      
![                                                              ]*      
![  POINTS WITH VP = 0.5                                        ]*      
![  ====================                                        ]*      
![                                                              ]*      
![      .00000000       .90233333       .00000000               ]*      
![      .00000000       .60155556       .28963786               ]*      
![      .00000000       .30077778       .51243621               ]*      
![      .00000000       .00000000       .60155556               ]*      
![      .00000000      -.28963786       .51243621               ]*      
![      .00000000      -.51243621       .28963786               ]*      
![      .00000000      -.60155556       .00000000               ]*      
![      .00000000      -.51243621      -.28963786               ]*      
![      .00000000      -.28963786      -.51243621               ]*      
![      .00000000       .00000000      -.60155556               ]*      
![      .00000000       .30077778      -.51243621               ]*      
![      .00000000       .60155556      -.28963786               ]*      
![      .00000000       .90233333       .00000000               ]*      
![                                                              ]*      
![  POINTS WITH WP = 0.5                                        ]*      
![  ====================                                        ]*      
![                                                              ]*      
![    -1.00000000       .00000000       .00000000               ]*      
![     -.89977778      -.15540741       .00000000               ]*      
![     -.78422222      -.30059259       .00000000               ]*      
![     -.63800000      -.42533333       .00000000               ]*      
![     -.45088889      -.52040329       .00000000               ]*      
![     -.23311111      -.58055967       .00000000               ]*      
![      .00000000      -.60155556       .00000000               ]*      
![      .23311111      -.58055967       .00000000               ]*      
![      .45088889      -.52040329       .00000000               ]*      
![      .63800000      -.42533333       .00000000               ]*      
![      .78422222      -.30059259       .00000000               ]*      
![      .89977778      -.15540741       .00000000               ]*      
![     1.00000000       .00000000       .00000000               ]*      
!                                                                *      
!                                                                *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author      : Hartmut Turowski                                *      
!  Date        : 5.18.1990                                       *      
!  Source code : FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!..                                                                     
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!..                                                                     
!..   number of patches in each direction                               
!..                                                                     
      PARAMETER (M = 4, N = 4, IANZ = 13) 
!..                                                                     
!..   declarations for Bezier                                           
!..                                                                     
      DIMENSION B (3, 0:3 * M, 0:3 * N), D (3, 0:M, 0:N), WERT (3, 0:M, &
      0:N), PUNKTE (3, IANZ)                                            
!..                                                                     
!..   declarations for modified Bezier                                  
!..                                                                     
      DIMENSION BM (3, 0:3 * M, 0:3 * N), DM (3, 0:M, 0:N) 
!..                                                                     
!..   initialize for modified Bezier                                    
!..                                                                     
!..   initialize left edge (semi circle in  X-Y plane)                  
!..                                                                     
      DATA ( (DM (K, 0, I), K = 1, 3), I = 0, N) / - 1.0D0, 0.0D0,      &
      0.0D0, - 0.707D0, 0.707D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.707D0,   &
      0.707D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0 /                             
!..                                     X        Y        Y             
!..                                                                     
!..   initialize right edge (identical semi circle in  X-Y plane)       
!..   (full sphere from 0 degrees to 360 degrees)                       
!..                                                                     
      DATA ( (DM (K, M, I), K = 1, 3), I = 0, N) / - 1.0D0, 0.0D0,      &
      0.0D0, - 0.707D0, 0.707D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.707D0,   &
      0.707D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0 /                             
!..                                       X        Y        Y           
!..                                                                     
!..   initialize upper boundary of sphere (north pole)                  
!..                                                                     
      DATA ( (DM (K, J, 0), J = 1, M - 1), K = 1, 3) / 3 * - 1.0D0, 3 * &
      0.0D0, 3 * 0.0D0 /                                                
!..                                      X         Y         Y          
!..                                                                     
!..   ditto south pole                                                  
!..                                                                     
      DATA ( (DM (K, J, N), J = 1, M - 1), K = 1, 3) / 3 * + 1.0D0, 3 * &
      0.0D0, 3 * 0.0D0 /                                                
!..                                      X         Y         Z          
!..                                                                     
!..   initialize interpolation points                                   
!..                                                                     
      DATA ( ( (DM (K, J, I), K = 1, 3), I = 1, M - 1), J = 1, N - 1)   &
      / - 0.707D0, 0.0D0, 0.707D0, 0.0D0, 0.0D0, 1.0D0, 0.707D0, 0.0D0, &
      0.707D0, - 0.707D0, - 0.707D0, 0.0D0, 0.0D0, - 1.0D0, 0.0D0,      &
      0.707D0, - 0.707D0, 0.0D0, - 0.707D0, 0.0D0, - 0.707D0, 0.0D0,    &
      0.0D0, - 1.0D0, 0.707D0, 0.0D0, - 0.707D0 /                       
!..                                   X         Y        Z              
!..                                                                     
!..   set up array for Bezier                                           
!..                                                                     
      DO 30 I = 0, N 
         DO 30 J = 0, M 
            DO 30 K = 1, 3 
               D (K, J, I) = DM (K, J, I) 
   30 CONTINUE 
!..                                                                     
!..   coordinates of the Bezier points (modified method)                
!..                                                                     
                                                                        
      ITYP = 1 
      EPS = 1.0D-3 
      CALL BEZIER (BM, DM, WERT, ITYP, M, N, EPS) 
!..                                                                     
!..   compute coordinates of Bezier points                              
!..                                                                     
      ITYP = 0 
      CALL BEZIER (B, D, WERT, ITYP, M, N, EPS) 
!..                                                                     
!..   compute intermediate points along the parameter lines             
!..   WP = 0.5 and VP = 0.5.                                            
!..   VP = 0.5 corresponds to a constant x value of 0.                  
!..   (X from -1 to +1) We interpolate a circle in the Y-Z plane        
!..                                                                     
      WRITE ( * , 100) ' MODIFIED ' 
      WRITE ( * , 200) 'POINTS WITH VP = 0.5' 
      VP = 0.5D0 
      CALL CALCVP (BM, M, N, VP, IANZ, PUNKTE) 
      DO 40 I = 1, IANZ 
         WRITE ( *, 300) (PUNKTE (K, I), K = 1, 3) 
   40 END DO 
!..                                                                     
!..   for WP = 0.5 we interpolate a semicircle in the 3rd and 4th       
!..   quadrant of the X-Y plane                                         
!..                                                                     
      WRITE ( * , 200) 'POINTS WITH WP = 0.5' 
      WP = 0.5D0 
      CALL CALCWP (BM, M, N, WP, IANZ, PUNKTE) 
      DO 50 I = 1, IANZ 
         WRITE ( *, 300) (PUNKTE (K, I), K = 1, 3) 
   50 END DO 
!..                                                                     
!..   results for the bicubic method                                    
!..                                                                     
      WRITE ( * , 100) ' ' 
      WRITE ( * , 200) 'POINTS WITH VP = 0.5' 
      VP = 0.5D0 
      CALL CALCVP (B, M, N, VP, IANZ, PUNKTE) 
      DO 60 I = 1, IANZ 
         WRITE ( *, 300) (PUNKTE (K, I), K = 1, 3) 
   60 END DO 
      WRITE ( * , 200) 'POINTS WITH WP = 0.5' 
      WP = 0.5D0 
      CALL CALCWP (B, M, N, WP, IANZ, PUNKTE) 
      DO 70 I = 1, IANZ 
         WRITE ( *, 300) (PUNKTE (K, I), K = 1, 3) 
   70 END DO 
      STOP 
  100 FORMAT (1X, 'C[', T66, ']*', /,                                   &
     &        1X, 'C[', T66, ']*', /,                                   &
     &        1X, 'C[ RESULTS FOR THE', A, 'BEZIER METHOD',             &
     &            T66, ']*')                                            
  200 FORMAT (1X, 'C[', T66, ']*', /,                                   &
     &        1X, 'C[  ', A, T66, ']*', /,                              &
     &        1X, 'C[  ', 20('='), T66, ']*', /,                        &
     &        1X, 'C[', T66, ']*')                                      
  300 FORMAT (1X, 'C[ ', 3(F14.8,2X), T66, ']*') 
      END PROGRAM TEST                              
