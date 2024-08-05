      PROGRAM TEST
C
C*****************************************************************
C----------------------------------------------------------------*
C                 DOUBLE PRECISION Version                       *
C----------------------------------------------------------------*
C                                                                *
C     Test program for the subroutines  BEZIER, CALCWP,          *
C                                       CALCVP and CALCP         *
C----------------------------------------------------------------*
C  required subroutines     : BEZIER, CALCWP, CALCVP, CALCP      *
C----------------------------------------------------------------*
C                                                                *
C  We read in points on the sphere with radius 1 and center at   *
C  (0,0,0).                                                      * 
C  We use the bicubic Bezier method and the modified bicubic     *
C  Bezier method to generate the bezier points for the given     *
C  nodes. Subsequently we interpolate along the parameter lines  *
C  VP=0.5 and WP=0.5.                                            *
C                                                                *
C  Results with the MICROSOFT FORTRAN 5.0 compiler:              *
C                                                                *
C[                                                              ]*
C[                                                              ]*
C[ RESULTS FOR THE MODIFIED BEZIER METHOD                       ]*
C[                                                              ]*
C[  POINTS WITH VP = 0.5                                        ]*
C[  ====================                                        ]*
C[                                                              ]*
C[      .00021295       .99986308      -.00068808               ]*
C[      .00006434       .68105293       .48127421               ]*
C[      .00000997       .36764110       .85182067               ]*
C[      .00000212       .00000608       .99999418               ]*
C[      .00000116      -.44353105       .85184865               ]*
C[     -.00000020      -.83287884       .48148197               ]*
C[     -.00000114     -1.00000323       .00000307               ]*
C[     -.00000110      -.83288142      -.48147852               ]*
C[     -.00000043      -.44353562      -.85185068               ]*
C[      .00000026       .00000072     -1.00000067               ]*
C[      .00000308       .36763341      -.85186954               ]*
C[      .00002102       .68102530      -.48161338               ]*
C[      .00006977       .99977855      -.00044204               ]*
C[                                                              ]*
C[  POINTS WITH WP = 0.5                                        ]*
C[  ====================                                        ]*
C[                                                              ]*
C[     -.99991353       .00090358      -.00016300               ]*
C[     -.93297203      -.25803785      -.00005425               ]*
C[     -.84299278      -.49959550      -.00001620               ]*
C[     -.70699570      -.70698796      -.00001085               ]*
C[     -.50966382      -.86505024      -.00000711               ]*
C[     -.26633273      -.96508299      -.00000138               ]*
C[     -.00000114     -1.00000323       .00000307               ]*
C[      .26633187      -.96508878       .00000384               ]*
C[      .50966585      -.86506049       .00000212               ]*
C[      .70700000      -.70699999      -.00000001               ]*
C[      .84300234      -.49964594       .00000030               ]*
C[      .93301618      -.25836517       .00000910               ]*
C[     1.00005365      -.00018156       .00003320               ]*
C[                                                              ]*
C[                                                              ]*
C[ RESULTS FOR THE BEZIER METHOD                                ]*
C[                                                              ]*
C[  POINTS WITH VP = 0.5                                        ]*
C[  ====================                                        ]*
C[                                                              ]*
C[      .00000000       .90233333       .00000000               ]*
C[      .00000000       .60155556       .28963786               ]*
C[      .00000000       .30077778       .51243621               ]*
C[      .00000000       .00000000       .60155556               ]*
C[      .00000000      -.28963786       .51243621               ]*
C[      .00000000      -.51243621       .28963786               ]*
C[      .00000000      -.60155556       .00000000               ]*
C[      .00000000      -.51243621      -.28963786               ]*
C[      .00000000      -.28963786      -.51243621               ]*
C[      .00000000       .00000000      -.60155556               ]*
C[      .00000000       .30077778      -.51243621               ]*
C[      .00000000       .60155556      -.28963786               ]*
C[      .00000000       .90233333       .00000000               ]*
C[                                                              ]*
C[  POINTS WITH WP = 0.5                                        ]*
C[  ====================                                        ]*
C[                                                              ]*
C[    -1.00000000       .00000000       .00000000               ]*
C[     -.89977778      -.15540741       .00000000               ]*
C[     -.78422222      -.30059259       .00000000               ]*
C[     -.63800000      -.42533333       .00000000               ]*
C[     -.45088889      -.52040329       .00000000               ]*
C[     -.23311111      -.58055967       .00000000               ]*
C[      .00000000      -.60155556       .00000000               ]*
C[      .23311111      -.58055967       .00000000               ]*
C[      .45088889      -.52040329       .00000000               ]*
C[      .63800000      -.42533333       .00000000               ]*
C[      .78422222      -.30059259       .00000000               ]*
C[      .89977778      -.15540741       .00000000               ]*
C[     1.00000000       .00000000       .00000000               ]*
C                                                                *
C                                                                *
C                                                                *
C*****************************************************************
C                                                                *
C  Author      : Hartmut Turowski                                *
C  Date        : 5.18.1990                                       *
C  Source code : FORTRAN 77                                      *
C                                                                *
C*****************************************************************
C..
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C..
C..   number of patches in each direction
C..
      PARAMETER (M=4,N=4,IANZ=13)
C..
C..   declarations for Bezier
C..
      DIMENSION B(3,0:3*M,0:3*N),D(3,0:M,0:N),
     F          WERT(3,0:M,0:N),
     F          PUNKTE (3,IANZ)
C..
C..   declarations for modified Bezier
C..
      DIMENSION BM(3,0:3*M,0:3*N),DM(3,0:M,0:N)
C..
C..   initialize for modified Bezier
C..
C..   initialize left edge (semi circle in  X-Y plane)
C..
      DATA ((DM(K,0,I),K=1,3),I=0,N)
C..                                     X        Y        Y
     F                               /-1.0D0,   0.0D0,   0.0D0,
     F                                -0.707D0, 0.707D0, 0.0D0,
     F                                 0.0D0,   1.0D0,   0.0D0,
     F                                 0.707D0, 0.707D0, 0.0D0,
     F                                 1.0D0,   0.0D0,   0.0D0/
C..
C..   initialize right edge (identical semi circle in  X-Y plane)
C..   (full sphere from 0 degrees to 360 degrees)
C..
      DATA ((DM(K,M,I),K=1,3),I=0,N)
C..                                       X        Y        Y
     F                                 /-1.0D0,   0.0D0,   0.0D0,
     F                                  -0.707D0, 0.707D0, 0.0D0,
     F                                   0.0D0,   1.0D0,   0.0D0,
     F                                   0.707D0, 0.707D0, 0.0D0,
     F                                   1.0D0,   0.0D0,   0.0D0/
C..
C..   initialize upper boundary of sphere (north pole)
C..
      DATA ((DM(K,J,0),J=1,M-1),K=1,3)
C..                                      X         Y         Y
     F                             / 3*-1.0D0,  3*0.0D0,  3*0.0D0/
C..
C..   ditto south pole
C..
      DATA ((DM(K,J,N),J=1,M-1),K=1,3)
C..                                      X         Y         Z
     F                             / 3*+1.0D0,  3*0.0D0,  3*0.0D0/
C..
C..   initialize interpolation points
C..
      DATA (((DM(K,J,I),K=1,3),I=1,M-1),J=1,N-1)
C..                                   X         Y        Z
     F                             /-0.707D0,  0.0D0,   0.707D0,
     F                               0.0D0,    0.0D0,   1.0D0,
     F                               0.707D0,  0.0D0,   0.707D0,
     F                              -0.707D0, -0.707D0, 0.0D0,
     F                               0.0D0,   -1.0D0,   0.0D0,
     F                               0.707D0, -0.707D0, 0.0D0,
     F                              -0.707D0,  0.0D0,  -0.707D0,
     F                               0.0D0,    0.0D0,  -1.0D0,
     F                               0.707D0,  0.0D0,  -0.707D0/
C..
C..   set up array for Bezier
C..
       DO 30 I = 0,N
          DO 30 J = 0,M
             DO 30 K = 1,3
                D(K,J,I) =  DM(K,J,I)
30     CONTINUE
C..
C..   coordinates of the Bezier points (modified method)
C..

      ITYP = 1
      EPS  = 1.0D-3
      CALL BEZIER(BM,DM,WERT,ITYP,M,N,EPS)
C..
C..   compute coordinates of Bezier points
C..
      ITYP = 0
      CALL BEZIER(B,D,WERT,ITYP,M,N,EPS)
C..
C..   compute intermediate points along the parameter lines
C..   WP = 0.5 and VP = 0.5.
C..   VP = 0.5 corresponds to a constant x value of 0.
C..   (X from -1 to +1) We interpolate a circle in the Y-Z plane
C.. 
      WRITE (*,100)' MODIFIED '
      WRITE (*,200)'POINTS WITH VP = 0.5'
      VP = 0.5D0
      CALL CALCVP (BM,M,N,VP,IANZ,PUNKTE)
      DO 40 I = 1, IANZ
         WRITE (*,300)(PUNKTE (K,I),K=1,3)
40    CONTINUE
C..
C..   for WP = 0.5 we interpolate a semicircle in the 3rd and 4th
C..   quadrant of the X-Y plane
C.. 
      WRITE (*,200)'POINTS WITH WP = 0.5'
      WP = 0.5D0
      CALL CALCWP (BM,M,N,WP,IANZ,PUNKTE)
      DO 50 I = 1, IANZ
         WRITE (*,300)(PUNKTE (K,I),K=1,3)
50    CONTINUE
C..
C..   results for the bicubic method
C..
      WRITE (*,100)' '
      WRITE (*,200)'POINTS WITH VP = 0.5'
      VP = 0.5D0
      CALL CALCVP (B,M,N,VP,IANZ,PUNKTE)
      DO 60 I = 1, IANZ
         WRITE (*,300)(PUNKTE (K,I),K=1,3)
60    CONTINUE
      WRITE (*,200)'POINTS WITH WP = 0.5'
      WP = 0.5D0
      CALL CALCWP (B,M,N,WP,IANZ,PUNKTE)
      DO 70 I = 1, IANZ
         WRITE (*,300)(PUNKTE (K,I),K=1,3)
70    CONTINUE
      STOP
100   FORMAT (1X, 'C[', T66, ']*', /,
     +        1X, 'C[', T66, ']*', /,
     +        1X, 'C[ RESULTS FOR THE', A, 'BEZIER METHOD',
     +            T66, ']*')
200   FORMAT (1X, 'C[', T66, ']*', /,
     +        1X, 'C[  ', A, T66, ']*', /,
     +        1X, 'C[  ', 20('='), T66, ']*', /,
     +        1X, 'C[', T66, ']*')
300   FORMAT (1X, 'C[ ', 3(F14.8,2X), T66, ']*')
      END
