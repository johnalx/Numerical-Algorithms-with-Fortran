      PROGRAM TEST
C***************************************************************************
C                                                                          *
C  Test progrmm for the subroutine   SHPFLG.                               *
C  Computes points on the surface which is generated via the local         *
C  Shepard interpolation with Franke-Little weights.                       *
C                                                                          *
C--------------------------------------------------------------------------*
C                                                                          *
C     The results are displayed on the screen for the following parameters.*
C     For IBM compatible PCs and the MICROSOFT FORTRAN 5.0 Compiler the    *
C     following results are computed:                                      *
C                                                                          *
C   GIVEN:                                                                 *
C   ======                                                                 *
C   RR   = 1.4D0  : Radius, we use the nodes inside this radius to compute *
C                   the new point                                          *
C   DMUE = 3      : Exponent used in a neighborhood of the interpolation   *
C                   point                                                  *
C   N    = 11     : Index of last node                                     *
C                                                                          *
C           FX           |            FY         |           FZ            *
C  ====================================================================    *
C  -.80000000000000D+00  |   .15000000000000D+01 |  .10535653752853D+01    *
C   .20000000000000D+00  |   .13000000000000D+01 |  .15066519173319D+01    *
C   .10000000000000D+01  |   .14000000000000D+01 |  .10198039027186D+01    *
C  -.12000000000000D+01  |   .30000000000000D+00 |  .15716233645502D+01    *
C  -.20000000000000D+00  |   .20000000000000D+00 |  .19798989873223D+01    *
C   .10000000000000D+01  |   .10000000000000D+01 |  .14142135623731D+01    *
C  -.80000000000000D+00  |  -.30000000000000D+00 |  .18083141320025D+01    *
C   .50000000000000D+00  |  -.40000000000000D+00 |  .18947295321496D+01    *
C   .15000000000000D+01  |  -.30000000000000D+00 |  .12884098726725D+01    *
C  -.11000000000000D+01  |  -.14000000000000D+01 |  .91104335791443D+00    *
C   .50000000000000D+00  |  -.17000000000000D+01 |  .92736184954957D+00    *
C   .13000000000000D+01  |  -.14000000000000D+01 |  .59160797830996D+00    *
C                                                                          *
C--------------------------------------------------------------------------*
C                                                                          *
C     Required subroutines: SHPFLG                                         *
C                                                                          *
C--------------------------------------------------------------------------*
C                                                                          *
C[                                                                        ]*
C[  WANTED TO FIND: POINTS THAT LIE ON THE SURFACE GENERATED BY LOCAL     ]*
C[  ==============  SHEPARD INTERPOLATION using FRANKE-LITTLE WEIGHTS     ]*
C[                                                                        ]*
C[                                                                        ]*
C[  SOLUTION: POINTS ON THE INTERPOLATING SURFACE                         ]*
C[  ========                                                              ]*
C[                                                                        ]*
C[                                                                        ]*
C[           X           |          Y           |          Z              ]*
C[  =====================+======================+======================   ]*
C[  -.20000000000000D+01 | -.20000000000000D+01 |  .91104335791443D+00    ]*
C[  -.20000000000000D+01 | -.15000000000000D+01 |  .91104335791443D+00    ]*
C[  -.20000000000000D+01 | -.10000000000000D+01 |  .91105896483399D+00    ]*
C[  -.20000000000000D+01 | -.50000000000000D+00 |  .15752968640045D+01    ]*
C[  -.20000000000000D+01 |  .00000000000000D+00 |  .15777783203363D+01    ]*
C[  -.20000000000000D+01 |  .50000000000000D+00 |  .15716233645502D+01    ]*
C[  -.20000000000000D+01 |  .10000000000000D+01 |  .15584304129828D+01    ]*
C[  -.20000000000000D+01 |  .15000000000000D+01 |  .10535653752853D+01    ]*
C[  -.20000000000000D+01 |  .20000000000000D+01 |  .10535653752853D+01    ]*
C[  -.15000000000000D+01 | -.20000000000000D+01 |  .91104335791443D+00    ]*
C[  -.15000000000000D+01 | -.15000000000000D+01 |  .91104451660051D+00    ]*
C[  -.15000000000000D+01 | -.10000000000000D+01 |  .10065108032420D+01    ]*
C[  -.15000000000000D+01 | -.50000000000000D+00 |  .16173454310709D+01    ]*
C[  -.15000000000000D+01 |  .00000000000000D+00 |  .16235974052526D+01    ]*
C[  -.15000000000000D+01 |  .50000000000000D+00 |  .15769221293391D+01    ]*
C[  -.15000000000000D+01 |  .10000000000000D+01 |  .13764718450711D+01    ]*
C[  -.15000000000000D+01 |  .15000000000000D+01 |  .10600328846103D+01    ]*
C[  -.15000000000000D+01 |  .20000000000000D+01 |  .10535653752853D+01    ]*
C[  -.10000000000000D+01 | -.20000000000000D+01 |  .91104335791443D+00    ]*
C[  -.10000000000000D+01 | -.15000000000000D+01 |  .91381332756624D+00    ]*
C[  -.10000000000000D+01 | -.10000000000000D+01 |  .11261605397241D+01    ]*
C[  -.10000000000000D+01 | -.50000000000000D+00 |  .17240245688862D+01    ]*
C[  -.10000000000000D+01 |  .00000000000000D+00 |  .17126344272204D+01    ]*
C[  -.10000000000000D+01 |  .50000000000000D+00 |  .16176386235409D+01    ]*
C[  -.10000000000000D+01 |  .10000000000000D+01 |  .12371852366067D+01    ]*
C[  -.10000000000000D+01 |  .15000000000000D+01 |  .10570103286106D+01    ]*
C[  -.10000000000000D+01 |  .20000000000000D+01 |  .10535662570312D+01    ]*
C[  -.50000000000000D+00 | -.20000000000000D+01 |  .91450198747067D+00    ]*
C[  -.50000000000000D+00 | -.15000000000000D+01 |  .91966078659882D+00    ]*
C[  -.50000000000000D+00 | -.10000000000000D+01 |  .13318236264308D+01    ]*
C[  -.50000000000000D+00 | -.50000000000000D+00 |  .18162139868411D+01    ]*
C[  -.50000000000000D+00 |  .00000000000000D+00 |  .18654440961765D+01    ]*
C[  -.50000000000000D+00 |  .50000000000000D+00 |  .18367225642311D+01    ]*
C[  -.50000000000000D+00 |  .10000000000000D+01 |  .13469882454728D+01    ]*
C[  -.50000000000000D+00 |  .15000000000000D+01 |  .11378321031252D+01    ]*
C[  -.50000000000000D+00 |  .20000000000000D+01 |  .11044349476463D+01    ]*
C[   .00000000000000D+00 | -.20000000000000D+01 |  .92726730727288D+00    ]*
C[   .00000000000000D+00 | -.15000000000000D+01 |  .93641713370211D+00    ]*
C[   .00000000000000D+00 | -.10000000000000D+01 |  .15246948684704D+01    ]*
C[   .00000000000000D+00 | -.50000000000000D+00 |  .19017423918353D+01    ]*
C[   .00000000000000D+00 |  .00000000000000D+00 |  .19462663531505D+01    ]*
C[   .00000000000000D+00 |  .50000000000000D+00 |  .18987685363452D+01    ]*
C[   .00000000000000D+00 |  .10000000000000D+01 |  .15230471425614D+01    ]*
C[   .00000000000000D+00 |  .15000000000000D+01 |  .14301535272748D+01    ]*
C[   .00000000000000D+00 |  .20000000000000D+01 |  .13866930526380D+01    ]*
C[   .50000000000000D+00 | -.20000000000000D+01 |  .91195808771492D+00    ]*
C[   .50000000000000D+00 | -.15000000000000D+01 |  .90487614412990D+00    ]*
C[   .50000000000000D+00 | -.10000000000000D+01 |  .13859374052259D+01    ]*
C[   .50000000000000D+00 | -.50000000000000D+00 |  .18754619352543D+01    ]*
C[   .50000000000000D+00 |  .00000000000000D+00 |  .18856801190485D+01    ]*
C[   .50000000000000D+00 |  .50000000000000D+00 |  .16307750376598D+01    ]*
C[   .50000000000000D+00 |  .10000000000000D+01 |  .13835745139658D+01    ]*
C[   .50000000000000D+00 |  .15000000000000D+01 |  .13334128767157D+01    ]*
C[   .50000000000000D+00 |  .20000000000000D+01 |  .12805404678918D+01    ]*
C[   .10000000000000D+01 | -.20000000000000D+01 |  .78781806235019D+00    ]*
C[   .10000000000000D+01 | -.15000000000000D+01 |  .70853511489148D+00    ]*
C[   .10000000000000D+01 | -.10000000000000D+01 |  .95970693321142D+00    ]*
C[   .10000000000000D+01 | -.50000000000000D+00 |  .15410475592105D+01    ]*
C[   .10000000000000D+01 |  .00000000000000D+00 |  .15523722764766D+01    ]*
C[   .10000000000000D+01 |  .50000000000000D+00 |  .13823636139404D+01    ]*
C[   .10000000000000D+01 |  .10000000000000D+01 |  .14142135623731D+01    ]*
C[   .10000000000000D+01 |  .15000000000000D+01 |  .11418202252375D+01    ]*
C[   .10000000000000D+01 |  .20000000000000D+01 |  .10912268465121D+01    ]*
C[   .15000000000000D+01 | -.20000000000000D+01 |  .62206258109425D+00    ]*
C[   .15000000000000D+01 | -.15000000000000D+01 |  .60581666154434D+00    ]*
C[   .15000000000000D+01 | -.10000000000000D+01 |  .80163462755432D+00    ]*
C[   .15000000000000D+01 | -.50000000000000D+00 |  .12680019254906D+01    ]*
C[   .15000000000000D+01 |  .00000000000000D+00 |  .13051684319971D+01    ]*
C[   .15000000000000D+01 |  .50000000000000D+00 |  .13355947578821D+01    ]*
C[   .15000000000000D+01 |  .10000000000000D+01 |  .12661692545215D+01    ]*
C[   .15000000000000D+01 |  .15000000000000D+01 |  .11464323984347D+01    ]*
C[   .15000000000000D+01 |  .20000000000000D+01 |  .10538675114362D+01    ]*
C[   .20000000000000D+01 | -.20000000000000D+01 |  .59160797830996D+00    ]*
C[   .20000000000000D+01 | -.15000000000000D+01 |  .59369634531229D+00    ]*
C[   .20000000000000D+01 | -.10000000000000D+01 |  .89051041398281D+00    ]*
C[   .20000000000000D+01 | -.50000000000000D+00 |  .12698037936415D+01    ]*
C[   .20000000000000D+01 |  .00000000000000D+00 |  .12884098726725D+01    ]*
C[   .20000000000000D+01 |  .50000000000000D+00 |  .13119836351226D+01    ]*
C[   .20000000000000D+01 |  .10000000000000D+01 |  .12782003739880D+01    ]*
C[   .20000000000000D+01 |  .15000000000000D+01 |  .11249959143536D+01    ]*
C[   .20000000000000D+01 |  .20000000000000D+01 |  .10198039027186D+01    ]*
C                                                                          *
C***************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FX(0:11),FY(0:11),F(0:11),W(0:11),R(0:11),
     F  EPS(0:11)
      DATA FX /  -.8D+00, .2D+00, .1D+01, -.12D+01, -.2D+00, .1D+01,
     F  -.8D+00, .5D+00, .15D+01, -.11D+01, .5D+00, .13D+01 /
      DATA FY /  .15D+01, .13D+01, .14D+01, .3D+00, .2D+00, .1D+01,
     F  -.3D+00, -.4D+00, -.3D+00, -.14D+01, -.17D+01, -.14D+01 /
      DATA F / .10535653752853D+01, .15066519173319D+01,
     F .10198039027186D+01, .15716233645502D+01, .19798989873223D+01,
     F .14142135623731D+01, .18083141320025D+01, .18947295321496D+01,
     F .12884098726725D+01, .91104335791443D+00, .92736184954957D+00,
     F .59160797830996D+00 /

      DMUE=3.D0
      RR=1.4D0
      N=11
      WRITE(*,890)
      WRITE(*,910)
      DO 10 X=-2.D0,2.D0,0.5D0
        DO 20 Y=-2.D0,2.D0,0.5D0
           CALL SHPFLG(X,Y,FX,FY,F,W,R,EPS,N,DMUE,RR,PHI,IERR)
           IF(IERR.EQ.1) THEN
              WRITE(*,*)IERR, 'SUM of weights is ZERO'
              GOTO 20
           ENDIF
           WRITE(*,900)X,Y,PHI
20      CONTINUE
10    CONTINUE
      STOP
890   FORMAT(1X,'C[',T76,']*',
     +    /,1X,'C[',2X,'WANTED TO FIND: POINTS THAT LIE ON THE SURFACE'
     +    ,' GENERATED BY LOCAL',T76,']*',
     +    /,1X,'C[',2X,'==============  SHEPARD INTERPOLATION using',
     +    ' FRANKE-LITTLE WEIGHTS',T76,']*',
     +    / ,1X,'C[',T76,']*',
     +    / ,1X,'C[',T76,']*',
     +    /,1X,'C[',2X,
     +    'SOLUTION: POINTS ON THE INTERPOLATING SURFACE ',T76,']*',
     +    /,1X,'C[',2X,'========',T76,']*',
     +    /,1X,'C[',T76,']*')
900   FORMAT(1X,'C[',2X,2(D20.14,1X,'|',1X),D20.14,T76,']*')
910   FORMAT(1X,'C[',T76,']*',
     +    /,1X,'C[',11X,'X',11X,'|',10X,'Y',11X,'|',10X,'Z',T76,']*',
     +    /,1X,'C[',2X,21('='),'+',22('='),'+',22('='),T76,']*')
      END
