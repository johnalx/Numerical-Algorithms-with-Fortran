      PROGRAM TEST 
!                                                                       
!***********************************************************************
!                Test program for the subroutines                      *
!                    ISPLTR, ISPLPE, TSPTAB                            *
!----------------------------------------------------------------------*
!     required subroutines    :  ISPLTR, ISPLPE, TSPTAB,               *
!                                CYTSY,  CYTSYP, CYTSYS                *
!----------------------------------------------------------------------*
!                                                                      *
!     This test program computes a transformed interpolating cubic     *
!     spline.                                                          *
!                                                                      *
![  GIVEN:   N + 1 POINTS X(I),Y(I), I=0(1)N, (N=12).                 ]*
![  ======                                                            ]*
![                                                                    ]*
![   I        X(I)         Y(I)                                       ]*
![  ----------------------------                                      ]*
![   0      32.0000      15.0000                                      ]*
![   1      16.0000      32.0000                                      ]*
![   2      -6.0000      31.0000                                      ]*
![   3     -18.0000      18.0000                                      ]*
![   4     -20.0000      10.0000                                      ]*
![   5     -20.0000      -7.0000                                      ]*
![   6     -13.0000     -22.0000                                      ]*
![   7      -4.0000     -26.0000                                      ]*
![   8       4.0000     -31.0000                                      ]*
![   9      12.0000     -28.0000                                      ]*
![  10      21.0000     -25.0000                                      ]*
![  11      33.0000     -11.0000                                      ]*
![  12      32.0000      15.0000                                      ]*
![                                                                    ]*
![                                                                    ]*
![  TO FIND:   A) TABLE OF THE TRANSFORMED COORDINATES PHI(I), R(I),  ]*
![  ========      I=0(1)N.                                            ]*
![                                                                    ]*
![             B) COEFFICIENTS OF A TRANSFORMED INTERPOLATING CUBIC   ]*
![                SPLINE S(PHI) FOR THE NODES PHI(I), R(I)            ]*
![                                                                    ]*
![             C) TABLE OF VALUES WITH XW=XW(S(PHI)), YW=YW(S(PHI))   ]*
![                                                                    ]*
![                                                                    ]*
![  SOLUTION:                                                         ]*
![  =========                                                         ]*
![                                                                    ]*
![  A) TRANSFORMED COORDINATES:                                       ]*
![     ------------------------                                       ]*
![                                                                    ]*
![  SHIFT IN X DIRECTION: PX =   6.500000                             ]*
![  SHIFT IN Y DIRECTION: PY =    .500000                             ]*
![  ROTATION OF COORDINATE SYSTEM (IN RAD): PHID =  .51703195D+00     ]*
![                                                                    ]*
![  TRANSFORMED NODES:                                                ]*
![                                                                    ]*
![   I         PHI(I)                  R(I)                           ]*
![  ------------------------------------------------                  ]*
![   0    .00000000000000D+00    .29334280287745D+02                  ]*
![   1    .76085197677783D+00    .32901367752724D+02                  ]*
![   2    .14427212544824D+01    .32962099447699D+02                  ]*
![   3    .20043112151154D+01    .30108138434649D+02                  ]*
![   4    .22803420181158D+01    .28151376520518D+02                  ]*
![   5    .29003665996868D+01    .27540878707841D+02                  ]*
![   6    .34812663292810D+01    .29774149861919D+02                  ]*
![   7    .38181079682508D+01    .28504385627478D+02                  ]*
![   8    .41161579569237D+01    .31599050618650D+02                  ]*
![   9    .43859959338072D+01    .29025850547400D+02                  ]*
![  10    .47123889803847D+01    .29334280287745D+02                  ]*
![  11    .53567161451619D+01    .28887713651309D+02                  ]*
![  12    .62831853071796D+01    .29334280287745D+02                  ]*
![                                                                    ]*
![                                                                    ]*
![  B) SPLINE COEFFICIENTS:                                           ]*
![     --------------------                                           ]*
![                                                                    ]*
![   I        A(I)           B(I)           C(I)           D(I)       ]*
![  --------------------------------------------------------------    ]*
![   0    .293343D+02    .377997D+01    .425528D+01   -.402375D+01    ]*
![   1    .329014D+02    .326726D+01   -.492915D+01    .393250D+00    ]*
![   2    .329621D+02   -.290629D+01   -.412471D+01    .446295D+00    ]*
![   3    .301081D+02   -.711682D+01   -.337281D+01    .125850D+02    ]*
![   4    .281514D+02   -.610214D+01    .704877D+01    .194339D+01    ]*
![   5    .275409D+02    .487997D+01    .106636D+02   -.214256D+02    ]*
![   6    .297741D+02   -.442092D+01   -.266748D+02    .849311D+02    ]*
![   7    .285044D+02    .651816D+01    .591502D+02   -.154950D+03    ]*
![   8    .315991D+02    .483085D+00   -.793987D+02    .156644D+03    ]*
![   9    .290259D+02   -.814966D+01    .474064D+02   -.598738D+02    ]*
![  10    .293343D+02    .366108D+01   -.112208D+02    .692679D+01    ]*
![  11    .288877D+02   -.217152D+01    .216856D+01    .750782D+00    ]*
!                                                                      *
!     The transformed coordinates PHIN(I), R(I), the shift coordinates *
!     PX, PY and the rotation of the system PHID are computed in ISPLTR*
!     These, the spline coefficients and the table of values for the   *
!     points XW(I), YW(I) on the curve are put onto TAPE3.             *
!                                                                      *
!***********************************************************************
!                                                                      *
!     Author      :  GÅnter Palm                                       *
!     Date        :  5.20.1988                                         *
!     Source code :  FORTRAN 77                                        *
!                                                                      *
!***********************************************************************
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      PARAMETER (N = 12, MT = 100) 
!                                                                       
      DIMENSION X (0:N), Y (0:N), PHIN (0:N), R (0:N) 
      DIMENSION B (0:N), C (0:N), D (0:N) 
      DIMENSION XW (0:N + MT + 2), YW (0:N + MT + 2) 
      DIMENSION HILF (5 * N + 1) 
!                                                                       
      DATA X / 32.0D0, 16.0D0, - 6.0D0, - 18.0D0, - 20.0D0, - 20.0D0,   &
      - 13.0D0, - 4.0D0, 4.0D0, 12.0D0, 21.0D0, 33.0D0, 32.0D0 /        
      DATA Y / 15.0D0, 32.0D0, 31.0D0, 18.0D0, 10.0D0, - 7.0D0, -       &
      22.0D0, - 26.0D0, - 31.0D0, - 28.0D0, - 25.0D0, - 11.0D0, 15.0D0 /
      DATA JV / - 1 / 
!                                                                       
      CALL ISPLTR (N, X, Y, JV, PX, PY, R, B, C, D, PHIN, PHID, HILF,   &
      IFEHL)                                                            
      IF (IFEHL.NE.0) THEN 
         WRITE ( *, 1000) IFEHL 
         STOP 
      ENDIF 
      OPEN (3, FILE = 'TAPE3') 
      REWIND (3) 
      WRITE (3, 2000) N 
      DO 100 I = 0, N, 1 
         WRITE (3, 2100) I, X (I), Y (I) 
  100 END DO 
      WRITE (3, 2200) 
      WRITE (3, 2250) PX, PY, PHID 
      DO 150 I = 0, N, 1 
         WRITE (3, 2300) I, PHIN (I), R (I) 
  150 END DO 
      WRITE (3, 2400) 
      DO 200 I = 0, N - 1, 1 
         WRITE (3, 2500) I, R (I), B (I), C (I), D (I) 
  200 END DO 
!                                                                       
      WRITE ( *, 1050) 
   10 WRITE ( *, 1051) 
      READ ( *, *, END = 999) MARGIN 
      IF (MARGIN.EQ.1) THEN 
         GOTO 999 
      ELSEIF (MARGIN.EQ.2) THEN 
         PANF = PHIN (0) 
         PEND = PHIN (N) 
      ELSEIF (MARGIN.EQ.3) THEN 
         WRITE ( *, 1100) 
         READ ( *, *, END = 999) PANF 
         WRITE ( *, 1200) 
         READ ( *, *, END = 999) PEND 
      ELSE 
         GOTO 10 
      ENDIF 
      CALL TSPTAB (N, MT, PANF, PEND, PHIN, R, B, C, D, PHID, PX, PY,   &
      NT, XW, YW, IFEHL)                                                
      IF (IFEHL.NE.0) THEN 
         WRITE ( *, 1300) 
         STOP 
      ENDIF 
      WRITE (3, 2700) 
      DO 300 I = 0, NT, 1 
         WRITE (3, 2800) I, XW (I), YW (I) 
  300 END DO 
      WRITE ( *, 2900) 
  999 STOP 
 1000 FORMAT(' ERROR WHEN COMPOUTING COEFFICIENTS.',/,                  &
     &       ' ERROR PARAMETER =',I4)                                   
 1050 FORMAT(///,' COEFFICIENTS COMPUTED.',/,                           &
     &' RESULTS ON TAPE3.',//,' TABLE OF VALUES:',/,                    &
     &' EXIT 1: NO TABLE OF VALUES (PROGRAM STOP)',/,                   &
     &' EXIT 2: END POINTS SET BY PROGRAM',/,                           &
     &'         STARTING VALUE = PHIN(0) = 0.0',/,                      &
     &'         FINAL VALUE    = PHIN(N) = 2*PI',/,                     &
     &' EXIT 3: SET UP TABLE OF VALUES')                                
 1051 FORMAT(' SECLECT BY ENTERING 1, 2 OR 3: ') 
 1100 FORMAT(' STARTING VALUE (PANF) ? ') 
 1200 FORMAT(' FINAL VALUE (PEND) ? ') 
 1300 FORMAT(' ERRONEOUS INPUT.',/,                                     &
     &       ' PANF EXCEEDS PEND.')                                     
 2000 FORMAT('C[  GIVEN:   N + 1 POINTS X(I),Y(I), I=0(1)N, ',          &
     &          '(N=',I2,').',T71,']*',/,                               &
     &       'C[  ======',T71,']*',/,                                   &
     &       'C[',T71,']*',/,                                           &
     &       'C[',3X,'I',8X,'X(I)',9X,'Y(I)',T71,']*',/,                &
     &       'C[  ',28('-'),T71,']*')                                   
 2100 FORMAT('C[  ',I2,4X,F9.4,4X,F9.4,T71,']*') 
 2200 FORMAT('C[',T71,']*',/,'C[',T71,']*',/,                           &
     &       'C[  TO FIND:   A) TABLE OF THE TRANSFORMED ',             &
     &          'COORDINATES PHI(I), R(I),',T71,']*',/,                 &
     &       'C[  ========',6X,'I=0(1)N.',T71,']*',/,                   &
     &       'C[',T71,']*',/,                                           &
     &       'C[',13X,'B) COEFFICIENTS OF A TRANSFORMED ',              &
     &          'INTERPOLATING CUBIC',T71,']*',/,                       &
     &       'C[',16X,'SPLINE S(PHI) FOR THE NODES PHI(I), R(I)',       &
     &          T71,']*',/,                                             &
     &       'C[',T71,']*',/,                                           &
     &       'C[',13X,'C) TABLE OF VALUES WITH ',                       &
     &          'XW=XW(S(PHI)), YW=YW(S(PHI))',T71,']*',/,              &
     &       'C[',T71,']*',/,'C[',T71,']*',/,                           &
     &       'C[  SOLUTION:',T71,']*',/,'C[  ',9('='),T71,']*',/,       &
     &       'C[',T71,']*',/,                                           &
     &       'C[  A) TRANSFORMED COORDINATES:',T71,']*',/,              &
     &       'C[',5X,24('-'),T71,']*',/,'C[',T71,']*')                  
 2250 FORMAT('C[  SHIFT IN X DIRECTION: PX = ',F10.6,T71,']*',/,        &
     &       'C[  SHIFT IN Y DIRECTION: PY = ',F10.6,T71,']*',/,        &
     &       'C[  ROTATION OF COORDINATE SYSTEM (IN RAD): PHID = ',     &
     &          D14.8,T71,']*',/,                                       &
     &       'C[',T71,']*',/,                                           &
     &       'C[  TRANSFORMED NODES:',T71,']*',/,                       &
     &       'C[',T71,']*',/,                                           &
     &       'C[   I',9X,'PHI(I)',18X,'R(I)',T71,']*',/,                &
     &       'C[  ',48('-'),T71,']*')                                   
 2300 FORMAT('C[  ',I2,2(3X,D20.14),T71,']*') 
 2400 FORMAT('C[',T71,']*',/,                                           &
     &       'C[',T71,']*',/,                                           &
     &       'C[  B) SPLINE COEFFICIENTS:',T71,']*',/,                  &
     &       'C[',5X,20('-'),T71,']*',/,                                &
     &       'C[',T71,']*',/,                                           &
     &       'C[   I',8X,'A(I)',11X,'B(I)',11X,'C(I)',11X,'D(I)',       &
     &          T71,']*',/,                                             &
     &       'C[  ',62('-'),T71,']*')                                   
 2500 FORMAT('C[  ',I2,4(3X,D12.6),T71,']*') 
 2700 FORMAT('C[',T71,']*',/,                                           &
     &       'C[',T71,']*',/,                                           &
     &       'C[  C) TABLE OF VALUES FOR SPLINE:',T71,']*',/,           &
     &       'C[  ',30('-'),T71,']*',/,                                 &
     &       'C[',T71,']*',/,                                           &
     &       'C[   I',11X,'XW(I)',19X,'YW(I)',T71,']*',/,               &
     &       'C[  ',51('-'),T71,']*')                                   
 2800 FORMAT('C[  ',I3,2(4X,D20.14),T71,']*') 
 2900 FORMAT(//,' TABLE OF VALUES COMPUTED.',/,                         &
     &          ' ALL RESULTS ON TAPE3.',//)                            
      END PROGRAM TEST                              
