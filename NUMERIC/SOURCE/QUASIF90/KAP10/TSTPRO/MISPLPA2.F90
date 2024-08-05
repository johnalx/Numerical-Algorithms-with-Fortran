      PROGRAM TEST 
!                                                                       
!***********************************************************************
!                 Test program for the subroutines                     *
!                      ISPLPA, ISPLPE, PSPTAB                          *
!----------------------------------------------------------------------*
!     required subroutines    :  ISPLPA, ISPLPE, PSPTAB,               *
!                                CYTSYS, CYTSY,  CYTSYP                *
!----------------------------------------------------------------------*
!                                                                      *
!     We compute a parametric periodic interpolating cubic spline.     *
!                                                                      *
![  GIVEN:   N + 1 POINTS X(I),Y(I), I=0(1)N, (N=12).                 ]*
![  ======                                                            ]*
![                                                                    ]*
![   I        X(I)         Y(I)                                       ]*
![  ----------------------------                                      ]*
![   0      32.0000       1.5000                                      ]*
![   1      25.0000       2.5000                                      ]*
![   2      16.0000       3.2000                                      ]*
![   3       5.0000       3.0000                                      ]*
![   4      -6.0000       3.1000                                      ]*
![   5     -18.0000       1.8000                                      ]*
![   6     -16.0000       -.7000                                      ]*
![   7     -13.0000      -2.2000                                      ]*
![   8      -6.0000      -3.4000                                      ]*
![   9       4.0000      -3.1000                                      ]*
![  10      21.0000      -2.5000                                      ]*
![  11      33.0000      -1.1000                                      ]*
![  12      32.0000       1.5000                                      ]*
![                                                                    ]*
![                                                                    ]*
![  TO FIND:   A) TABLE OF VALUES FOR THE PARAMETER T(I), I=0(1)N     ]*
![  ========                                                          ]*
![             B) COEFFICIENTS OF A PARAMETRIC, PERIODIC, CUBIC       ]*
![                INTERPOLATING SPLINE WITH COMPONENT FUNCTIONS       ]*
![                SX AND SY FOR THE NODES X(I), Y(I).                 ]*
![                                                                    ]*
![             C) TABLE OF VALUES FOR SX(T), SY(T)                    ]*
![                                                                    ]*
![                                                                    ]*
![  SOLUTION:                                                         ]*
![  =========                                                         ]*
![                                                                    ]*
![  A) PARAMETER VALUES T:                                            ]*
![     -------------------                                            ]*
![                                                                    ]*
![   I           T(I)                                                 ]*
![  -------------------------                                         ]*
![   0    .00000000000000D+00                                         ]*
![   1    .70710678118655D+01                                         ]*
![   2    .16098248988733D+02                                         ]*
![   3    .27100067020314D+02                                         ]*
![   4    .38100521556377D+02                                         ]*
![   5    .50170732822134D+02                                         ]*
![   6    .53372294940850D+02                                         ]*
![   7    .56726396907100D+02                                         ]*
![   8    .63828509268926D+02                                         ]*
![   9    .73833008256881D+02                                         ]*
![  10    .90843593196853D+02                                         ]*
![  11    .10292498384528D+03                                         ]*
![  12    .10571066150072D+03                                         ]*
![                                                                    ]*
![                                                                    ]*
![  B) SPLINE COEFFICIENTS:                                           ]*
![     --------------------                                           ]*
![                                                                    ]*
![   I       AX(I)          BX(I)          CX(I)          DX(I)       ]*
![  --------------------------------------------------------------    ]*
![   0    .320000D+02   -.645071D+00   -.831890D-01    .486712D-02    ]*
![   1    .250000D+02   -.109147D+01    .200583D-01   -.106254D-02    ]*
![   2    .160000D+02   -.989092D+00   -.871695D-02    .703566D-03    ]*
![   3    .500000D+01   -.925418D+00    .145046D-01   -.193453D-02    ]*
![   4   -.600000D+01   -.130860D+01   -.493376D-01    .624566D-02    ]*
![   5   -.180000D+02    .230161D+00    .176822D+00   -.167387D-01    ]*
![   6   -.160000D+02    .847657D+00    .160517D-01   -.628348D-03    ]*
![   7   -.130000D+02    .934128D+00    .972905D-02   -.348986D-03    ]*
![   8   -.600000D+01    .101951D+01    .229344D-02   -.428692D-03    ]*
![   9    .400000D+01    .936679D+00   -.105731D-01    .838239D-03    ]*
![  10    .210000D+02    .130463D+01    .322037D-01   -.479879D-02    ]*
![  11    .330000D+02   -.185340D-01   -.141724D+00    .700433D-02    ]*
![                                                                    ]*
![                                                                    ]*
![   I       AY(I)          BY(I)          CY(I)          DY(I)       ]*
![  --------------------------------------------------------------    ]*
![   0    .150000D+01    .758479D+00   -.146000D+00    .830643D-02    ]*
![   1    .250000D+01   -.603148D-01    .302056D-01   -.165435D-02    ]*
![   2    .320000D+01    .805893D-01   -.145967D-01    .510756D-03    ]*
![   3    .300000D+01   -.551250D-01    .226105D-02    .325121D-03    ]*
![   4    .310000D+01    .112649D+00    .129905D-01   -.258872D-02    ]*
![   5    .180000D+01   -.705206D+00   -.807486D-01    .178399D-01    ]*
![   6   -.700000D+00   -.673674D+00    .905976D-01   -.688120D-02    ]*
![   7   -.220000D+01   -.298167D+00    .213568D-01   -.445576D-03    ]*
![   8   -.340000D+01   -.622353D-01    .118632D-01   -.264396D-03    ]*
![   9   -.310000D+01    .957449D-01    .392773D-02   -.439887D-03    ]*
![  10   -.250000D+01   -.152486D+00   -.185205D-01    .337161D-02    ]*
![  11   -.110000D+01    .876368D+00    .103681D+00   -.298768D-01    ]*
!                                                                      *
!     The parameter values T(I), I=0, .., N, are computed by the sub-  *
!     routine ISPLPA. They are put onto TAPE3 together with the spline *
!     coefficients and a table of values for the points XW, YW on the  *
!     curve.                                                           *
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
      DIMENSION X (0:N), Y (0:N), T (0:N) 
      DIMENSION BX (0:N), CX (0:N), DX (0:N) 
      DIMENSION BY (0:N), CY (0:N), DY (0:N) 
      DIMENSION XW (0:N + MT + 2), YW (0:N + MT + 2) 
      DIMENSION ALPHA (2), BETA (2) 
      DIMENSION HILF (5 * N + 1) 
!                                                                       
      DATA X / 32.0D0, 25.0D0, 16.0D0, 5.0D0, - 6.0D0, - 18.0D0,        &
      - 16.0D0, - 13.0D0, - 6.0D0, 4.0D0, 21.0D0, 33.0D0, 32.0D0 /      
      DATA Y / 1.5D0, 2.5D0, 3.2D0, 3.0D0, 3.1D0, 1.8D0, - 0.7D0,       &
      - 2.2D0, - 3.4D0, - 3.1D0, - 2.5D0, - 1.1D0, 1.5D0 /              
      DATA JT / 1 / 
      DATA IR / 4 / 
!                                                                       
      CALL ISPLPA (N, X, Y, T, JT, IR, ALPHA, BETA, BX, CX, DX, BY, CY, &
      DY, HILF, IFEHL)                                                  
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
      DO 150 I = 0, N, 1 
         WRITE (3, 2300) I, T (I) 
  150 END DO 
      WRITE (3, 2400) 
      DO 200 I = 0, N - 1, 1 
         WRITE (3, 2500) I, X (I), BX (I), CX (I), DX (I) 
  200 END DO 
      WRITE (3, 2600) 
      DO 250 I = 0, N - 1, 1 
         WRITE (3, 2500) I, Y (I), BY (I), CY (I), DY (I) 
  250 END DO 
!                                                                       
      WRITE ( *, 1050) T (0), T (N) 
   10 WRITE ( *, 1051) 
      READ ( *, *, END = 999) MARGIN 
      IF (MARGIN.EQ.1) THEN 
         GOTO 999 
      ELSEIF (MARGIN.EQ.2) THEN 
         TANF = T (0) 
         TEND = T (N) 
      ELSEIF (MARGIN.EQ.3) THEN 
         WRITE ( *, 1100) 
         READ ( *, *, END = 999) TANF 
         WRITE ( *, 1200) 
         READ ( *, *, END = 999) TEND 
      ELSE 
         GOTO 10 
      ENDIF 
      CALL PSPTAB (N, MT, TANF, TEND, T, X, BX, CX, DX, Y, BY, CY, DY,  &
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
 1000 FORMAT(' ERROR WHEN COMPUTING COEFFICIENTS.',/,                   &
     &       ' ERROR PARAMETER =',I4)                                   
 1050 FORMAT(///,' COEFFICIENTS COMPUTED.',/,                           &
     &' RESULTS ON TAPE3.',//,' TABLE OF VALUES:',/,                    &
     &' EXIT 1: NO TABLE OF VALUES (PROGRAM STOP)',/,                   &
     &' EXIT 2: END POINTS ARE SET BY PROGRAMM ',/,                     &
     &'         FIRST PARAMETER POINT = T(0) = ',F6.2,/,                &
     &'         LAST PARAMETER POINT  = T(N) = ',F6.2,/,                &
     &' EXIT 3: SET UP PARAMETER END POINTS')                           
 1051 FORMAT(' CHOOSE BY ENTERING  1, 2 OR 3: ') 
 1100 FORMAT(' FIRST PARAMETER VALUE (TANF) ? ') 
 1200 FORMAT(' LAST PARAMETER VALUE (TEND) ? ') 
 1300 FORMAT(' ERRONEOUS INPUT.',/,                                     &
     &       ' TANF EXCEEDS TEND.')                                     
 2000 FORMAT('C[  GIVEN:   N + 1 POINTS X(I),Y(I), I=0(1)N, ',          &
     &          '(N=',I2,').',T71,']*',/,                               &
     &       'C[  ======',T71,']*',/,                                   &
     &       'C[',T71,']*',/,                                           &
     &       'C[',3X,'I',8X,'X(I)',9X,'Y(I)',T71,']*',/,                &
     &       'C[  ',28('-'),T71,']*')                                   
 2100 FORMAT('C[  ',I2,4X,F9.4,4X,F9.4,T71,']*') 
 2200 FORMAT('C[',T71,']*',/,'C[',T71,']*',/,                           &
     &       'C[  TO FIND:   A) TABLE OF PARAMETER VALUES ',            &
     &       'T(I), ','I=0(1)N',T71,']*',/,                             &
     &       'C[  ========',T71,']*',/,                                 &
     &       'C[',13X,'B) COEFFICIENTS OF A PARAMETRIC, PERIODIC, ',    &
     &       'CUBIC ',T71,']*',/,                                       &
     &       'C[',16X,'INTERPOLATING SPLINE WITH COMPONENT FUNCTIONS',  &
     &          T71,']*',/,                                             &
     &       'C[',16X,'SX AND SY FOR THE NODES X(I), Y(I).',T71,']*',/, &
     &       'C[',T71,']*',/,                                           &
     &       'C[',13X,'C) TABLE OF VALUES WITH SX(T), SY(T) ',          &
     &          T71,']*',/,'C[',T71,']*',/,'C[',T71,']*',/,             &
     &       'C[  SOLUTION:',T71,']*',/,'C[  ',9('='),T71,']*',/,       &
     &       'C[',T71,']*',/,'C[  A) PARAMETER VALUES T:',T71,']*',/,   &
     &       'C[',5X,19('-'),T71,']*',/,'C[',T71,']*',/,                &
     &       'C[   I',11X,'T(I)',T71,']*',/,'C[  ',25('-'),T71,']*')    
 2300 FORMAT('C[  ',I2,3X,D20.14,T71,']*') 
 2400 FORMAT('C[',T71,']*',/,                                           &
     &       'C[',T71,']*',/,                                           &
     &       'C[  B) SPLINE COEFFICIENTS:',T71,']*',/,                  &
     &       'C[',5X,20('-'),T71,']*',/,                                &
     &       'C[',T71,']*',/,                                           &
     &       'C[   I',7X,'AX(I)',10X,'BX(I)',10X,'CX(I)',10X,'DX(I)',   &
     &          T71,']*',/,                                             &
     &       'C[  ',62('-'),T71,']*')                                   
 2500 FORMAT('C[  ',I2,4(3X,D12.6),T71,']*') 
 2600 FORMAT('C[',T71,']*',/,                                           &
     &       'C[',T71,']*',/,                                           &
     &       'C[   I',7X,'AY(I)',10X,'BY(I)',10X,'CY(I)',10X,'DY(I)',   &
     &          T71,']*',/,                                             &
     &       'C[  ',62('-'),T71,']*')                                   
 2700 FORMAT('C[',T71,']*',/,                                           &
     &       'C[',T71,']*',/,                                           &
     &       'C[  C) TABLE OF VALUES:',T71,']*',/,                      &
     &       'C[',5X,16('-'),T71,']*',/,                                &
     &       'C[',T71,']*',/,                                           &
     &       'C[   I',12X,'SX(T(I))',16X,'SY(T(I))',T71,']*',/,         &
     &       'C[  ',51('-'),T71,']*')                                   
 2800 FORMAT('C[  ',I3,2(4X,D20.14),T71,']*') 
 2900 FORMAT(//,' TABLE COMPUTED.',/,                                   &
     &          ' ALL RESULTS ON TAPE3.',//)                            
      END PROGRAM TEST                              
