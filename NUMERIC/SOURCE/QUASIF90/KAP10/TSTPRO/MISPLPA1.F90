!                                                                       
!***********************************************************************
!                Test program FOR THE SubroutineS                      *
!                 ISPLPA, ISPL1D, ISPL2D, PSPTAB                       *
!----------------------------------------------------------------------*
!     REQUIRED SUBROUTINES    :  ISPLPA, ISPL1D, ISPL2D, PSPTAB,       *
!                                TRDSY,  TRDSYP, TRDSYS                *
!----------------------------------------------------------------------*
!                                                                      *
!     We compute a parametric interpolating cubic spline with given    *
!     end point derivatives.                                           *
!                                                                      *
![  GIVEN:   N + 1 POINTS X(I),Y(I), I=0(1)N, (N= 8).                 ]*
![  ======                                                            ]*
![                                                                    ]*
![   I        X(I)         Y(I)                                       ]*
![  ----------------------------                                      ]*
![   0       1.0000       1.0000                                      ]*
![   1       1.5000       2.0000                                      ]*
![   2       2.0000       2.5000                                      ]*
![   3       2.5000       2.0000                                      ]*
![   4       2.5000       1.0000                                      ]*
![   5       2.0000       1.5000                                      ]*
![   6       2.5000       2.0000                                      ]*
![   7       3.0000       3.0000                                      ]*
![   8       4.0000       3.0000                                      ]*
![                                                                    ]*
![                                                                    ]*
![  TO FIND:   A) TABLE OF VALUES FOR THE PARAMETERS T(I), I=0(1)N    ]*
![  ========                                                          ]*
![             B) COEFFICIENTS OF A PARAMETRIC INTERPOLATING CUBIC    ]*
![                SPLINE WITH COORDINATE FUNCTIONS SX AND SY FOR THE  ]*
![                NODES X(I), Y(I).                                   ]*
![                END POINT CONDITIONS FOR SX AND SY:                 ]*
![                1. DERIVATIVE WITH RESPECT TO PARAMETER GIVEN:      ]*
![                FOR SX: DERIVATIVE AT  X(0) =     .1000             ]*
![                        DITTO AT       X(N) =    1.0000             ]*
![                FOR SY: DERIVATIVE AT  X(0) =    1.0000             ]*
![                        DITTO AT       X(N) =   -1.0000             ]*
![                                                                    ]*
![             C) TABLE OF VALUES WITH SX(T), SY(T)                   ]*
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
![   1    .11180339887499D+01                                         ]*
![   2    .18251407699364D+01                                         ]*
![   3    .25322475511230D+01                                         ]*
![   4    .35322475511230D+01                                         ]*
![   5    .42393543323095D+01                                         ]*
![   6    .49464611134961D+01                                         ]*
![   7    .60644951022460D+01                                         ]*
![   8    .70644951022460D+01                                         ]*
![                                                                    ]*
![                                                                    ]*
![  B) SPLINE COEFFICIENTS:                                           ]*
![     --------------------                                           ]*
![                                                                    ]*
![   I       AX(I)          BX(I)          CX(I)          DX(I)       ]*
![  --------------------------------------------------------------    ]*
![   0    .100000D+01    .100000D+00    .432645D+00   -.109198D+00    ]*
![   1    .150000D+01    .657929D+00    .663821D-01    .447679D-02    ]*
![   2    .200000D+01    .758523D+00    .758788D-01   -.210141D+00    ]*
![   3    .250000D+01    .550620D+00   -.369897D+00   -.180723D+00    ]*
![   4    .250000D+01   -.731343D+00   -.912066D+00    .133833D+01    ]*
![   5    .200000D+01   -.137064D-01    .192696D+01   -.128350D+01    ]*
![   6    .250000D+01    .786169D+00   -.795763D+00    .440588D+00    ]*
![   7    .300000D+01    .658993D+00    .682014D+00   -.341007D+00    ]*
![                                                                    ]*
![                                                                    ]*
![   I       AY(I)          BY(I)          CY(I)          DY(I)       ]*
![  --------------------------------------------------------------    ]*
![   0    .100000D+01    .100000D+01   -.240006D+00    .130209D+00    ]*
![   1    .200000D+01    .951616D+00    .196730D+00   -.767237D+00    ]*
![   2    .250000D+01    .789791D-01   -.143082D+01    .451320D+00    ]*
![   3    .200000D+01   -.126753D+01   -.473430D+00    .740962D+00    ]*
![   4    .100000D+01    .849490D-02    .174946D+01   -.107688D+01    ]*
![   5    .150000D+01    .867277D+00   -.534955D+00    .436200D+00    ]*
![   6    .200000D+01    .765036D+00    .390365D+00   -.245640D+00    ]*
![   7    .300000D+01    .716768D+00   -.433537D+00   -.283232D+00    ]*
![                                                                    ]*
![                                                                    ]*
![  TO FIND:   A) TABLE OF VALUES FOR THE PARAMETERS T(I), I=0(1)N    ]*
![  ========                                                          ]*
![             B) COEFFICIENTS OF A PARAMETRIC INTERPOLATING CUBIC    ]*
![                SPLINE WITH COORDINATE FUNCTIONS SX AND SY FOR THE  ]*
![                NODES X(I), Y(I).                                   ]*
![                END POINT CONDITIONS FOR SX AND SY:                 ]*
![                2. DERIVATIVE WITH RESPECT TO PARAMETER GIVEN:      ]*
![                FOR SX: DERIVATIVE AT  X(0) =     .0000             ]*
![                        DITTO AT       X(N) =     .0000             ]*
![                FOR SY: DERIVATIVE AT  X(0) =     .0000             ]*
![                        DITTO AT       X(N) =     .0000             ]*
![                                                                    ]*
![             C) TABLE OF VALUES WITH SX(T), SY(T)                   ]*
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
![   1    .11180339887499D+01                                         ]*
![   2    .18251407699364D+01                                         ]*
![   3    .25322475511230D+01                                         ]*
![   4    .35322475511230D+01                                         ]*
![   5    .42393543323095D+01                                         ]*
![   6    .49464611134961D+01                                         ]*
![   7    .60644951022460D+01                                         ]*
![   8    .70644951022460D+01                                         ]*
![                                                                    ]*
![                                                                    ]*
![  B) SPLINE COEFFICIENTS:                                           ]*
![     --------------------                                           ]*
![                                                                    ]*
![   I       AX(I)          BX(I)          CX(I)          DX(I)       ]*
![  --------------------------------------------------------------    ]*
![   0    .100000D+01    .370426D+00    .000000D+00    .614298D-01    ]*
![   1    .150000D+01    .600788D+00    .206042D+00   -.787499D-01    ]*
![   2    .200000D+01    .774050D+00    .389881D-01   -.189025D+00    ]*
![   3    .250000D+01    .545651D+00   -.361994D+00   -.183656D+00    ]*
![   4    .250000D+01   -.729307D+00   -.912964D+00    .133553D+01    ]*
![   5    .200000D+01   -.171435D-01    .192011D+01   -.126695D+01    ]*
![   6    .250000D+01    .797881D+00   -.767496D+00    .405935D+00    ]*
![   7    .300000D+01    .603965D+00    .594052D+00   -.198017D+00    ]*
![                                                                    ]*
![                                                                    ]*
![   I       AY(I)          BY(I)          CY(I)          DY(I)       ]*
![  --------------------------------------------------------------    ]*
![   0    .100000D+01    .850027D+00    .000000D+00    .355201D-01    ]*
![   1    .200000D+01    .983228D+00    .119138D+00   -.720728D+00    ]*
![   2    .250000D+01    .706219D-01   -.140976D+01    .438241D+00    ]*
![   3    .200000D+01   -.126572D+01   -.480109D+00    .745824D+00    ]*
![   4    .100000D+01    .115390D-01    .175736D+01   -.109415D+01    ]*
![   5    .150000D+01    .855599D+00   -.563681D+00    .500181D+00    ]*
![   6    .200000D+01    .808705D+00    .497362D+00   -.376276D+00    ]*
![   7    .300000D+01    .509804D+00   -.764707D+00    .254902D+00    ]*
![                                                                    ]*
![                                                                    ]*
![  TO FIND:   A) TABLE OF VALUES FOR THE PARAMETERS T(I), I=0(1)N    ]*
![  ========                                                          ]*
![             B) COEFFICIENTS OF A PARAMETRIC INTERPOLATING CUBIC    ]*
![                SPLINE WITH COORDINATE FUNCTIONS SX AND SY FOR THE  ]*
![                NODES X(I), Y(I).                                   ]*
![                END POINT CONDITIONS FOR SX AND SY:                 ]*
![                1ST DERIVATIVE (DY/DX) GIVEN:                       ]*
![                AT X(0) :    1.5000                                 ]*
![                AT X(N) :    -.6000                                 ]*
![                                                                    ]*
![             C) TABLE OF VALUES WITH SX(T), SY(T)                   ]*
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
![   1    .11180339887499D+01                                         ]*
![   2    .18251407699364D+01                                         ]*
![   3    .25322475511230D+01                                         ]*
![   4    .35322475511230D+01                                         ]*
![   5    .42393543323095D+01                                         ]*
![   6    .49464611134961D+01                                         ]*
![   7    .60644951022460D+01                                         ]*
![   8    .70644951022460D+01                                         ]*
![                                                                    ]*
![                                                                    ]*
![  B) SPLINE COEFFICIENTS:                                           ]*
![     --------------------                                           ]*
![                                                                    ]*
![   I       AX(I)          BX(I)          CX(I)          DX(I)       ]*
![  --------------------------------------------------------------    ]*
![   0    .100000D+01    .554700D+00   -.294852D+00    .177734D+00    ]*
![   1    .150000D+01    .561894D+00    .301286D+00   -.135658D+00    ]*
![   2    .200000D+01    .784490D+00    .135115D-01   -.173874D+00    ]*
![   3    .250000D+01    .542786D+00   -.355332D+00   -.187454D+00    ]*
![   4    .250000D+01   -.730241D+00   -.917695D+00    .134408D+01    ]*
![   5    .200000D+01   -.119308D-01    .193354D+01   -.129636D+01    ]*
![   6    .250000D+01    .777964D+00   -.816459D+00    .465663D+00    ]*
![   7    .300000D+01    .698542D+00    .745422D+00   -.443965D+00    ]*
![                                                                    ]*
![                                                                    ]*
![   I       AY(I)          BY(I)          CY(I)          DY(I)       ]*
![  --------------------------------------------------------------    ]*
![   0    .100000D+01    .832050D+00    .287402D-01    .241955D-01    ]*
![   1    .200000D+01    .987048D+00    .109894D+00   -.715297D+00    ]*
      PROGRAM TEST 
![   2    .250000D+01    .695164D-01   -.140748D+01    .437232D+00    ]*
![   3    .200000D+01   -.126511D+01   -.479972D+00    .745086D+00    ]*
![   4    .100000D+01    .102007D-01    .175529D+01   -.108854D+01    ]*
![   5    .150000D+01    .859743D+00   -.553852D+00    .477992D+00    ]*
![   6    .200000D+01    .793466D+00    .460122D+00   -.330777D+00    ]*
![   7    .300000D+01    .581917D+00   -.649338D+00    .674210D-01    ]*
!                                                                      *
!     The parameter values T(I), I = 0, .., N, are computed in ISPLPA  *
!     and are stored together with the spline coefficients and a table *
!     of values for XW and YW on TAPE3.                                *
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
      PARAMETER (N = 8, MT = 100) 
!                                                                       
      DIMENSION X (0:N), Y (0:N), T (0:N) 
      DIMENSION BX (0:N), CX (0:N), DX (0:N) 
      DIMENSION BY (0:N), CY (0:N), DY (0:N) 
      DIMENSION XW (0:N + MT + 2), YW (0:N + MT + 2) 
      DIMENSION ALPHA (1:2, 1:3), BETA (1:2, 1:3) 
      DIMENSION HILF (5 * N + 1) 
!                                                                       
      DATA X / 1.0D0, 1.5D0, 2.0D0, 2.5D0, 2.5D0, 2.0D0, 2.5D0, 3.0D0,  &
      4.0D0 /                                                           
      DATA Y / 1.0D0, 2.0D0, 2.5D0, 2.0D0, 1.0D0, 1.5D0, 2.0D0, 3.0D0,  &
      3.0D0 /                                                           
      DATA ALPHA (1, 1) / 0.1D0 /, ALPHA (2, 1) / 1.0D0 / 
      DATA BETA (1, 1) / 1.0D0 /, BETA (2, 1) / - 1.0D0 / 
      DATA ALPHA (1, 2) / 0.0D0 /, ALPHA (2, 2) / 0.0D0 / 
      DATA BETA (1, 2) / 0.0D0 /, BETA (2, 2) / 0.0D0 / 
      DATA ALPHA (1, 3) / 1.5D0 /, BETA (1, 3) / - 0.6D0 / 
      DATA JT / 1 / 
!                                                                       
   10 WRITE ( *, 1000) 
      READ ( *, *, END = 999) IAB 
      IF (IAB.LT.1.OR.IAB.GT.3) GOTO 10 
      CALL ISPLPA (N, X, Y, T, JT, IAB, ALPHA (1, IAB), BETA (1, IAB),  &
      BX, CX, DX, BY, CY, DY, HILF, IFEHL)                              
      IF (IFEHL.NE.0) THEN 
         WRITE ( *, 1100) IFEHL 
         STOP 
      ENDIF 
      OPEN (3, FILE = 'TAPE3') 
      REWIND (3) 
      WRITE (3, 2000) N 
      DO 100 I = 0, N, 1 
         WRITE (3, 2100) I, X (I), Y (I) 
  100 END DO 
      WRITE (3, 2200) 
      IF (IAB.LE.2) THEN 
         WRITE (3, 2250) IAB, ALPHA (1, IAB), BETA (1, IAB), ALPHA (2,  &
         IAB), BETA (2, IAB)                                            
      ELSE 
         WRITE (3, 2260) ALPHA (1, IAB), BETA (1, IAB) 
      ENDIF 
      WRITE (3, 2270) 
      DO 150 I = 0, N, 1 
         WRITE (3, 2275) I, T (I) 
  150 END DO 
      WRITE (3, 2280) 
      DO 200 I = 0, N - 1, 1 
         WRITE (3, 2300) I, X (I), BX (I), CX (I), DX (I) 
  200 END DO 
      WRITE (3, 2290) 
      DO 250 I = 0, N - 1, 1 
         WRITE (3, 2300) I, Y (I), BY (I), CY (I), DY (I) 
  250 END DO 
!                                                                       
      WRITE ( *, 1150) T (0), T (N) 
   20 WRITE ( *, 1151) 
      READ ( *, *, END = 999) MARGIN 
      IF (MARGIN.EQ.1) THEN 
         GOTO 999 
      ELSEIF (MARGIN.EQ.2) THEN 
         TANF = T (0) 
         TEND = T (N) 
      ELSEIF (MARGIN.EQ.3) THEN 
         WRITE ( *, 1200) 
         READ ( *, *, END = 999) TANF 
         WRITE ( *, 1300) 
         READ ( *, *, END = 999) TEND 
      ELSE 
         GOTO 20 
      ENDIF 
      CALL PSPTAB (N, MT, TANF, TEND, T, X, BX, CX, DX, Y, BY, CY, DY,  &
      NT, XW, YW, IFEHL)                                                
      IF (IFEHL.NE.0) THEN 
         WRITE ( *, 1400) 
         STOP 
      ENDIF 
      WRITE (3, 2400) 
      DO 300 I = 0, NT, 1 
         WRITE (3, 2500) I, XW (I), YW (I) 
  300 END DO 
      WRITE ( *, 2600) 
  999 STOP 
 1000 FORMAT(' GIVEN END POINT CONDITION:',/,                           &
     & ' 1 - 1. END POINT DERIVATIVE WITH RESPECT TO PARAMETER GIVEN',/,&
     & ' 2 - 2. END POINT DERIVATIVE WITH RESPECT TO PARAMETER GIVEN',/,&
     & ' 3 - 1. END POINT DERIVATIVE (DY/DX) GIVEN ? ')                 
 1100 FORMAT(' ERROR WHEN COMPUTING COEFFICIENTS.',/,                   &
     &       ' ERROR PARAMETER =',I4)                                   
 1150 FORMAT(///,' COEFFICIENTS FOUND.',/,                              &
     &' RESULTS ON TAPE3.',//,' TABLE OF VALUES:',/,                    &
     &' EXIT 1: NO TABLE OF VALUES (PROGRAM STOP)',/,                   &
     &' EXIT 2: PARAMETER END POINTS ARE SET BY PROGRAM',/,             &
     &'         STARTING PARAMETER VALUE = T(0) = ',F6.2,/,             &
     &'         FINAL PARAMETER VALUE    = T(N) = ',F6.2,/,             &
     &' EXIT 3: SET PARAMETER END POINTS')                              
 1151 FORMAT(' SELECT BY ENTERING  1, 2 OR 3: ') 
 1200 FORMAT(' PARAMETER START (TANF) ? ') 
 1300 FORMAT(' FINAL PARAMETER VALUE (TEND) ? ') 
 1400 FORMAT(' ERRONEOUS INPUT.',/,                                     &
     &       ' TANF EXCEEDS TEND.')                                     
 2000 FORMAT('C[  GIVEN:   N + 1 POINTS X(I),Y(I), I=0(1)N, ',          &
     &          '(N=',I2,').',T71,']*',/,                               &
     &       'C[  ======',T71,']*',/,                                   &
     &       'C[',T71,']*',/,                                           &
     &       'C[',3X,'I',8X,'X(I)',9X,'Y(I)',T71,']*',/,                &
     &       'C[  ',28('-'),T71,']*')                                   
 2100 FORMAT('C[  ',I2,4X,F9.4,4X,F9.4,T71,']*') 
 2200 FORMAT('C[',T71,']*',/,'C[',T71,']*',/,                           &
     &       'C[  TO FIND:   A) TABLE OF VALUES FOR THE PARAMETERS',    &
     &       ' T(I), I=0(1)N',T71,']*',/,                               &
     &       'C[  ========',T71,']*',/,                                 &
     &       'C[             B) COEFFICIENTS OF A PARAMETRIC',          &
     &       ' INTERPOLATING CUBIC',T71,']*',/,                         &
     &       'C[                SPLINE WITH COORDINATE FUNCTIONS SX',   &
     &       ' AND SY FOR THE',T71,']*',/,                              &
     &       'C[                NODES X(I), Y(I).',T71,']*',/,          &
     &       'C[                END POINT CONDITIONS FOR SX AND SY:',   &
     &       T71,']*')                                                  
 2250 FORMAT('C[',16X,I1,'. DERIVATIVE WITH RESPECT TO',                &
     &       ' PARAMETER GIVEN:',T71,']*',/,                            &
     &       'C[                FOR SX: DERIVATIVE AT  X(0) = ',        &
     &        F9.4,T71,']*',/,                                          &
     &       'C[                        DITTO AT       X(N) = ',        &
     &        F9.4,T71,']*',/,                                          &
     &       'C[                FOR SY: DERIVATIVE AT  X(0) = ',        &
     &        F9.4,T71,']*',/,                                          &
     &       'C[                        DITTO AT       X(N) = ',        &
     &        F9.4,T71,']*')                                            
 2260 FORMAT('C[',16X,'1ST DERIVATIVE (DY/DX) GIVEN:',T71,']*',/,       &
     &       'C[',16X,'AT X(0) : ',F9.4,T71,']*',/,                     &
     &       'C[',16X,'AT X(N) : ',F9.4,T71,']*')                       
 2270 FORMAT('C[',T71,']*',/,                                           &
     &       'C[',13X,'C) TABLE OF VALUES WITH SX(T), ',                &
     &          'SY(T)',T71,']*',/,                                     &
     &       'C[',T71,']*',/,                                           &
     &       'C[',T71,']*',/,                                           &
     &       'C[  SOLUTION:',T71,']*',/,                                &
     &       'C[  ',9('='),T71,']*',/,                                  &
     &       'C[',T71,']*',/,                                           &
     &       'C[  A) PARAMETER VALUES T:',T71,']*',/,                   &
     &       'C[',5X,19('-'),T71,']*',/,                                &
     &       'C[',T71,']*',/,                                           &
     &       'C[   I',11X,'T(I)',T71,']*',/,                            &
     &       'C[  ',25('-'),T71,']*')                                   
 2275 FORMAT('C[  ',I2,3X,D20.14,T71,']*') 
 2280 FORMAT('C[',T71,']*',/,                                           &
     &       'C[',T71,']*',/,                                           &
     &       'C[  B) SPLINE COEFFICIENTS:',T71,']*',/,                  &
     &       'C[',5X,20('-'),T71,']*',/,                                &
     &       'C[',T71,']*',/,                                           &
     &       'C[   I',7X,'AX(I)',10X,'BX(I)',10X,'CX(I)',10X,'DX(I)',   &
     &          T71,']*',/,                                             &
     &       'C[  ',62('-'),T71,']*')                                   
 2290 FORMAT('C[',T71,']*',/,                                           &
     &       'C[',T71,']*',/,                                           &
     &       'C[   I',7X,'AY(I)',10X,'BY(I)',10X,'CY(I)',10X,'DY(I)',   &
     &          T71,']*',/,                                             &
     &       'C[  ',62('-'),T71,']*')                                   
 2300 FORMAT('C[  ',I2,4(3X,D12.6),T71,']*') 
 2400 FORMAT('C[',T71,']*',/,                                           &
     &       'C[',T71,']*',/,                                           &
     &       'C[  C) TABLE OF VALUES:',T71,']*',/,                      &
     &       'C[',5X,16('-'),T71,']*',/,                                &
     &       'C[',T71,']*',/,                                           &
     &       'C[   I',12X,'SX(T(I))',16X,'SY(T(I))',T71,']*',/,         &
     &       'C[  ',51('-'),T71,']*')                                   
 2500 FORMAT('C[  ',I3,2(4X,D20.14),T71,']*') 
 2600 FORMAT(//,' TABLE OF VALUES DONE.',/,                             &
     &          ' ALL RESULTS ON TAPE3.',//)                            
      END PROGRAM TEST                              
