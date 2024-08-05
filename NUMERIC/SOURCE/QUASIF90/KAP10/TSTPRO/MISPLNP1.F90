      PROGRAM TEST 
!                                                                       
!***********************************************************************
!                 Test program for the subroutines                     *
!               ISPLNP, ISPL1D, ISPL2D, ISPL3D, SPTAB                  *
!----------------------------------------------------------------------*
!    required subroutines    :  ISPLNP, ISPL1D, I SPL2D, ISPL3D,       *
!                                SPTAB,  TRDSY,  TRDSYP, TRDSYS        *
!----------------------------------------------------------------------*
!                                                                      *
!     We determine a non parametric interpolating spline of degree     *
!     three for given first, second and third derivatives at the end-  *
!     points.                                                          *
!                                                                      *
![  GIVEN:   N + 1 POINTS X(I),Y(I), I=0(1)N, (N= 8).                 ]*
![  ======                                                            ]*
![                                                                    ]*
![   I        X(I)         Y(I)                                       ]*
![  ----------------------------                                      ]*
![   0       1.0000       4.0000                                      ]*
![   1       2.0000       5.0000                                      ]*
![   2       3.0000       4.5000                                      ]*
![   3       3.5000       3.2000                                      ]*
![   4       4.8000       2.7000                                      ]*
![   5       5.7000       2.0000                                      ]*
![   6       7.0000       1.0000                                      ]*
![   7       8.5000       1.5000                                      ]*
![   8       9.2000       2.5000                                      ]*
![                                                                    ]*
![                                                                    ]*
![  TO FIND:   A) COEFFICIENTS OF A NON PARAMETRIC INTERPOLATING      ]*
![  ========      CUBIC SPLINE FUNCTION S(X) FOR THE NODES X(I),Y(I). ]*
![                ENDPOINT CONDITION FOR S(X):                        ]*
![                1. DERIVATIVE AT  X(0) =    1.0000                  ]*
![                1. DERIVATIVE AT  X(N) =    1.5000                  ]*
![                                                                    ]*
![             B) TABLE OF VALUES FOR S(X)                            ]*
![                                                                    ]*
![                                                                    ]*
![  SOLUTION:                                                         ]*
![  =========                                                         ]*
![                                                                    ]*
![  A) SPLINE COEFFICIENTS:                                           ]*
![     --------------------                                           ]*
![                                                                    ]*
![   I        A(I)           B(I)           C(I)           D(I)       ]*
![  --------------------------------------------------------------    ]*
![   0    .400000D+01    .100000D+01    .312014D+00   -.312014D+00    ]*
![   1    .500000D+01    .687986D+00   -.624027D+00   -.563959D+00    ]*
![   2    .450000D+01   -.225195D+01   -.231590D+01    .323959D+01    ]*
![   3    .320000D+01   -.213816D+01    .254348D+01   -.918924D+00    ]*
![   4    .270000D+01   -.184055D+00   -.104032D+01    .422926D+00    ]*
![   5    .200000D+01   -.102893D+01    .101575D+00    .755342D-01    ]*
![   6    .100000D+01   -.381878D+00    .396158D+00    .537665D-01    ]*
![   7    .150000D+01    .116952D+01    .638107D+00   -.382905D+00    ]*
![                                                                    ]*
![                                                                    ]*
![  TO FIND:   A) COEFFICIENTS OF A NON PARAMETRIC INTERPOLATING      ]*
![  ========      CUBIC SPLINE FUNCTION S(X) FOR THE NODES X(I),Y(I). ]*
![                ENDPOINT CONDITION FOR S(X):                        ]*
![                2. DERIVATIVE AT  X(0) =     .0000                  ]*
![                2. DERIVATIVE AT  X(N) =     .0000                  ]*
![                                                                    ]*
![             B) TABLE OF VALUES FOR S(X)                            ]*
![                                                                    ]*
![                                                                    ]*
![  SOLUTION:                                                         ]*
![  =========                                                         ]*
![                                                                    ]*
![  A) SPLINE COEFFICIENTS:                                           ]*
![     --------------------                                           ]*
![                                                                    ]*
![   I        A(I)           B(I)           C(I)           D(I)       ]*
![  --------------------------------------------------------------    ]*
![   0    .400000D+01    .117958D+01    .000000D+00   -.179578D+00    ]*
![   1    .500000D+01    .640843D+00   -.538735D+00   -.602108D+00    ]*
![   2    .450000D+01   -.224295D+01   -.234506D+01    .326192D+01    ]*
![   3    .320000D+01   -.214157D+01    .254782D+01   -.920245D+00    ]*
![   4    .270000D+01   -.182876D+00   -.104114D+01    .422371D+00    ]*
![   5    .200000D+01   -.103056D+01    .992662D-01    .782733D-01    ]*
![   6    .100000D+01   -.375621D+00    .404532D+00    .454028D-01    ]*
![   7    .150000D+01    .114444D+01    .608845D+00   -.289926D+00    ]*
![                                                                    ]*
![                                                                    ]*
![  TO FIND:   A) COEFFICIENTS OF A NON PARAMETRIC INTERPOLATING      ]*
![  ========      CUBIC SPLINE FUNCTION S(X) FOR THE NODES X(I),Y(I). ]*
![                ENDPOINT CONDITION FOR S(X):                        ]*
![                3. DERIVATIVE AT  X(0) =   -1.2000                  ]*
![                3. DERIVATIVE AT  X(N) =   -1.7000                  ]*
![                                                                    ]*
![             B) TABLE OF VALUES FOR S(X)                            ]*
![                                                                    ]*
![                                                                    ]*
![  SOLUTION:                                                         ]*
![  =========                                                         ]*
![                                                                    ]*
![  A) SPLINE COEFFICIENTS:                                           ]*
![     --------------------                                           ]*
![                                                                    ]*
![   I        A(I)           B(I)           C(I)           D(I)       ]*
![  --------------------------------------------------------------    ]*
![   0    .400000D+01    .115189D+01    .481108D-01   -.200000D+00    ]*
![   1    .500000D+01    .648111D+00   -.551889D+00   -.596222D+00    ]*
![   2    .450000D+01   -.224433D+01   -.234055D+01    .325844D+01    ]*
![   3    .320000D+01   -.214106D+01    .254710D+01   -.919994D+00    ]*
![   4    .270000D+01   -.182963D+00   -.104087D+01    .422188D+00    ]*
![   5    .200000D+01   -.103062D+01    .990322D-01    .784901D-01    ]*
![   6    .100000D+01   -.375192D+00    .405144D+00    .448047D-01    ]*
![   7    .150000D+01    .114267D+01    .606765D+00   -.283333D+00    ]*
!                                                                      *
!     The results (coefficients and table of values) are sent to TAPE3.*
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
      DIMENSION X (0:N), Y (0:N) 
      DIMENSION B (0:N), C (0:N), D (0:N) 
      DIMENSION XW (0:N + MT + 2), YW (0:N + MT + 2) 
      DIMENSION ALPHA (3), BETA (3) 
      DIMENSION HILF (5 * N + 1) 
!                                                                       
      DATA X / 1.0D0, 2.0D0, 3.0D0, 3.5D0, 4.8D0, 5.7D0, 7.0D0, 8.5D0,  &
      9.2D0 /                                                           
      DATA Y / 4.0D0, 5.0D0, 4.5D0, 3.2D0, 2.7D0, 2.0D0, 1.0D0, 1.5D0,  &
      2.5D0 /                                                           
      DATA ALPHA / 1.0D0, 0.0D0, - 1.2D0 / 
      DATA BETA / 1.5D0, 0.0D0, - 1.7D0 / 
!                                                                       
   10 WRITE ( *, 1000) 
      READ ( *, *, END = 999) IAB 
      IF (IAB.LT.1.OR.IAB.GT.3) GOTO 10 
      CALL ISPLNP (N, X, Y, IAB, ALPHA (IAB), BETA (IAB), B, C, D, HILF,&
      IFEHL)                                                            
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
      WRITE (3, 2200) IAB, ALPHA (IAB), IAB, BETA (IAB) 
      DO 200 I = 0, N - 1, 1 
         WRITE (3, 2300) I, Y (I), B (I), C (I), D (I) 
  200 END DO 
!                                                                       
      WRITE ( *, 1150) X (0), X (N) 
   20 WRITE ( *, 1151) 
      READ ( *, *, END = 999) MARGIN 
      IF (MARGIN.EQ.1) THEN 
         GOTO 999 
      ELSEIF (MARGIN.EQ.2) THEN 
         XANF = X (0) 
         XEND = X (N) 
      ELSEIF (MARGIN.EQ.3) THEN 
         WRITE ( *, 1200) 
         READ ( *, *, END = 999) XANF 
         WRITE ( *, 1300) 
         READ ( *, *, END = 999) XEND 
      ELSE 
         GOTO 20 
      ENDIF 
      CALL SPTAB (N, MT, XANF, XEND, X, Y, B, C, D, NT, XW, YW, IFEHL) 
      WRITE (3, 2400) 
      DO 300 I = 0, NT, 1 
         WRITE (3, 2500) I, XW (I), YW (I) 
  300 END DO 
      WRITE ( *, 2600) 
  999 STOP 
 1000 FORMAT(' GIVEN ENDPOINT DERIVATIVE (1ST, 2ND OR 3RD) ? ') 
 1100 FORMAT(' ERROR IN COMPUTING COEFFICIENTS.',/,                     &
     &       ' ERROR PARAMETER =',I4)                                   
 1150 FORMAT(///,' COEFFICIENTS COMPUTED ; NO ERRROR.',/,               &
     &' RESULTS ON  TAPE3.',//,' TABLE OF VALUES:',/,                   &
     &' EXIT 1: NO TABLE OF VALUES (PROGRAM STOP)',/,                   &
     &' EXIT 2: END POINTS OF TABLE OF VALUES SET IN PROGRAM',/,        &
     &'         LEFT END POINT = X(0) = ',F6.2,/,                       &
     &'         RIGHT END POINT = X(N) = ',F6.2,/,                      &
     &' EXIT 3: PUT IN END POINTS')                                     
 1151 FORMAT(' CHOOSE BY ENTERING  1, 2 OR 3: ') 
 1200 FORMAT(' LEFT END POINT ? ') 
 1300 FORMAT(' RIGHT END POINT ? ') 
 1400 FORMAT(' ERRONEOUS INPUT.',/,                                     &
     &       ' XANF EXCEEEDS XEND.')                                    
 2000 FORMAT('C[  GIVEN:   N + 1 POINTS X(I),Y(I), I=0(1)N, ',          &
     &       '(N=',I2,').',T71,']*',/,                                  &
     &       'C[  ',6('='),T71,']*',/,                                  &
     &       'C[',T71,']*',/,                                           &
     &       'C[',3X,'I',8X,'X(I)',9X,'Y(I)',T71,']*',/,                &
     &       'C[  ',28('-'),T71,']*')                                   
 2100 FORMAT('C[  ',I2,4X,F9.4,4X,F9.4,T71,']*') 
 2200 FORMAT('C[',T71,']*',/,'C[',T71,']*',/,                           &
     &       'C[  TO FIND:   A) COEFFICIENTS OF A NON PARAMETRIC',      &
     &       ' INTERPOLATING      ]*',/,                                &
     &       'C[  ========      CUBIC SPLINE FUNCTION S(X) FOR THE',    &
     &       ' NODES X(I),Y(I).',T71,']*',/,                            &
     &       'C[',16X,'ENDPOINT CONDITION FOR S(X):',T71,']*',/,        &
     &       'C[',16X,I1,'. DERIVATIVE AT  X(0) = ',F9.4,T71,']*',/,    &
     &       'C[',16X,I1,'. DERIVATIVE AT  X(N) = ',F9.4,T71,']*',/,    &
     &       'C[',T71,']*',/,                                           &
     &       'C[',13X,'B) TABLE OF VALUES FOR S(X)',T71,']*',/,         &
     &       'C[',T71,']*',/,'C[',T71,']*',/,                           &
     &       'C[  SOLUTION:',T71,']*',/,'C[  ',9('='),T71,']*',/,       &
     &       'C[',T71,']*',/,                                           &
     &       'C[  A) SPLINE COEFFICIENTS:',T71,']*',/,                  &
     &       'C[',5x,20('-'),T71,']*',/,'C[',T71,']*',/,                &
     &       'C[   I',8X,'A(I)',11X,'B(I)',11X,'C(I)',11X,'D(I)',       &
     &          T71,']*',/,'C[  ',62('-'),T71,']*')                     
 2300 FORMAT('C[  ',I2,4(3X,D12.6),T71,']*') 
 2400 FORMAT('C[',T71,']*',/,                                           &
     &       'C[',T71,']*',/,                                           &
     &       'C[  B) TABLE OF VALUES FOR SPLINE:',T71,']*',/,           &
     &       'C[  ',30('-'),T71,']*',/,                                 &
     &       'C[',T71,']*',/,                                           &
     &       'C[   I',11X,'XW(I)',19X,'YW(I)',T71,']*',/,               &
     &       'C[  ',51('-'),T71,']*')                                   
 2500 FORMAT('C[  ',I3,2(4X,D20.14),T71,']*') 
 2600 FORMAT(//,' TABLE OF VALUES COMPUTED.',/,                         &
     &          ' RESULTS ON TAPE3.',//)                                
      END PROGRAM TEST                              
