      PROGRAM TEST 
!                                                                       
!***********************************************************************
!                Test program for the subroutines                      *
!                    ISPLNP, ISPLPE, SPTAB                             *
!----------------------------------------------------------------------*
!     required subroutines    :  ISPLNP, ISPLPE, SPTAB,                *
!                                CYTSY,  CYTSYP, CYTSYS                *
!----------------------------------------------------------------------*
!                                                                      *
!     We compute a non parametric interpolating cubic spline.          *
!                                                                      *
![  GIVEN:   N + 1 POINTS X(I),Y(I), I=0(1)N, (N=10).                 ]*
![  ======                                                            ]*
![                                                                    ]*
![   I        X(I)         Y(I)                                       ]*
![  ----------------------------                                      ]*
![   0        .0000        .0000                                      ]*
![   1       1.5000       2.0000                                      ]*
![   2       2.5000        .0000                                      ]*
![   3       5.5000      -4.0000                                      ]*
![   4       8.0000      -2.0000                                      ]*
![   5      10.0000        .0000                                      ]*
![   6      12.0000       2.0000                                      ]*
![   7      14.5000       4.0000                                      ]*
![   8      17.5000        .0000                                      ]*
![   9      18.5000      -2.0000                                      ]*
![  10      20.0000        .0000                                      ]*
![                                                                    ]*
![                                                                    ]*
![  TO FIND:   A) COEFFICIENTS OF A NON PARAMETRIC INTERPOLATING      ]*
![  ========      CUBIC SPLINE FUNCTION S(X) FOR THE NODES X(I),Y(I). ]*
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
![   0    .000000D+00    .236488D+01    .000000D+00   -.458467D+00    ]*
![   1    .200000D+01   -.729769D+00   -.206310D+01    .792871D+00    ]*
![   2    .000000D+00   -.247736D+01    .315511D+00    .219438D-01    ]*
![   3   -.400000D+01    .818741D-02    .513005D+00   -.785119D-01    ]*
![   4   -.200000D+01    .110111D+01   -.758347D-01    .126391D-01    ]*
![   5    .000000D+00    .949444D+00    .830939D-17    .126391D-01    ]*
![   6    .200000D+01    .110111D+01    .758347D-01   -.785119D-01    ]*
![   7    .400000D+01    .818741D-02   -.513005D+00    .219438D-01    ]*
![   8    .000000D+00   -.247736D+01   -.315511D+00    .792871D+00    ]*
![   9   -.200000D+01   -.729769D+00    .206310D+01   -.458467D+00    ]*
!                                                                      *
!     The results (coefficients and table of values) are sent to TAPE3.*
!                                                                      *
!***********************************************************************
!                                                                      *
!     Author      :  GÅnter Palm                                       *
!     Date        :  5.20..1988                                        *
!     Source code :  FORTRAN 77                                        *
!                                                                      *
!***********************************************************************
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      PARAMETER (N = 10, MT = 100) 
!                                                                       
      DIMENSION X (0:N), Y (0:N) 
      DIMENSION B (0:N), C (0:N), D (0:N) 
      DIMENSION XW (0:N + MT + 2), YW (0:N + MT + 2) 
      DIMENSION HILF (5 * N + 1) 
!                                                                       
      DATA X / 0.0D0, 1.5D0, 2.5D0, 5.5D0, 8.0D0, 10.0D0, 12.0D0,       &
      14.5D0, 17.5D0, 18.5D0, 20.0D0 /                                  
      DATA Y / 0.0D0, 2.0D0, 0.0D0, - 4.0D0, - 2.0D0, 0.0D0, 2.0D0,     &
      4.0D0, 0.0D0, - 2.0D0, 0.0D0 /                                    
!                                                                       
      IR = 4 
      CALL ISPLNP (N, X, Y, IR, 0.0D0, 0.0D0, B, C, D, HILF, IFEHL) 
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
      DO 200 I = 0, N - 1, 1 
         WRITE (3, 2300) I, Y (I), B (I), C (I), D (I) 
  200 END DO 
!                                                                       
      WRITE ( *, 1150) X (0), X (N) 
   10 WRITE ( *, 1151) 
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
         GOTO 10 
      ENDIF 
      CALL SPTAB (N, MT, XANF, XEND, X, Y, B, C, D, NT, XW, YW, IFEHL) 
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
 1100 FORMAT(' ERROR WITH COMPUTING COEFFICIENTS.',/,                   &
     &       ' ERROR PARAMETER =',I4)                                   
 1150 FORMAT(///,' COEFFICIENTS COMPUTED.',/,                           &
     &' RESULTS ON TAPE3.',//,' TABLE OF VALUES:',/,                    &
     &' EXIT 1: NO TABLE OF VALUES (PROGRAM STOP)',/,                   &
     &' EXIT 2: END POINTS SET BY PROGRAM',/,                           &
     &'         LEFT END POINT = X(0) = ',F6.2,/,                       &
     &'         RIGHT ENDPOINT = X(N) = ',F6.2,/,                       &
     &' EXIT 3: ENTER END POINTS')                                      
 1151 FORMAT(' CHOOSE BY ENTERING 1, 2 OR 3: ') 
 1200 FORMAT(' LEFT END POINT ? ') 
 1300 FORMAT(' RIGHT END POINT ? ') 
 1400 FORMAT(' ERRONEOUS INPUT.',/,                                     &
     &       ' XANF EXCEEDS XEND.')                                     
 2000 FORMAT('C[  GIVEN:   N + 1 POINTS X(I),Y(I), I=0(1)N, ',          &
     &          '(N=',I2,').',T71,']*',/,                               &
     &       'C[  ======',T71,']*',/,                                   &
     &       'C[',T71,']*',/,                                           &
     &       'C[',3X,'I',8X,'X(I)',9X,'Y(I)',T71,']*',/,                &
     &       'C[  ',28('-'),T71,']*')                                   
 2100 FORMAT('C[  ',I2,4X,F9.4,4X,F9.4,T71,']*') 
 2200 FORMAT('C[',T71,']*',/,'C[',T71,']*',/,                           &
     &       'C[  TO FIND:   A) COEFFICIENTS OF A NON PARAMETRIC',      &
     &       ' INTERPOLATING',T71,']*',/,                               &
     &       'C[  ========      CUBIC SPLINE FUNCTION S(X) FOR THE',    &
     &       ' NODES X(I),Y(I).',T71,']*',/,                            &
     &       'C[',T71,']*',/,                                           &
     &       'C[',13X,'B) TABLE OF VALUES FOR S(X)',T71,']*',/,         &
     &       'C[',T71,']*',/,'C[',T71,']*',/,                           &
     &       'C[  SOLUTION:',T71,']*',/,'C[  ',9('='),T71,']*',/,       &
     &       'C[',T71,']*',/,                                           &
     &       'C[  A) SPLINE COEFFICIENTS:',T71,']*',/,                  &
     &       'C[',5X,20('-'),T71,']*',/,'C[',T71,']*',/,                &
     &       'C[   I',8X,'A(I)',11X,'B(I)',11X,'C(I)',11X,'D(I)',       &
     &          T71,']*',/,'C[  ',62('-'),T71,']*')                     
 2300 FORMAT('C[  ',I2,4(3X,D12.6),T71,']*') 
 2400 FORMAT('C[',T71,']*',/,                                           &
     &       'C[',T71,']*',/,                                           &
     &       'C[  B) TABLE OF VALUES FOR SPLINE FUNCTION:',T71,']*',/,  &
     &       'C[  ',39('-'),T71,']*',/,                                 &
     &       'C[',T71,']*',/,                                           &
     &       'C[   I',11X,'XW(I)',19X,'YW(I)',T71,']*',/,               &
     &       'C[  ',51('-'),T71,']*')                                   
 2500 FORMAT('C[  ',I3,2(4X,D20.14),T71,']*') 
 2600 FORMAT(//,' TABLE OF VALUES DONE.',/,                             &
     &          ' ALL RESULTS ON TAPE3.',//)                            
      END PROGRAM TEST                              
