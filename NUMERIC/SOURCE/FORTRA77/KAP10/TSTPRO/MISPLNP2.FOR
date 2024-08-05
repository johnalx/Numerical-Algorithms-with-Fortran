      PROGRAM TEST
C
C***********************************************************************
C                Test program for the subroutines                      *
C                    ISPLNP, ISPLNK, SPTAB                             *
C----------------------------------------------------------------------*
C     required subroutines    :  ISPLNP, ISPLNK, SPTAB,                *
C                                TRDIG,  TRDIGP, TRDIGS                *
C----------------------------------------------------------------------*
C                                                                      *
C     We test a non parametric interpolating cubic spline with not a   *
C     node condition.                                                  *
C                                                                      *
C[  GIVEN:   N + 1 POINTS X(I),Y(I), I=0(1)N, (N=10).                 ]*
C[  ======                                                            ]*
C[                                                                    ]*
C[   I        X(I)         Y(I)                                       ]*
C[  ----------------------------                                      ]*
C[   0        .0000        .0000                                      ]*
C[   1       1.5000       2.0000                                      ]*
C[   2       2.5000        .0000                                      ]*
C[   3       5.5000      -4.0000                                      ]*
C[   4       8.0000      -2.0000                                      ]*
C[   5      10.0000        .0000                                      ]*
C[   6      12.0000       2.0000                                      ]*
C[   7      14.5000       4.0000                                      ]*
C[   8      17.5000        .0000                                      ]*
C[   9      18.5000      -2.0000                                      ]*
C[  10      20.0000        .0000                                      ]*
C[                                                                    ]*
C[                                                                    ]*
C[  TO FIND:   A) COEFFICIENTS OF A NON PARAMETRIC INTERPOLATING      ]*
C[  ========      CUBIC SPLINE FUNCTION  S(X) FOR THE NODES X(I),Y(I) ]*
C[                WITH NOT-A-NODE CONDITION.                          ]*
C[                                                                    ]*
C[             B) TABLE OF VALUES FOR S(X)                            ]*
C[                                                                    ]*
C[                                                                    ]*
C[  SOLUTION:                                                         ]*
C[  =========                                                         ]*
C[                                                                    ]*
C[  A) SPLINE COEFFICIENTS:                                           ]*
C[     --------------------                                           ]*
C[                                                                    ]*
C[   I        A(I)           B(I)           C(I)           D(I)       ]*
C[  --------------------------------------------------------------    ]*
C[   0    .000000D+00    .503922D+01   -.315294D+01    .454902D+00    ]*
C[   1    .200000D+01   -.134902D+01   -.110588D+01    .454902D+00    ]*
C[   2    .000000D+00   -.238039D+01    .258823D+00    .300655D-01    ]*
C[   3   -.400000D+01   -.156845D-01    .529413D+00   -.812556D-01    ]*
C[   4   -.200000D+01    .110784D+01   -.800044D-01    .130430D-01    ]*
C[   5    .000000D+00    .944335D+00   -.174634D-02    .147893D-01    ]*
C[   6    .200000D+01    .111482D+01    .869897D-01   -.851674D-01    ]*
C[   7    .400000D+01   -.471186D-01   -.551766D+00    .410093D-01    ]*
C[   8    .000000D+00   -.225046D+01   -.182683D+00    .433147D+00    ]*
C[   9   -.200000D+01   -.131639D+01    .111676D+01    .433147D+00    ]*
C                                                                      *
C     The results (coefficients and table of values) are sent to TAPE3.*
C                                                                      *
C***********************************************************************
C                                                                      *
C     Author      :  GÅnter Palm                                       *
C     Date        :  5.20.1988                                         *
C     Source code :  FORTRAN 77                                        *
C                                                                      *
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (N=10,MT=100)
C
      DIMENSION X(0:N),Y(0:N)
      DIMENSION B(0:N),C(0:N),D(0:N)
      DIMENSION XW(0:N+MT+2),YW(0:N+MT+2)
      DIMENSION HILF(5*N+1)
C
      DATA X / 0.0D0,  1.5D0,  2.5D0,  5.5D0, 8.0D0, 10.0D0, 12.0D0,
     +        14.5D0, 17.5D0, 18.5D0, 20.0D0/
      DATA Y / 0.0D0, 2.0D0,  0.0D0, -4.0D0, -2.0D0, 0.0D0, 2.0D0,
     +         4.0D0, 0.0D0, -2.0D0,  0.0D0/
C
      IR = 5
      CALL ISPLNP (N,X,Y,IR,0.0D0,0.0D0,B,C,D,HILF,IFEHL)
      IF (IFEHL .NE. 0) THEN
         WRITE (*,1100) IFEHL
         STOP
      ENDIF
      OPEN (3,FILE='TAPE3')
      REWIND (3)
      WRITE (3,2000) N
      DO 100 I = 0,N,1
         WRITE (3,2100) I,X(I),Y(I)
  100 CONTINUE
      WRITE (3,2200)
      DO 200 I = 0,N-1,1
         WRITE (3,2300) I,Y(I),B(I),C(I),D(I)
  200 CONTINUE
C
      WRITE (*,1150) X(0),X(N)
   10 WRITE (*,1151)
      READ (*,*,END=999) MARGIN
      IF (MARGIN .EQ. 1) THEN
         GOTO 999
      ELSEIF (MARGIN .EQ. 2) THEN
         XANF = X(0)
         XEND = X(N)
      ELSEIF (MARGIN .EQ. 3) THEN
         WRITE (*,1200)
         READ (*,*,END=999) XANF
         WRITE (*,1300)
         READ (*,*,END=999) XEND
      ELSE
         GOTO 10
      ENDIF
      CALL SPTAB (N,MT,XANF,XEND,X,Y,B,C,D,NT,XW,YW,IFEHL)
      IF (IFEHL .NE. 0) THEN
         WRITE (*,1400)
         STOP
      ENDIF
      WRITE (3,2400)
      DO 300 I = 0,NT,1
         WRITE (3,2500) I,XW(I),YW(I)
  300 CONTINUE
      WRITE (*,2600)
  999 STOP
 1100 FORMAT(' ERROR IN COMPUTING COEFFICIENTS.',/,
     +       ' ERROR PARAMETER =',I4)
 1150 FORMAT(///,' NO ERROR WITH COEFFICIENTS.',/,
     +' RESULTS ON TAPE3.',//,' TABLE OF VALUES:',/,
     +' EXIT 1: NO TABLE OF VALUES (PROGRAM STOP)',/,
     +' EXIT 2: END POINTS SET BY PROGRAM',/,
     +'         LEFT  ENDPOINT = X(0) = ',F6.2,/,
     +'         RIGHT ENDPOINT = X(N) = ',F6.2,/,
     +' EXIT 3: SET END POINTS')
 1151 FORMAT(' CHOOSE BY SELECTING 1, 2 OR 3: ')
 1200 FORMAT(' LEFT ENDPOINT ? ')
 1300 FORMAT(' RIGHT ENDPOINT ? ')
 1400 FORMAT(' ERRONEOUS INPUT.',/,
     +       ' XANF EXCEEDS XEND.')
 2000 FORMAT('C[  GIVEN:   N + 1 POINTS X(I),Y(I), I=0(1)N, ',
     +       '(N=',I2,').',T71,']*',/,
     +       'C[  ======',T71,']*',/,
     +       'C[',T71,']*',/,
     +       'C[',3X,'I',8X,'X(I)',9X,'Y(I)',T71,']*',/,
     +       'C[  ',28('-'),T71,']*')
 2100 FORMAT('C[  ',I2,4X,F9.4,4X,F9.4,T71,']*')
 2200 FORMAT('C[',T71,']*',/,'C[',T71,']*',/,
     +       'C[  TO FIND:   A) COEFFICIENTS OF A NON PARAMETRIC',
     +       ' INTERPOLATING',T71,']*',/,
     +       'C[  ========      CUBIC SPLINE FUNCTION  S(X) FOR THE',
     +       ' NODES X(I),Y(I)',T71,']*',/,
     +       'C[                WITH NOT-A-NODE CONDITION.',T71,']*',/,
     +       'C[',T71,']*',/,
     +       'C[',13X,'B) TABLE OF VALUES FOR S(X)',T71,']*',/,
     +       'C[',T71,']*',/,'C[',T71,']*',/,
     +       'C[  SOLUTION:',T71,']*',/,'C[  ',9('='),T71,']*',/,
     +       'C[',T71,']*',/,
     +       'C[  A) SPLINE COEFFICIENTS:',T71,']*',/,
     +       'C[',5X,20('-'),T71,']*',/,'C[',T71,']*',/,
     +       'C[   I',8X,'A(I)',11X,'B(I)',11X,'C(I)',11X,'D(I)',
     +          T71,']*',/,'C[  ',62('-'),T71,']*')
 2300 FORMAT('C[  ',I2,4(3X,D12.6),T71,']*')
 2400 FORMAT('C[',T71,']*',/,
     +       'C[',T71,']*',/,
     +       'C[  B) TABLE OF VALUES FOR SPLINE FUNCTION:',T71,']*',/,
     +       'C[  ',39('-'),T71,']*',/,
     +       'C[',T71,']*',/,
     +       'C[   I',11X,'XW(I)',19X,'YW(I)',T71,']*',/,
     +       'C[  ',51('-'),T71,']*')
 2500 FORMAT('C[  ',I3,2(4X,D20.14),T71,']*')
 2600 FORMAT(//,' TABLE OF VALUES DONE.',/,
     +          ' ALL RESULTS ON TAPE3.',//)
      END
