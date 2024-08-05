      PROGRAM TEST
C
C     Test program for the SUBROUTINE AKIMA1        7.4.93
C                                                  (Andreas Sohn)
C
C*****************************************************************
C                                                                *
C     Test results;                                              *
C                                                                *
C[                                                              ]*
C[ TEST EXAMPLE:                                                ]*
C[ =============                                                ]*
C[                                                              ]*
C[ NODES WITHOUT ROUNDING OF CORNERS:                           ]*
C[                                                              ]*
C[  I      X(I)        Y(I)                                     ]*
C[ --------------------------                                   ]*
C[  0   -70.00000    30.00000                                   ]*
C[  1   -60.00000    40.00000                                   ]*
C[  2   -50.00000    60.00000                                   ]*
C[  3   -40.00000    70.00000                                   ]*
C[  4   -30.00000    70.00000                                   ]*
C[  5   -20.00000    70.00000                                   ]*
C[  6   -10.00000    70.00000                                   ]*
C[  7      .00000    50.00000                                   ]*
C[  8    10.00000    30.00000                                   ]*
C[  9    20.00000    10.00000                                   ]*
C[ 10    30.00000    30.00000                                   ]*
C[ 11    40.00000    50.00000                                   ]*
C[ 12    50.00000    70.00000                                   ]*
C[ 13    60.00000    70.00000                                   ]*
C[ 14    70.00000    60.00000                                   ]*
C[ 15    80.00000    35.00000                                   ]*
C[ 16    90.00000    10.00000                                   ]*
C[ 17   110.00000    30.00000                                   ]*
C[                                                              ]*
C[ AKIMA SPLINE COEFFICIENTS:                                   ]*
C[                                                              ]*
C[  I        FN(I)         B(I)          C(I)          D(I)     ]*
C[ ----------------------------------------------------------   ]*
C[  0     30.000000      1.000000       .000000       .000000   ]*
C[  1     40.000000      1.000000       .250000      -.015000   ]*
C[  2     60.000000      1.500000       .000000      -.005000   ]*
C[  3     70.000000       .000000       .000000       .000000   ]*
C[  4     70.000000       .000000       .000000       .000000   ]*
C[  5     70.000000       .000000       .000000       .000000   ]*
C[  6     70.000000     -2.000000       .000000       .000000   ]*
C[  7     50.000000     -2.000000       .000000       .000000   ]*
C[  8     30.000000     -2.000000       .000000       .000000   ]*
C[  9     10.000000      2.000000       .000000       .000000   ]*
C[ 10     30.000000      2.000000       .000000       .000000   ]*
C[ 11     50.000000      2.000000       .000000       .000000   ]*
C[ 12     70.000000      2.000000      -.342857       .014286   ]*
C[ 13     70.000000      -.571429       .064286      -.010714   ]*
C[ 14     60.000000     -2.500000       .000000       .000000   ]*
C[ 15     35.000000     -2.500000       .000000       .000000   ]*
C[ 16     10.000000      1.000000       .000000       .000000   ]*
C[                                                              ]*
C[                                                              ]*
C[ NODES AFTER ROUNDING OF CORNERS:                             ]*
C[                                                              ]*
C[  I      X(I)        Y(I)                                     ]*
C[ --------------------------                                   ]*
C[  0   -70.00000    30.00000                                   ]*
C[  1   -60.00000    40.00000                                   ]*
C[  2   -50.00000    60.00000                                   ]*
C[  3   -40.00000    70.00000                                   ]*
C[  4   -30.00000    70.00000                                   ]*
C[  5   -20.00000    70.00000                                   ]*
C[  6   -15.00000    70.00000                                   ]*
C[  7    -5.00000    60.00000                                   ]*
C[  8      .00000    50.00000                                   ]*
C[  9    10.00000    30.00000                                   ]*
C[ 10    15.00000    20.00000                                   ]*
C[ 11    25.00000    20.00000                                   ]*
C[ 12    30.00000    30.00000                                   ]*
C[ 13    40.00000    50.00000                                   ]*
C[ 14    50.00000    70.00000                                   ]*
C[ 15    60.00000    70.00000                                   ]*
C[ 16    70.00000    60.00000                                   ]*
C[ 17    80.00000    35.00000                                   ]*
C[ 18    85.00000    22.50000                                   ]*
C[ 19    95.00000    15.00000                                   ]*
C[ 20   110.00000    30.00000                                   ]*
C[                                                              ]*
C[ AKIMA SPLINE COEFFICIENTS:                                   ]*
C[                                                              ]*
C[  I        FN(I)         B(I)          C(I)          D(I)     ]*
C[ ----------------------------------------------------------   ]*
C[  0     30.000000      1.000000       .000000       .000000   ]*
C[  1     40.000000      1.000000       .250000      -.015000   ]*
C[  2     60.000000      1.500000       .000000      -.005000   ]*
C[  3     70.000000       .000000       .000000       .000000   ]*
C[  4     70.000000       .000000       .000000       .000000   ]*
C[  5     70.000000       .000000       .000000       .000000   ]*
C[  6     70.000000       .000000      -.100000       .000000   ]*
C[  7     60.000000     -2.000000       .000000       .000000   ]*
C[  8     50.000000     -2.000000       .000000       .000000   ]*
C[  9     30.000000     -2.000000       .000000       .000000   ]*
C[ 10     20.000000     -2.000000       .200000       .000000   ]*
C[ 11     20.000000      2.000000       .000000       .000000   ]*
C[ 12     30.000000      2.000000       .000000       .000000   ]*
C[ 13     50.000000      2.000000       .000000       .000000   ]*
C[ 14     70.000000      2.000000      -.342857       .014286   ]*
C[ 15     70.000000      -.571429       .064286      -.010714   ]*
C[ 16     60.000000     -2.500000       .000000       .000000   ]*
C[ 17     35.000000     -2.500000       .000000       .000000   ]*
C[ 18     22.500000     -2.500000       .175000       .000000   ]*
C[ 19     15.000000      1.000000       .000000       .000000   ]*
C                                                                *
C*****************************************************************
C
      PARAMETER (N=17, NK=26)
      DOUBLE PRECISION XN(0:NK), FN(0:NK), B(0:NK-1), C(0:NK-1),
     +                 D(0:NK-1), HELP(-2:NK+1, 1:4), BETA
C
C the declaration below must be 'uncommented', in case that output 
C in form of a table of values is desired 
C see below!
C
C     DOUBLE PRECISION X, XA, XE, XSCH, SPVAL
C
      INTEGER NI, IERR, ISWTCH
C
C     initialize data, change if desired
C
      DATA (XN(I),FN(I), I=0,N)
     F            / -70.0D0 ,  30.0D0 ,
     F              -60.0D0 ,  40.0D0 ,
     F              -50.0D0 ,  60.0D0 ,
     F              -40.0D0 ,  70.0D0 ,
     F              -30.0D0 ,  70.0D0 ,
     F              -20.0D0 ,  70.0D0 ,
     F              -10.0D0 ,  70.0D0 ,
     F                0.0D0 ,  50.0D0 ,
     F               10.0D0 ,  30.0D0 ,
     F               20.0D0 ,  10.0D0 ,
     F               30.0D0 ,  30.0D0 ,
     F               40.0D0 ,  50.0D0 ,
     F               50.0D0 ,  70.0D0 ,
     F               60.0D0 ,  70.0D0 ,
     F               70.0D0 ,  60.0D0 ,
     F               80.0D0 ,  35.0D0 ,
     F               90.0D0 ,  10.0D0 ,
     F              110.0D0 ,  30.0D0 /
C
      NI      = 17
      ISWTCH  = 1
C
C     output of nodes
C
      WRITE (*,2000)
      WRITE (*,2005)
      DO 10, I = 0,NI
         WRITE(*,2010) I,XN(I),FN(I)
   10 CONTINUE
C
C     Compute Akima coefficients with rounding
C
      BETA = -0.5D+00
      CALL AKIMA1 (NI,XN,FN,NK,BETA,B,C,D,ISWTCH,HELP,IERR)
C
C     output
C
      IF (IERR .EQ. 0) THEN
         WRITE(*,2020)
         DO 20,I = 0, NI-1
   20       WRITE(*,2030)I,FN(I),B(I),C(I),D(I)
C
C when table of values is desired, 'uncomment' these instructions
C see above! 
C
C        OPEN(1, FILE='ORUNDUNG')
C        REWIND(1)
C        XSCH = (XN(NI)-XN(0))/150.0D0
C        XA = XN(0)  - 0.5D0*XSCH
C        XE = XN(NI) + XSCH
C        DO 111 X=XA, XE+0.5D0*XSCH, XSCH
C           WRITE(1,'(1X,2(F10.5,2X))') X, SPVAL(X,NI,XN,FN,B,C,D)
C 111       CONTINUE
C        CLOSE(1)
C
      ELSE
         WRITE(*,2040)
         IF (IERR .EQ. -1) THEN
            WRITE(*,2050)
         ELSEIF (IERR .EQ. -2) THEN
            WRITE(*,2060)
         ELSEIF (IERR .EQ. -3) THEN
            WRITE(*,2070)
         ELSEIF (IERR .EQ. -4) THEN
            WRITE(*,2080)
         ENDIF
      ENDIF
C
C     compute coefficients with rounding of corners
C
      BETA = 0.5D+00
      CALL AKIMA1 (NI,XN,FN,NK,BETA,B,C,D,ISWTCH,HELP,IERR)
C
C     output with rounded corners
C
      WRITE (*,2090)
      WRITE (*,2005)
      DO 30, I = 0,NI
         WRITE(*,2010) I,XN(I),FN(I)
   30 CONTINUE
C
C    output
C
      IF (IERR .EQ. 0) THEN
         WRITE(*,2020)
         DO 40,I = 0, NI-1
   40       WRITE(*,2030)I,FN(I),B(I),C(I),D(I)
C
C in case of desired table of values 'uncomment' these lines:
C
C        OPEN(2, FILE='RUNDUNG')
C        REWIND(2)
C        XSCH = (XN(NI)-XN(0))/150.0D0
C        XA = XN(0)  - 0.5D0*XSCH
C        XE = XN(NI) + XSCH
C        DO 222 X=XA, XE+0.5D0*XSCH, XSCH
C           WRITE(2,'(1X,2(F10.5,2X))') X, SPVAL(X,NI,XN,FN,B,C,D)
C 222       CONTINUE
C        CLOSE(2)
C
      ELSE
         WRITE(*,2040)
         IF (IERR .EQ. -1) THEN
            WRITE(*,2050)
         ELSEIF (IERR .EQ. -2) THEN
            WRITE(*,2060)
         ELSEIF (IERR .EQ. -3) THEN
            WRITE(*,2070)
         ELSEIF (IERR .EQ. -4) THEN
            WRITE(*,2080)
         ENDIF
      ENDIF
 2000 FORMAT ( 1X,'C[', T66, ']*', /,
     +         1X,'C[ TEST EXAMPLE:', T66, ']*', /,
     +         1X,'C[ ', 13('='), T66, ']*', /,
     +         1X,'C[', T66, ']*', /,
     +         1X,'C[ NODES WITHOUT ROUNDING OF CORNERS:', T66, ']*')
 2005 FORMAT ( 1X,'C[', T66, ']*', /,
     +         1X,'C[  I', 5X, ' X(I)', 7X, ' Y(I)', T66, ']*', /,
     +         1X,'C[ ', 26('-'), T66, ']*')
 2010 FORMAT ( 1X,'C[ ', I2, 3X, F9.5, 3X, F9.5, T66, ']*')
 2020 FORMAT ( 1X,'C[', T66, ']*', /,
     +         1X,'C[ AKIMA SPLINE COEFFICIENTS:', T66, ']*', /,
     +         1X,'C[', T66, ']*', /,
     +         1X,'C[  I', 7X, ' FN(I)', 8X, ' B(I)', 9X, ' C(I)', 9X,
     +            ' D(I)', T66, ']*', /,
     +         1X,'C[ ', 58('-'), T66, ']*')
 2030 FORMAT ( 1X,'C[ ', I2, 4X, F10.6, 4X, F10.6, 4X, F10.6, 4X,
     +             F10.6, T66, ']*')
 2040 FORMAT ( 1X,'C[', T66, ']*', /,
     +         1X,1X,'C[ ERROR IN AKIMA!', T66, ']*')
 2050 FORMAT ( 1X,'C[ LESS THAN  4 NODES!', T66, ']*')
 2060 FORMAT ( 1X,'C[ NODES NOT MONOTONIC!',
     +            T66, ']*')
 2070 FORMAT ( 1X,'C[ FIRST AND LAST NODE DO NOT COINCIDE!',
     +            T66, ']*')
 2080 FORMAT ( 1X,'C[ NK TOO SMALL!', T66, ']*')
 2090 FORMAT ( 1X,'C[', T66, ']*', /,
     +         1X,'C[', T66, ']*', /,
     +         1X,'C[ NODES AFTER ROUNDING OF CORNERS:', T66, ']*')
      STOP
      END
