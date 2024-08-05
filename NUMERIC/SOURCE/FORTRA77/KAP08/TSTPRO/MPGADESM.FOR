C                                                         4.26.87
      PROGRAM TEST
C
C*****************************************************************************
C                                                                            *
C  Testprogram for the  SUBROUTINE GADESM.                                   *
C                                                                            *
C  Test example: Numerische Mathematik fÅr Ingenieure;                       *
C                G. Engeln-MÅllges/F. Reutter, 4th ed. 1984;                 *
C                Example 6.3 (III), p. 213 - 215.                            *
C                                                                            *
C[   I *  X(I)  *  F(I)  *  W(I)  *                                         ]*
C[  -------------------------------                                         ]*
C[   0 *    .02 *  50.00 *   1.00 *                                         ]*
C[   1 *    .10 *  10.00 *   1.00 *                                         ]*
C[   2 *    .50 *   1.00 *   1.00 *                                         ]*
C[   3 *   1.00 *    .00 *   1.00 *                                         ]*
C[                                                                          ]*
C[  MAX. DEGREE OF THE LEAST SQUARES POLYNOMIAL, N =  3                     ]*
C[                                                                          ]*
C[  IFEHL =  0                                                              ]*
C[                                                                          ]*
C[   J *         C(J)         *                                             ]*
C[  ---------------------------                                             ]*
C[   0 *      62.981434240377 *                                             ]*
C[   1 *    -680.869756236137 *                                             ]*
C[   2 *    1609.739229025772 *                                             ]*
C[   3 *    -991.850907030012 *                                             ]*
C[                                                                          ]*
C[  TABLE OF VALUES                                                         ]*
C[  ---------------                                                         ]*
C[                                                                          ]*
C[    .00        62.981434240377                                            ]*
C[    .10         9.999999420817                                            ]*
C[    .20       -16.737755566964                                            ]*
C[    .30       -23.182936431181                                            ]*
C[    .40       -15.286648880048                                            ]*
C[    .50         1.000001378218                                            ]*
C[    .60        19.725908635403                                            ]*
C[    .70        34.939967183292                                            ]*
C[    .80        40.691071313669                                            ]*
C[    .90        31.028115318319                                            ]*
C[   1.00         -.000006510973                                            ]*
C[   1.10       -58.344399882423                                            ]*
C[   1.20      -149.956170504246                                            ]*
C[   1.30      -280.786424084656                                            ]*
C[   1.40      -456.786266331869                                            ]*
C[   1.50      -683.906802954101                                            ]*
C[   1.60      -968.099139659566                                            ]*
C[   1.70     -1315.314382156479                                            ]*
C[   1.80     -1731.503636153057                                            ]*
C[   1.90     -2222.618007357514                                            ]*
C[   2.00     -2794.608601478065                                            ]*
C                                                                            *
C*****************************************************************************
C
      PARAMETER (IA=4,N=3,M=3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(0:M),F(0:M),W(0:M),C(0:N),A(IA,N+1),B(N+1),Y(N+1),
     +     Z(N+1)
      DATA X,F,W/.02D0,.1D0,.5D0,1.D0,50.D0,10.D0,1.D0,0.D0,4*1.D0/
      WRITE(*,900)
      DO 10 I=0,M
         WRITE(*,910)I,X(I),F(I),W(I)
   10 CONTINUE
      WRITE(*,980)'MAX. DEGREE OF THE LEAST SQUARES POLYNOMIAL, N = ',N
      CALL GADESM(M,X,F,W,IA,N,C,A,B,Y,Z,IFEHL)
      WRITE(*,980)'IFEHL = ',IFEHL
      IF(IFEHL.EQ.0)THEN
         WRITE(*,920)
         DO 20 I=0,N
            WRITE(*,930)I,C(I)
   20    CONTINUE
         WRITE(*,940)
         DO 30 XX=0.,2.05,.1
            WRITE(*,950)XX,APPROX(N,C,XX)
   30    CONTINUE
      ELSE
         WRITE(*,970)'***** ERROR ***** SEE PROGRAM DESCRIPTION !'
      ENDIF
      STOP
  900 FORMAT(1X,'C[',2X,' I *  X(I)  *  F(I)  *  W(I)  *',T78,']*',/,
     +       1X,'C[',2X,31('-'),T78,']*')
  910 FORMAT(1X,'C[',3X,I1,' *',3(1X,F6.2,' *'),T78,']*')
  920 FORMAT(1X,'C[',T78,']*',/,
     +       1X,'C[',2X,' J *         C(J)         *',T78,']*',/,
     +       1X,'C[',2X,27('-'),T78,']*')
  930 FORMAT(1X,'C[',3X,I1,' * ',F20.12,' *',T78,']*')
  940 FORMAT(1X,'C[',T78,']*',/,
     +       1X,'C[',2X,'TABLE OF VALUES',T78,']*',/,
     +       1X,'C[',2X,'---------------',T78,']*',/,
     +       1X,'C[',T78,']*')
  950 FORMAT(1X,'C[',2X,F5.2,3X,F20.12,T78,']*')
  970 FORMAT(1X,'C[',T78,']*',/,
     +       1X,'C[',2X,A,T78,']*')
  980 FORMAT(1X,'C[',T78,']*',/,
     +       1X,'C[',2X,A,I2,T78,']*')
      END
C
C
      DOUBLE PRECISION FUNCTION APPROX(N,C,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION C(0:N)
      APPROX=0.D0
      DO 10 I=N,0,-1
         APPROX=APPROX*X+C(I)
   10 CONTINUE
      RETURN
      END
