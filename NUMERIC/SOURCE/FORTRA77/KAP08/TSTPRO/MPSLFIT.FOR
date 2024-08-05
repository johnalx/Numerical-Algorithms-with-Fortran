C                                                         89.01.14
      PROGRAM TEST
C
C*****************************************************************************
C                                                                            *
C  Testprogram for the SUBROUTINES SLFIT, SLPRE.                             *
C                                                                            *
C  Test example: Numerische Mathematik fÅr Ingenieure;                       *
C                G. Engeln-MÅllges/F. Reutter, 5th ed. 1987;                 *
C                EXAMPLE 6.3 (II), p. 213 - 215.                             *
C                                                                            *
C  Results from a PC:                                                        *
C                                                                            *
C[   I *  X(I)  *  F(I)  *  W(I)  *                                         ]*
C[  -------------------------------                                         ]*
C[   0 *    .02 *  50.00 *   1.00 *                                         ]*
C[   1 *    .10 *  10.00 *   1.00 *                                         ]*
C[   2 *    .50 *   1.00 *   1.00 *                                         ]*
C[   3 *   1.00 *    .00 *   1.00 *                                         ]*
C[                                                                          ]*
C[  NUMBER OF FUNCTIONS, N+1 =  3                                           ]*
C[                                                                          ]*
C[  IFEHL =  0                                                              ]*
C[                                                                          ]*
C[   J *         C(J)         *                                             ]*
C[  ---------------------------                                             ]*
C[   0 *      39.678891125667 *                                             ]*
C[   1 *    -136.551407186668 *                                             ]*
C[   2 *      97.982953935049 *                                             ]*
C[                                                                          ]*
C[  TABLE OF VALUES:                                                        ]*
C[  ----------------                                                        ]*
C[                                                                          ]*
C[    .00        39.678891125667                                            ]*
C[    .10        27.003579946350                                            ]*
C[    .20        16.287927845735                                            ]*
C[    .30         7.531934823821                                            ]*
C[    .40          .735600880607                                            ]*
C[    .50        -4.101073983905                                            ]*
C[    .60        -6.978089769717                                            ]*
C[    .70        -7.895446476827                                            ]*
C[    .80        -6.853144105237                                            ]*
C[    .90        -3.851182654945                                            ]*
C[   1.00         1.110437874047                                            ]*
C[   1.10         8.031717481741                                            ]*
C[   1.20        16.912656168135                                            ]*
C[   1.30        27.753253933231                                            ]*
C[   1.40        40.553510777027                                            ]*
C[   1.50        55.313426699524                                            ]*
C[   1.60        72.033001700722                                            ]*
C[   1.70        90.712235780622                                            ]*
C[   1.80       111.351128939222                                            ]*
C[   1.90       133.949681176523                                            ]*
C[   2.00       158.507892492525                                            ]*
C[                                                                          ]*
C[  AVERAGE LEAST SQUARES ERROR =  .22038885324468E+02                      ]*
C                                                                            *
C*****************************************************************************
C
      PARAMETER (IA=3,N=2,M=3,IANZ=21)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(0:M),Y(0:M),W(0:M),A(0:IA,0:N+1),D(0:N),C(0:N),
     +     FW(1:IANZ),PSI(1:IANZ)
      EXTERNAL FKT
      DATA X,Y,W/.02D0,.1D0,.5D0,1.D0,50.D0,10.D0,1.D0,0.D0,4*1.D0/
      WRITE(*,900)
      DO 10 I=0,M
         WRITE(*,910)I,X(I),Y(I),W(I)
   10 CONTINUE
      WRITE(*,980)'NUMBER OF FUNCTIONS, N+1 = ',N+1
      XX=0.
      DO 20 I=1,IANZ
         PSI(I)=XX
         XX=XX+.1D0
   20 CONTINUE
      CALL SLFIT(X,Y,W,0,FKT,PSI,IA,M,N,IANZ,A,D,C,FW,QUADFE,IFEHL)
      WRITE(*,980)'IFEHL = ',IFEHL
      IF(IFEHL.EQ.0)THEN
         WRITE(*,920)
         DO 30 I=0,N
            WRITE(*,930)I,C(I)
   30    CONTINUE
         WRITE(*,940)
         DO 40 I=1,IANZ
            WRITE(*,950)PSI(I),FW(I)
   40    CONTINUE
         WRITE(*,960)'AVERAGE LEAST SQUARES ERROR = ',QUADFE
      ELSE
         WRITE(*,970)'***** ERROR ***** SEE PROGRAM DESCRIPTION!'
      END IF
      STOP
  900 FORMAT(1X,'C[',2X,' I *  X(I)  *  F(I)  *  W(I)  *',T78,']*',/,
     +       1X,'C[',2X,31('-'),T78,']*')
  910 FORMAT(1X,'C[',3X,I1,' *',3(1X,F6.2,' *'),T78,']*')
  920 FORMAT(1X,'C[',T78,']*',/,
     +       1X,'C[',2X,' J *         C(J)         *',T78,']*',/,
     +       1X,'C[',2X,27('-'),T78,']*')
  930 FORMAT(1X,'C[',3X,I1,' * ',F20.12,' *',T78,']*')
  940 FORMAT(1X,'C[',T78,']*',/,
     +       1X,'C[',2X,'TABLE OF VALUES:',T78,']*',/,
     +       1X,'C[',2X,'----------------',T78,']*',/,
     +       1X,'C[',T78,']*')
  950 FORMAT(1X,'C[',2X,F5.2,3X,F20.12,T78,']*')
  960 FORMAT(1X,'C[',T78,']*',/,
     +       1X,'C[',2X,A,E20.14,T78,']*')
  970 FORMAT(1X,'C[',T78,']*',/,
     +       1X,'C[',2X,A,T78,']*')
  980 FORMAT(1X,'C[',T78,']*',/,
     +       1X,'C[',2X,A,I2,T78,']*')
      END
C
C
      SUBROUTINE FKT(X,N,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION F(0:N)
      F(0)=1.D0
      F(1)=X
      F(2)=X*X
      RETURN
      END
