      PROGRAM TEST
C                                                ( Thomas Meuser )
C**************************************************************************
C                                                                         *
C     Testprogram for subroutine  PIVOT.                                  *
C     Compute a matrix inverse via exchange steps.                        *
C                                                                         *
C     The test example produces this output:                              *
C                                                                         *
C[  EXAMPLE:                                                             ]*
C[  ========                                                             ]*
C[                                                                       ]*
C[  INPUT PARAMETER:                                                     ]*
C[  IA= 4    N= 4                                                        ]*
C[                                                                       ]*
C[  INPUT MATRIX A(IA,N):                                                ]*
C[                                                                       ]*
C[     3.00000       1.00000       2.00000       5.00000                 ]*
C[     2.00000       1.00000       3.00000       7.00000                 ]*
C[     3.00000       1.00000       2.00000       4.00000                 ]*
C[     4.00000       1.00000       3.00000       2.00000                 ]*
C[                                                                       ]*
C[  INVERSE MATRIX B(IA,N):                                              ]*
C[                                                                       ]*
C[     2.50000       -.50000      -2.50000        .50000                 ]*
C[   -10.50000        .50000      13.50000      -2.50000                 ]*
C[     -.50000        .50000       -.50000        .50000                 ]*
C[     1.00000        .00000      -1.00000        .00000                 ]*
C[                                                                       ]*
C[  CONTROL PARAMETERS:                                                  ]*
C[  S1=    .00000    S2=    .00000                                       ]*
C[  IERR= 1                                                              ]*
C[                                                                       ]*
C[  AUX ARRAY MX(1:N):                                                   ]*
C[   1     2     3     4                                                 ]*
C[                                                                       ]*
C[  AUX ARRAY MY(1:N):                                                   ]*
C[   1     2     3     4                                                 ]*
C                                                                         *
C     Other test examples can be run for different data.                  *
C                                                                         *
C**************************************************************************
C
      PARAMETER (IA=4,N=4)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A(IA,N),B(IA,N),S1,S2,WERT
      INTEGER MX(N),MY(N)
C
C     Initialize data; change problem, if desired
C
      DATA A/3.D0,2.D0,3.D0,4.D0,4*1.D0,2.D0,3.D0,2.D0,3.D0,5.D0,7.D0
     +,4.D0,2.D0/
      WRITE(*,899)
      WRITE(*,900)IA,N
      WRITE(*,910)
      DO 10 I=1,IA
         WRITE(*,920)(A(I,K),K=1,N)
   10 CONTINUE
C
      CALL PIVOT(A,IA,N,B,S1,S2,IERR,MX,MY,WERT)
C
C     Output of result
C
      IF(IERR.EQ.1)THEN
         WRITE(*,930)
         DO 20 I=1,IA
            WRITE(*,920)(B(I,K),K=1,N)
   20    CONTINUE
         WRITE(*,940)S1,S2,IERR
         WRITE(*,950)(MX(I),I=1,N)
         WRITE(*,960)(MY(I),I=1,N)
      ELSE
         WRITE(*,970)WERT,IERR
      END IF
      STOP
C
  899 FORMAT (1X,'C[',2X,'EXAMPLE:',T75,']*',/,
     +        1X,'C[',2X,8('='),T75,']*',/,
     +        1X,'C[',T75,']*')
  900 FORMAT (1X,'C[',2X,'INPUT PARAMETER:',T75,']*',/,
     +        1X,'C[',2X,'IA=',I2,4X,'N=',I2,T75,']*',/,
     +        1X,'C[',T75,']*')
  910 FORMAT (1X,'C[',2X,'INPUT MATRIX A(IA,N):',T75,']*',/,
     +        1X,'C[',T75,']*')
  920 FORMAT( 1X,'C[',2X,4(F10.5,4X),T75,']*')
  930 FORMAT( 1X,'C[',T75,']*',/,
     +        1X,'C[',2X,'INVERSE MATRIX B(IA,N):',T75,']*',/,
     +        1X,'C[',T75,']*')
  940 FORMAT (1X,'C[',T75,']*',/,
     +        1X,'C[',2X,'CONTROL PARAMETERS:',T75,']*',/,
     +        1X,'C[',2X,'S1=',F10.5,4X,'S2=',F10.5,T75,']*',/,
     +        1X,'C[',2X,'IERR=',I2,T75,']*')
  950 FORMAT (1X,'C[',T75,']*',/,
     +        1X,'C[',2X,'AUX ARRAY MX(1:N):',T75,']*',/,
     +        1X,'C[',2X,4(I2,4X),T75,']*')
  960 FORMAT (1X,'C[',T75,']*',/,
     +        1X,'C[',2X,'AUX ARRAY MY(1:N):',T75,']*',/,
     +        1X,'C[',2X,4(I2,4X),T75,']*')
  970 FORMAT (1X,'C[',T75,']*',/,
     +        1X,'C[',2X,'PIVOT VALUE= ',D20.14,T75,']*',/,
     +        1X,'C[',2X,'IERR=',I2,T75,']*',/,
     +        1X,'C[',2X,'MATRIX A IS SINGULAR, AN INVERSE',T75,']*',/,
     +        1X,'C[',2X,'DOES NOT EXIST',T75,']*')
      END
