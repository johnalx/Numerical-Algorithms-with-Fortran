      PROGRAM TEST
C
C*****************************************************************************
C                                                                            *
C     Testprogram for the condition number subroutine CONDAE.                *
C                                                                            *
C     The test data produces the following output:                           *
C                                                                            *
C[  EXAMPLE:                                                                ]*
C[  ========                                                                ]*
C[  MATRIX AO:                                                              ]*
C[                                                                          ]*
C[   .10000D+01   .50000D+00   .33333D+00   .25000D+00   .20000D+00         ]*
C[   .50000D+00   .33333D+00   .25000D+00   .20000D+00   .16667D+00         ]*
C[   .33333D+00   .25000D+00   .20000D+00   .16667D+00   .14286D+00         ]*
C[   .25000D+00   .20000D+00   .16667D+00   .14286D+00   .12500D+00         ]*
C[   .20000D+00   .16667D+00   .14286D+00   .12500D+00   .11111D+00         ]*
C[                                                                          ]*
C[  RIGHT HAND SIDE Y:                                                      ]*
C[                                                                          ]*
C[   .22833D+01   .14500D+01   .10929D+01   .88452D+00   .74563D+00         ]*
C[                                                                          ]*
C[                                                                          ]*
C[  MATRIX A:                                                               ]*
C[                                                                          ]*
C[   .10000D+01   .50000D+00   .33333D+00   .25000D+00   .20000D+00         ]*
C[   .20000D+00   .66667D-01   .76190D-01   .75000D-01   .71111D-01         ]*
C[   .50000D+00   .12500D+01  -.11905D-01  -.18750D-01  -.22222D-01         ]*
C[   .33333D+00   .12500D+01   .53333D+00  -.41667D-03  -.84656D-03         ]*
C[   .25000D+00   .11250D+01   .20000D+00   .64286D+00  -.11338D-04         ]*
C[                                                                          ]*
C[  SOLUTION VECTOR X:                                                      ]*
C[                                                                          ]*
C[   .10000D+01   .10000D+01   .10000D+01   .10000D+01   .10000D+01         ]*
C[                                                                          ]*
C[                                                                          ]*
C[  ESTIMATE OF THE CONDITION NUMBER OF A=-.40544090000000D+32 ;   N=  5    ]*
C                                                                            *
C*****************************************************************************
C
      PARAMETER(N=5,IA=N)
      DOUBLE PRECISION A(IA,N),AO(IA,N),Y(N),X(N),D(N),Z(N),R(N)
      INTEGER IPIVOT(N)
C
      DO 10 I=1,N
         DO 10 K=1,N
            A(I,K)=1.0D+00/DBLE(I+K-1)
            AO(I,K)=A(I,K)
10    CONTINUE
      DO 20 I=1,N
         Y(I)=0.0D+00
         DO 20 K=1,N
            Y(I)=Y(I)+A(I,K)
20    CONTINUE
C
C Output matrix AO
C
      WRITE(*,2000)
      DO 30 I=1,5
        WRITE(*,2020)(AO(I,J),J=1,5)
30    CONTINUE
      WRITE(*,2030)(Y(I),I=1,5)
C
      CALL GAUSS(N,A,IA,Y,X,MARKE,D,IPIVOT)
C
C Output after SUBROUTINE GAUSS
C
      IF(MARKE.EQ.0) THEN
        WRITE(*,2050)'MATRIX SINGULAR',' N=',N
        STOP
      ELSE
        WRITE(*,2010)
        DO 40 I=1,5
          WRITE(*,2020)(A(I,J),J=1,5)
40      CONTINUE
        WRITE(*,2040)(X(I),I=1,5)
      ENDIF
C
      CALL CONDAE(N,AO,A,IA,IPIVOT,Y,X,Z,R,CONDA)
C
      WRITE(*,2060)'ESTIMATE OF THE CONDITION NUMBER OF A=',CONDA,
     +             ' ;   N=',N
      STOP
 2000 FORMAT (1X,'C[',2X,'EXAMPLE:',T78,']*',/,
     +        1X,'C[',2X,8('='),T78,']*',/,
     +        1X,'C[',2X,'MATRIX A0:',T78,']*',/,
     +        1X,'C[',T78,']*')
 2010 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[',2X,'MATRIX A:',T78,']*',/,
     +        1X,'C[',T78,']*')
 2020 FORMAT (1X,'C[',5(1X,D12.5),T78,']*')
 2030 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[',2X,'RIGHT HAND SIDE Y:',T78,']*',/,
     +        1X,'C[',T78,']*',/,
     +        1X,'C[',5(1X,D12.5),T78,']*',/,
     +        1X,'C[',T78,']*')
 2040 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[',2X,'SOLUTION VECTOR X:',T78,']*',/,
     +        1X,'C[',T78,']*',/,
     +        1X,'C[',5(1X,D12.5),T78,']*',/,
     +        1X,'C[',T78,']*',/,1X,'C[',T78,']*')
 2050 FORMAT (1X,'C[',2X,A,I3,T78,']*')
 2060 FORMAT (1X,'C[',2X,A,D20.14,A,I3,T78,']*')
      END
