      PROGRAM  TEST
C                                                ( Thomas Meuser )
C*****************************************************************************
C                                                                            *
C     Testprogram for the subroutines BAND, BANDP, BANDS.                    *
C     Solve a system of linear equations  A*X = B with a banded matrix A.    *
C                                                                            *
C     For the included test examples the following results are computed:     *
C                                                                            *
C[  EXAMPLE:                                                                ]*
C[  ========                                                                ]*
C[  COEFFICIENT MATRIX A:                                                   ]*
C[                                                                          ]*
C[   .50000D+01   .10000D+01   .10000D+01   .00000D+00   .00000D+00         ]*
C[   .20000D+01   .60000D+01   .20000D+01   .00000D+00   .00000D+00         ]*
C[   .20000D+01   .20000D+01   .70000D+01   .00000D+00   .20000D+01         ]*
C[   .00000D+00   .10000D+01   .30000D+01   .80000D+01   .20000D+01         ]*
C[   .00000D+00   .00000D+00   .30000D+01   .10000D+01   .90000D+01         ]*
C[                                                                          ]*
C[  RIGHT HAND SIDE                                                         ]*
C[                                                                          ]*
C[   .70000D+01   .10000D+02   .13000D+02   .14000D+02   .13000D+02         ]*
C[                                                                          ]*
C[                                                                          ]*
C[  MATRIX IN PACKED FORM:                                                  ]*
C[                                                                          ]*
C[                                                                          ]*
C[   .00000D+00   .00000D+00   .50000D+01   .10000D+01   .10000D+01         ]*
C[   .00000D+00   .20000D+01   .60000D+01   .20000D+01   .00000D+00         ]*
C[   .20000D+01   .20000D+01   .70000D+01   .00000D+00   .20000D+01         ]*
C[   .10000D+01   .30000D+01   .80000D+01   .20000D+01   .00000D+00         ]*
C[   .30000D+01   .10000D+01   .90000D+01   .00000D+00   .00000D+00         ]*
C[                                                                          ]*
C[  SOLUTION                                                                ]*
C[                                                                          ]*
C[   .10000D+01   .10000D+01   .10000D+01   .10000D+01   .10000D+01         ]*
C[  STOP. NO ERROR !                                                        ]*
C                                                                            *
C     Further tests can be performed by changing the data.                   *
C                                                                            *
C*****************************************************************************
C
      PARAMETER (N=5,LDAP=N)
      PARAMETER (MU=2,MO=2,MB=MU+MO+1+MU+1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A (N,N),AP(LDAP,MB)
      DOUBLE PRECISION B (N)
      INTEGER IP(N)
C
C   Initialize; for other test examples, change the input data.
C
      DATA ((A(J,I),I=1,N),J=1,N)/5.D0,1.D0,1.D0,0.D0,0.D0,2.D0,6.D0
     +,2.D0,0.D0,0.D0,2.D0,2.D0,7.D0,0.D0,2.D0,0.D0,1.D0,3.D0,8.D0,2.D0
     +,0.D0,0.D0,3.D0,1.D0,9.D0/
      DATA (B(I),I=1,N) /7.D0,10.D0,13.D0,14.D0,13.D0/
C
C   Output of test example
C
      WRITE (*,2000)
      DO 10 I=1,N
   10 WRITE (*,2010)  (A(I,J),J=1,N)
      WRITE (*,2020)  (B(I),I=1,N)
C
C   Compute the packed matrix AP to use in BAND.
C
      DO 20 I=1,N
         DO 20 K=MAX(1,I-MU),MIN(N,I+MO)
         AP(I,MU+1+K-I)=A(I,K)
   20 CONTINUE
C
C   Output of the condensed matrix AP
C
      WRITE (*,2023)
      DO 30 I=1,N
         WRITE (*,2025) (0.D0,J=1,MAX(0,MU+1-I)),
     +                  (AP(I,J),J=MAX(0,MU+1-I)+1,
     +                             MIN(MU+1+MO,N+MU+1-I)),
     +                  (0.D0,J=MIN(MU+1+MO,N+MU+1-I)+1,MU+1+MO)
   30 CONTINUE
C
      CALL BAND (AP,LDAP,MB,N,MU,MO,B,IFEHL,IP)
C
      IF (IFEHL.EQ.0) THEN
C
C   Output of a solution or an error message according to IFEHL
C
        WRITE (*,2030) (B(I),I=1,N)
        WRITE (*,2040) 'STOP. NO ERROR !'
        STOP
       ELSE
        WRITE (*,2040) 'STOP. ERROR --- MATRIX IS SINGULAR !'
        STOP
      ENDIF
C
 2000 FORMAT (1X,'C[',2X,'EXAMPLE:',T78,']*',/,
     +        1X,'C[',2X,8('='),T78,']*',/,
     +        1X,'C[',2X,'COEFFICIENT MATRIX A:',T78,']*',/,
     +        1X,'C[',T78,']*')
 2010 FORMAT (1X,'C[',5(1X,D12.5),T78,']*')
 2020 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[',2X,'RIGHT HAND SIDE',T78,']*',/,
     +        1X,'C[',T78,']*',/,
     +        1X,'C[',5(1X,D12.5),T78,']*',/,
     +        1X,'C[',T78,']*',/,1X,'C[',T78,']*')
 2023 FORMAT (1X,'C[',2X,'MATRIX IN PACKED FORM:',T78,']*',/,
     +        1X,'C[',T78,']*',/,1X,'C[',T78,']*')
 2025 FORMAT (1X,'C[',5(1X,D12.5),T78,']*')
 2030 FORMAT (1X,'C[',T78,']*',/,
     +        1X,'C[',2X,'SOLUTION',T78,']*',/,
     +        1X,'C[',T78,']*',/,
     +        1X,'C[',5(1X,D12.5),T78,']*')
2040  FORMAT (1X,'C[',2X,A,T78,']*')
      END
