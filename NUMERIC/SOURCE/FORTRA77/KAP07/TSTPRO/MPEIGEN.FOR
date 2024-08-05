      PROGRAM TEST
C                                                                     18/12/87
C                                              (Thomas Meuser)
C*****************************************************************************
C                                                                            *
C     Testprogram for the subroutine EIGEN.                                  *
C     Compute all eigenvalues and eigenvectors of a matrix A with the QR     *
C     method according to Martin, Parlett, Peters, Reinsch and Wilkinson.    *
C                                                                            *
C     The test gives the following results:                                  *
C                                                                            *
C[  TEST EXAMPLE:                                                           ]*
C[  =============                                                           ]*
C[  GIVEN MATRIX:                                                           ]*
C[                                                                          ]*
C[    -.20000D+01   .20000D+01   .20000D+01   .20000D+01                    ]*
C[    -.30000D+01   .30000D+01   .20000D+01   .20000D+01                    ]*
C[    -.20000D+01   .00000D+00   .40000D+01   .20000D+01                    ]*
C[    -.10000D+01   .00000D+00   .00000D+00   .50000D+01                    ]*
C[                                                                          ]*
C[  COMPUTED EIGENVALUES:                                                   ]*
C[                                                                          ]*
C[            REAL PART   I IMAGINARY PART                                  ]*
C[          ------------------------------                                  ]*
C[           .10000D+01   I     .00000D+00                                  ]*
C[           .20000D+01   I     .00000D+00                                  ]*
C[           .30000D+01   I     .00000D+00                                  ]*
C[           .40000D+01   I     .00000D+00                                  ]*
C[                                                                          ]*
C[  NORMALIZED EIGENVECTORS:                                                ]*
C[                                                                          ]*
C[     .10000D+01   .10000D+01   .10000D+01   .10000D+01                    ]*
C[     .75000D+00   .10000D+01   .10000D+01   .10000D+01                    ]*
C[     .50000D+00   .66667D+00   .10000D+01   .10000D+01                    ]*
C[     .25000D+00   .33333D+00   .50000D+00   .10000D+01                    ]*
C[  STOP. NO ERROR !                                                        ]*
C                                                                            *
C     Other tests with different data are possible.                          *
C                                                                            *
C*****************************************************************************
C
      PARAMETER (N=4,ND=N)
      INTEGER CNT (N),EIGEN,BASIS,RES,LOW,HIGH
      DOUBLE PRECISION MAT(ND,N),WERTR(N),WERTI(N),EIVEC(ND,N),SKAL(N),
     +                 D(N)
      LOGICAL VEC,ORTHO,EVNORM  
C
C     Initialize; change if desired
C
      DATA MAT/-2.0D+00,-3.0D+00,-2.0D+00,-1.0D+00,2.0D+00,3.0D+00,
     +          0.0D+00,0.0D+00,2.0D+00,2.0D+00,4.0D+00,0.0D+00,
     +          2.0D+00,2.0D+00,2.0D+00,5.0D+00/
      DATA BASIS/2/
      DATA VEC,ORTHO,EVNORM /.TRUE.,.FALSE.,.TRUE./
C
C     Print out test example
C
      WRITE(*,2000)
      DO 10 I=1,N
   10 WRITE(*,2010)(MAT(I,J),J=1,N)
C
C
      RES = EIGEN(VEC,ORTHO,EVNORM,BASIS,ND,N,MAT,SKAL,D,EIVEC,WERTR,
     +            WERTI,CNT,LOW,HIGH)
C
C
      IF (RES .GT. 0) THEN
          RES = RES-400
C
C     Output of error code in  RES (0/401/402/403)
C
      GOTO (30,40,50),RES
   30 WRITE (*,2050)
     +      'ERROR: SIZE N OF THE MATRIX IS SMALLER THAN 1 !'
      STOP
   40 WRITE (*,2050) 'ERROR: MAT IS THE ZERO MATRIX !'
      STOP
   50 WRITE (*,2050)
     +      'ERROR: MAX NUMBER OF STEPS FOR QR-METHOD IS REACHED!'
      STOP
      ELSE
C
C     Output
C
      WRITE(*,2015)
      WRITE(*,2020)(WERTR(I),WERTI(I),I=1,N)
      WRITE(*,2030)
      DO 20 I=1,N
   20 WRITE(*,2040)(EIVEC(I,J),J=1,N)
      WRITE (*,2050) 'STOP. NO ERROR !'
      STOP
      ENDIF
C
C
 2000 FORMAT(1X,'C[',2X,'TEST EXAMPLE:',T78,']*',/,
     +       1X,'C[',2X,13('='),T78,']*',/,
     +       1X,'C[',2X,'GIVEN MATRIX:',T78,']*',/,
     +       1X,'C[',T78,']*')
 2010 FORMAT(1X,'C[',2X,4(1X,D12.5),T78,']*')
 2015 FORMAT(1X,'C[',T78,']*',/,
     +       1X,'C[',2X,'COMPUTED EIGENVALUES:',T78,']*',/,
     +       1X,'C[',T78,']*',/,
     +       1X,'C[',12X,'REAL PART   I ','IMAGINARY PART',T78,']*',/,
     +       1X,'C[',10X,30('-'),T78,']*')
 2020 FORMAT(1X,'C[',7X,D14.5,'   I ',D14.5,T78,']*')
 2030 FORMAT(1X,'C[',T78,']*',/,
     +       1X,'C[',2X,'NORMALIZED EIGENVECTORS:',T78,']*',/,
     +       1X,'C[',T78,']*')
 2040 FORMAT(1X,'C[',2X,4(1X,D12.5),T78,']*')
 2050 FORMAT (1X,'C[',2X,A,T78,']*')
      END
