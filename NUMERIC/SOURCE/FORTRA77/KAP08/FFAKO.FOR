C[BA*)
C[LE*)
       SUBROUTINE FFAKO (M, FRE, FIM, HRE, HIM)
C[IX{FFAKO}*)
C
C*****************************************************************
C                                                                *
C  This program uses the fast Fourier transform (FFT) to compute *
C  the discrete values of the convolution                        *
C[BE*)
C                                                                *
C      Falt(j) = 1/M * (SUM k=0 to M-1) F(j-k)*H(k)              *
C                                                                *
C[BA*)
C  and of the discrete cyclic correlation                        *
C[BE*)
C                                                                *
C      Korr(j) = 1/M * (SUM k=0 to M-1) F(j+k)*CONJG(H(k))       *
C                                                                *
C[BA*)
C  of F and H at the given complex functional values             *
C  F(0), F(1), ... F(M-1) and  H(0), H(1), ... , H(M-1)  at      *
C  equidistant nodes in the period interval for                  *
C  j = 0, 1, ... , M-1 .                                         *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  M        : Number of nodes.                                   *
C             If M is a power of two ( M = 2**ITAU  for a        *
C             positive integer ITAU ), we make use of the        *
C             subroutineso  FFT, which runs for all values of M. *
C             If M  is not a power of two, then we call on FFTB, * 
C             which limits the available size for its auxiliary  *
C             vectors F1RE, F1IM, GRE, GIM  to  1 <= M <= 1366.  *
C             For larger M, the dimension statement inside FFTB  *
C             must be amended as indicated there.                *
C  FRE, FIM : DOUBLE PRECISION vectors for the real and          *
C             imaginary parts for M functional values of F.      * 
C  HRE, HIM : DOUBLE PRECISION vectors for the real and          *
C             imaginary parts for M functional values of H.      *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  FRE, FIM : DOUBLE PRECISION vectors with real and imaginary   *
C             parts of the discrete cyclic convolution Falt(j).  *
C  HRE, HIM : DOUBLE PRECISION vectors with real and imaginary   *
C             parts of the discrete cyclic correlation Korr(j).  *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Subroutines used: FFT or FFTB                                 *
C                                                                *
C*****************************************************************
C                                                                *
C  Author      : Klaus Niederdrenk                               *
C  Date        : 09.08.1994                                      *
C  Source code : FORTRAN 77                                      *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER M, ITAU
      DIMENSION FRE(0:M-1), FIM(0:M-1), HRE(0:M-1), HIM(0:M-1)
      LOGICAL POT2
C
C     complex multiplication: Real part
C
      CMLR(AR,AI,BR,BI)=AR*BR-AI*BI
C
C     complex multiplication: Imaginary part
C
      CMLI(AR,AI,BR,BI)=AR*BI+AI*BR
C
C     Check whether M is a power of two
C
      ITAU=NINT(LOG(REAL(M))/LOG(2.0))
      IF(2**ITAU .EQ. M) THEN
        POT2=.TRUE.
      ELSE
        POT2=.FALSE.
      ENDIF
C
C*****************************************************************
C     Compute the discrete Fourier transformation                *
C*****************************************************************
C
      IF( POT2 ) THEN
        CALL FFT(ITAU,FRE,FIM,0)
        CALL FFT(ITAU,HRE,HIM,0)
      ELSE
        CALL FFTB(M,FRE,FIM,0)
        CALL FFTB(M,HRE,HIM,0)
      ENDIF
C
C*****************************************************************
C     Conpute the transform of the corresponding discrete cyclic *
C     convolution and of the discrete cyclic correlation         * 
C*****************************************************************
C
      DO 10 K = 0, M-1
          H1    =CMLR(FRE(K),FIM(K),HRE(K),HIM(K)) 
          H2    =CMLI(FRE(K),FIM(K),HRE(K),HIM(K))
          H     =CMLR(FRE(K),FIM(K),HRE(K),-HIM(K))
          HIM(K)=CMLI(FRE(K),FIM(K),HRE(K),-HIM(K))
          HRE(K)=H
          FRE(K)=H1
          FIM(K)=H2
  10  CONTINUE
C
C*****************************************************************
C     Evaluate the functional values                             *
C*****************************************************************
C
      IF( POT2 ) THEN
        CALL FFT(ITAU,FRE,FIM,1)
        CALL FFT(ITAU,HRE,HIM,1)
      ELSE                                                  
        CALL FFTB(M,FRE,FIM,1)
        CALL FFTB(M,HRE,HIM,1)
      ENDIF
      RETURN
      END
