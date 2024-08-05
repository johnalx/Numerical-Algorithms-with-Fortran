C[BA*)
C[LE*)
      SUBROUTINE FFT (ITAU, FRE, FIM, IR)
C[IX{FFT}*)
C
C*****************************************************************
C                                                                *
C  If IR = 0, this programm uses the fast Fourier transform      *
C  ( F F T ) in order to determine the discrete Fourier          *
C  coefficients   F^(-M/2), ... , F^(M/2-1)  for  M = 2**ITAU    *
C  given real or complex functional values F(0), F(1), ... ,     *
C  F(M-1).                                                       *
C[BE*)
C  These coefficients define the discrete Fourier series         *
C                                                                *
C     (SUM from K=-M/2 to M/2-1)  F^(K)*EXP(I*K*OMEGA*X)         *
C                                                                *
C  ( I : imaginary unit with  I**2 = -1 ; OMEGA = 2*PI/L ,       *
C    L : period ) .                                              *
C[BA*)
C  For  IR = 1  the inverse Fourier transform is performed.      *
C[BE*)
C  The double precision vectors  FRE  and  FIM  contain the real *
C  and imaginary parts of  F  and  F^, respectively.             *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  ITAU   : the number of data points is  M = 2**ITAU            *
C  FRE, FIM : DOUBLE PRECISION vectors for the real and imaginary*
C             parts of M values.                                 *
C             Depending on IR, FRE, FIM  must contain the        *
C             following:                                         *
C             IR=0: FRE, FIM  contain the real and imaginary     *
C                   parts of the functional values               *
C             IR=1: FRE, FIM  contain the real and imaginary     *
C                   parts of the discrete Fourier coefficients   *
C                   as follows:                                  *
C                   FRE(K), FIM(K) : Real and imaginary parts of *
C                   F^(K)  for  K = 0 to M/2-1  and of F^(K-M)   *
C                   for  K = M/2 to M-1.                         *
C  IR     : determines the direction of the desired              *
C           transfornation                                       *
C           If IR=0, we compute the discrete Fourier             *
C                    coefficients                                *
C           If IR=1: we compute the function values              *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  FRE, FIM : DOUBLE PRECISION vectors for the real and          *
C             imaginary parts of M complex numbers.              *
C             Depending on IR the following values are returned: *
C             If IR=0: The vectors FRE and FIM contain the       *
C                      real and imaginary parts of the discrete  *
C                      Fourier coefficients as follows:          *
C                      FRE(K+M), FIM(K+M): The real and imaginary*
C                           parts of  F^(K) for  K = -M/2 to -1  *
C                      FRE(K), FIM(K) : The real and imaginary   *
C                           parts of  F^(K)  for  K = 0 to M/2-1 *
C             If IR=1: The functional values are returned in FRE *
C                      and FIM.                                  *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Required subroutines: none                                    *
C                                                                *
C  Reference: Niederdrenk, K., see [NIED84]                      *
C                                                                *
C*****************************************************************
C                                                                *
C  Author    : Klaus Niederdrenk                                 *
C  Date      : 03.22.1993                                        *
C  Source    : FORTRAN 77                                        *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER ITAU, IR, SIGMA, M
      DIMENSION FRE(0:2**ITAU-1), FIM(0:2**ITAU-1)
C
C     Define FUNCTION for complex multiplication: real part
C
      CMLR(AR,AI,BR,BI)=AR*BR-AI*BI
C
C     Define FUNCTION for complex multiplication: real part
C
      CMLI(AR,AI,BR,BI)=AR*BI+AI*BR
C
      M=2**ITAU
      FAKTOR=1.0D0/DBLE(M)
      PI=4.0D0*ATAN(1.0D0)
      VZ=DBLE(2*IR-1)
      EWPHI=VZ*2.0D0*PI*FAKTOR
      EWR=COS(EWPHI)
      EWI=SIN(EWPHI)
      IF (IR .EQ. 1) FAKTOR=1.0D0
C
C*****************************************************************
C     Restore using bit reversal                                 *
C     ( normalize, if IR=0 )                                     *
C*****************************************************************
C
      DO 30 J = 0, M-1
        K=J
        SIGMA=0
        DO 20 N = 1, ITAU
          KD2=K/2
          SIGMA=2*SIGMA+K-2*KD2
          K=KD2
  20    CONTINUE
        IF (SIGMA .LT. J) GOTO 30
        UR=FRE(J)
        FRE(J)=FRE(SIGMA)*FAKTOR
        FRE(SIGMA)=UR*FAKTOR
        UI=FIM(J)
        FIM(J)=FIM(SIGMA)*FAKTOR
        FIM(SIGMA)=UI*FAKTOR
  30  CONTINUE
C
C*****************************************************************
C     Execute the (inverse) transformation                       *
C*****************************************************************
C
C*****  N min 1 = 2**( N - 1 )         *****
C*****  N min 0 = 2**( N )             *****
C*****     EW   = root of unity        *****
C*****     W    = EW**(2**(ITAU-N))    *****
C*****    EPS   = EW**(L*2**(ITAU-N))  *****
C
      NMIN1=1
      DO 130 N = 1, ITAU
        WR=EWR
        WI=EWI
        DO 100 K = 1, ITAU-N
           H  =CMLR(WR,WI,WR,WI)
           WI =CMLI(WR,WI,WR,WI)
           WR =H
  100   CONTINUE
        EPSR=1.0D0
        EPSI=0.0D0
        NMIN0=NMIN1+NMIN1
        DO 120 L = 0, NMIN1-1
          DO 110 J = 0, M-NMIN0, NMIN0
C
C     U=F(J+L+NMIN1)*EPS
C
            UR=CMLR(FRE(J+L+NMIN1),FIM(J+L+NMIN1),EPSR,EPSI)
            UI=CMLI(FRE(J+L+NMIN1),FIM(J+L+NMIN1),EPSR,EPSI)
            FRE(J+L+NMIN1)=FRE(J+L)-UR
            FIM(J+L+NMIN1)=FIM(J+L)-UI
            FRE(J+L)=FRE(J+L)+UR
            FIM(J+L)=FIM(J+L)+UI
  110     CONTINUE
          H    =CMLR(EPSR,EPSI,WR,WI)
          EPSI =CMLI(EPSR,EPSI,WR,WI)
          EPSR =H
  120   CONTINUE
        NMIN1=NMIN0
  130 CONTINUE
      RETURN
      END
