C[BA*)
C[LE*)
      SUBROUTINE FFTB (M, FRE, FIM, IR)
C[IX{FFTB}*)
C
C*****************************************************************
C                                                                *
C  For IR = 0 this programm uses the fast Fourier transform      *
C  ( F F T ) in order to determine the discrete Fourier          *
C  coefficients  F^(-M/2), ... , F^(M/2-1) for an arbitrary      *
C  number M of given real or complex functional values F(0),     *
C  F(1), ... , F(M-1).                                           *
C[BE*)
C  These coefficients define the discrete Fourier series         *
C                                                                *
C     (SUM from K=-M/2 to M/2-1)  F^(K)*EXP(I*K*OMEGA*X)         *
C                                                                *
C  if M is even, or                                              *
C                                                                *
C     (SUM from K=-(M-1)/2 to (M-1)/2)  F^(K)*EXP(I*K*OMEGA*X) , *
C                                                                *
C  if  M  is odd.                                                *
C                                                                *
C  ( I : imaginary unit with  I**2 = -1 ; OMEGA = 2*PI/L ,       *
C    L : period ) .                                              *
C  For  IR = 1  the inverse Fourier transform is perfomed.       *
C  The double precision vectors  FRE  and  FIM  contain the real *
C  and imaginary parts of the values  F  and  F^.                *
C  The actual computations are performed using an FFT with radix *
C  2 and a discrete convolution.                                 *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  M      : the number of data points                            *
C           Due to the explicit declarations of the auxiliary    *
C           vectors F1RE, F1IM and GRE, GIM in the program, we   *
C           can only allow M  with 1 <= M <= 1366.               *
C           Larger values for M require that redimensioning of   *
C           F1RE, F1IM and GRe, GIM to admit 0 : 2**TAU-1        *
C           entries, so that  M <= (2**TAU + 1)/3.               *
C           For example:                                         *
C           for  M <=  2731 use  indices between 0  and  8191,   *
C           for  M <=  5462 use  indices between 0  and 16383,   *
C           for  M <= 10923 use  indices between 0  and 32767,   *
C           for  M <= 21846 use  indices between 0  and 65535,   *
C           etc.                                                 *
C  FRE, FIM : DOUBLE PRECISION vectors for the real and imaginary*
C             parts of M complex numbers:                        *
C             Depending on IR, FRE, FIM  must contain the        *
C             following:                                         *
C             IR=0: FRE, FIM  contain the real and imaginary     *
C                   parts of the functional values               *
C             IR=1: FRE, FIM  contain the real and imaginary     *
C                   parts of the discrete Fourier coefficients   *
C                   as follows:                                  *
C                   FRE(K), FIM(K) : Real and imaginary parts of *
C                   F^(K) for K = 0 to M/2-1   if M is even, and *
C                         for K = 0 to (M-1)/2 if M is odd,      *
C                   and of                                       *
C                   F^(K-M) for K = M/2 to M-1 if M is even, and *
C                           for K = (M+1)/2 to M-1 if M is odd.  *
C  IR     : determines the direction of the desired              *
C           transfornation:                                      *
C           If IR=0, we compute the discrete Fourier             *
C                    coefficients                                *
C           If IR=1: we compute the function values              *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  FRE, FIM : DOUBLE PRECISION vectors of real and imaginary     *
C             parts of M complex numbers.                        *
C             Depending on IR, these contain the following:      *
C             If IR=0: They contain the real and imaginary parts *
C                      of the discrete Fourier coefficients as   *
C                      follows:                                  *
C                FRE(K+M), FIM(K+M) : real and imaginary parts of*
C                F^(K) for K = -M/2, ..., -1    if M is even, or *
C                      for K = -(M-1)/2,..., -1 if M is odd;     *
C                F^(K) for K = 0, ..., M/2-1   if M is even, or  *
C                      for K = 0, ..., (M-1)/2 if M is odd.      *
C             If IR=1: the function values are returned.         *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Required subroutines: FFT                                     *
C                                                                *
C*****************************************************************
C                                                                *
C  Authors   : Klaus Niederdrenk                                 *
C  Date      : 03.22.1992                                        *
C  Source    : FORTRAN 77                                        *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER M, IR, ITAU
      DIMENSION FRE(0:M-1), FIM(0:M-1)
      DIMENSION F1RE(0:4095), F1IM(0:4095)
      DIMENSION GRE(0:4095), GIM(0:4095)
C
C     Define FUNCTION for complex multiplication: real part
C
      CMLR(AR,AI,BR,BI)=AR*BR-AI*BI
C
C     Define FUNCTION for complex multiplication: imaginary part
C
      CMLI(AR,AI,BR,BI)=AR*BI+AI*BR
C
      FAKTOR=1.0D0/DBLE(M)
      PI=4.0D0*ATAN(1.0D0)
      VZ=DBLE(2*IR-1)
C
C     EW1=EXP(CMPLX(0.0D0,VZ*PI*FAKTOR))
C
      EWPHI=VZ*PI*FAKTOR
      EW1R=COS(EWPHI)
      EW1I=SIN(EWPHI)
C
C     EW2=EW1*EW1
C
      EW2R=CMLR(EW1R,EW1I,EW1R,EW1I)
      EW2I=CMLI(EW1R,EW1I,EW1R,EW1I)
      IF(IR .EQ. 1) FAKTOR=1.0D0
C
C*****************************************************************
C     Determine the suitable power of two that shall become the  *
C     length of the auxiliary vectors F1RE, F1IM  and  GRE, GIM  *
C*****************************************************************
C
      ITAU=INT(LOG(3.0D0*DBLE(M)-2.0D0)/LOG(2.0D0))+1
      L=2**ITAU
      IF(L/2 .GE. 3*M-2) THEN
        L=L/2
        ITAU=ITAU-1
      ENDIF
C
C*****************************************************************
C     Initialize the vectors  F1RE, F1IM and  GRE, GIM           *
C*****************************************************************
C
C*****   EW1 = 2*M - th root of unity    *****
C*****   EWK = EW1 ** (J**2)             *****
C
      DO 10 J = 0, L-1
        F1RE(J)=0.0D0
        F1IM(J)=0.0D0
        GRE(J)=0.0D0
        GIM(J)=0.0D0
  10  CONTINUE
      F1RE(0)=FRE(0)
      F1IM(0)=FIM(0)
      GRE(M-1)=1.0D0
C
C     EWK=EW1
C
      EWKR=EW1R
      EWKI=EW1I
C
C     EW3=EW1
C
      EW3R=EW1R
      EW3I=EW1I
      DO 20 J = 1, M-1
C
C       F1(J)=F(J)*EWK
C
        F1RE(J)=CMLR(FRE(J),FIM(J),EWKR,EWKI)
        F1IM(J)=CMLI(FRE(J),FIM(J),EWKR,EWKI)
C
C       G(M-1+J)=CONJG(EWK)
C
        GRE(M-1+J)=EWKR
        GIM(M-1+J)=-EWKI
C
C       G(M-1-J)=G(M-1+J)
C
        GRE(M-1-J)=GRE(M-1+J)
        GIM(M-1-J)=GIM(M-1+J)
C
C       EW3=EW3*EW2
C
        H   =CMLR(EW3R,EW3I,EW2R,EW2I)
        EW3I=CMLI(EW3R,EW3I,EW2R,EW2I)
        EW3R=H
C
C       EWK=EWK*EW3
C
        H   =CMLR(EWKR,EWKI,EW3R,EW3I)
        EWKI=CMLI(EWKR,EWKI,EW3R,EW3I)
        EWKR=H
  20  CONTINUE
C
C*****************************************************************
C     Discrete convolution of the vectors F1RE, F1IM and GRE, GIM*
C     using the FFT for radix two (via subroutine FFT)           *
C*****************************************************************
C
      CALL FFT (ITAU, F1RE, F1IM, 0)
      CALL FFT (ITAU, GRE, GIM, 0)
      DO 30 K = 0, L-1
C
C       F1(K)=F1(K)*G(K)
C
        H      =CMLR(F1RE(K),F1IM(K),GRE(K),GIM(K))
        F1IM(K)=CMLI(F1RE(K),F1IM(K),GRE(K),GIM(K))
        F1RE(K)=H
  30  CONTINUE
      CALL FFT (ITAU, F1RE, F1IM, 1)
C
C*****************************************************************
C     Store needed values in  FRE, FIM                           *
C*****************************************************************
C
C*****   EW1 = 2*M - th root of unity    *****
C*****   EWK = EW1 ** (K**2)             *****
C
      FAKTL=FAKTOR*DBLE(L)
      FRE(0)=F1RE(M-1)*FAKTL
      FIM(0)=F1IM(M-1)*FAKTL
      EWKR=EW1R
      EWKI=EW1I
      EW3R=EW1R
      EW3I=EW1I
      DO 40 K = 1, M-1
C
C       F(K)=F1(K+M-1)*EWK*FAKTL
C
        FRE(K)=CMLR(F1RE(K+M-1),F1IM(K+M-1),EWKR*FAKTL,EWKI*FAKTL)
        FIM(K)=CMLI(F1RE(K+M-1),F1IM(K+M-1),EWKR*FAKTL,EWKI*FAKTL)
C
C       EW3=EW3*EW2
C
        H   =CMLR(EW3R,EW3I,EW2R,EW2I)
        EW3I=CMLI(EW3R,EW3I,EW2R,EW2I)
        EW3R=H
C
C       EWK=EWK*EW3
C
        H   =CMLR(EWKR,EWKI,EW3R,EW3I)
        EWKI=CMLI(EWKR,EWKI,EW3R,EW3I)
        EWKR=H
  40  CONTINUE
      RETURN
      END
