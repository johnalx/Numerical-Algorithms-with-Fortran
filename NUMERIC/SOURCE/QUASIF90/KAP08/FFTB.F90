      SUBROUTINE FFTB (M, FRE, FIM, IR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  For IR = 0 this programm uses the fast Fourier transform      *      
!  ( F F T ) in order to determine the discrete Fourier          *      
!  coefficients  F^(-M/2), ... , F^(M/2-1) for an arbitrary      *      
!  number M of given real or complex functional values F(0),     *      
!  F(1), ... , F(M-1).                                           *      
!  These coefficients define the discrete Fourier series         *      
!                                                                *      
!     (SUM from K=-M/2 to M/2-1)  F^(K)*EXP(I*K*OMEGA*X)         *      
!                                                                *      
!  if M is even, or                                              *      
!                                                                *      
!     (SUM from K=-(M-1)/2 to (M-1)/2)  F^(K)*EXP(I*K*OMEGA*X) , *      
!                                                                *      
!  if  M  is odd.                                                *      
!                                                                *      
!  ( I : imaginary unit with  I**2 = -1 ; OMEGA = 2*PI/L ,       *      
!    L : period ) .                                              *      
!  For  IR = 1  the inverse Fourier transform is perfomed.       *      
!  The double precision vectors  FRE  and  FIM  contain the real *      
!  and imaginary parts of the values  F  and  F^.                *      
!  The actual computations are performed using an FFT with radix *      
!  2 and a discrete convolution.                                 *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  M      : the number of data points                            *      
!           Due to the explicit declarations of the auxiliary    *      
!           vectors F1RE, F1IM and GRE, GIM in the program, we   *      
!           can only allow M  with 1 <= M <= 1366.               *      
!           Larger values for M require that redimensioning of   *      
!           F1RE, F1IM and GRe, GIM to admit 0 : 2**TAU-1        *      
!           entries, so that  M <= (2**TAU + 1)/3.               *      
!           For example:                                         *      
!           for  M <=  2731 use  indices between 0  and  8191,   *      
!           for  M <=  5462 use  indices between 0  and 16383,   *      
!           for  M <= 10923 use  indices between 0  and 32767,   *      
!           for  M <= 21846 use  indices between 0  and 65535,   *      
!           etc.                                                 *      
!  FRE, FIM : DOUBLE PRECISION vectors for the real and imaginary*      
!             parts of M complex numbers:                        *      
!             Depending on IR, FRE, FIM  must contain the        *      
!             following:                                         *      
!             IR=0: FRE, FIM  contain the real and imaginary     *      
!                   parts of the functional values               *      
!             IR=1: FRE, FIM  contain the real and imaginary     *      
!                   parts of the discrete Fourier coefficients   *      
!                   as follows:                                  *      
!                   FRE(K), FIM(K) : Real and imaginary parts of *      
!                   F^(K) for K = 0 to M/2-1   if M is even, and *      
!                         for K = 0 to (M-1)/2 if M is odd,      *      
!                   and of                                       *      
!                   F^(K-M) for K = M/2 to M-1 if M is even, and *      
!                           for K = (M+1)/2 to M-1 if M is odd.  *      
!  IR     : determines the direction of the desired              *      
!           transfornation:                                      *      
!           If IR=0, we compute the discrete Fourier             *      
!                    coefficients                                *      
!           If IR=1: we compute the function values              *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  FRE, FIM : DOUBLE PRECISION vectors of real and imaginary     *      
!             parts of M complex numbers.                        *      
!             Depending on IR, these contain the following:      *      
!             If IR=0: They contain the real and imaginary parts *      
!                      of the discrete Fourier coefficients as   *      
!                      follows:                                  *      
!                FRE(K+M), FIM(K+M) : real and imaginary parts of*      
!                F^(K) for K = -M/2, ..., -1    if M is even, or *      
!                      for K = -(M-1)/2,..., -1 if M is odd;     *      
!                F^(K) for K = 0, ..., M/2-1   if M is even, or  *      
!                      for K = 0, ..., (M-1)/2 if M is odd.      *      
!             If IR=1: the function values are returned.         *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: FFT                                     *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Authors   : Klaus Niederdrenk                                 *      
!  Date      : 03.22.1992                                        *      
!  Source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      INTEGER M, IR, ITAU 
      DIMENSION FRE (0:M - 1), FIM (0:M - 1) 
      DIMENSION F1RE (0:4095), F1IM (0:4095) 
      DIMENSION GRE (0:4095), GIM (0:4095) 
!                                                                       
!     Define FUNCTION for complex multiplication: real part             
!                                                                       
      CMLR (AR, AI, BR, BI) = AR * BR - AI * BI 
!                                                                       
!     Define FUNCTION for complex multiplication: imaginary part        
!                                                                       
      CMLI (AR, AI, BR, BI) = AR * BI + AI * BR 
!                                                                       
      FAKTOR = 1.0D0 / DBLE (M) 
      PI = 4.0D0 * ATAN (1.0D0) 
      VZ = DBLE (2 * IR - 1) 
!                                                                       
!     EW1=EXP(CMPLX(0.0D0,VZ*PI*FAKTOR))                                
!                                                                       
      EWPHI = VZ * PI * FAKTOR 
      EW1R = COS (EWPHI) 
      EW1I = SIN (EWPHI) 
!                                                                       
!     EW2=EW1*EW1                                                       
!                                                                       
      EW2R = CMLR (EW1R, EW1I, EW1R, EW1I) 
      EW2I = CMLI (EW1R, EW1I, EW1R, EW1I) 
      IF (IR.EQ.1) FAKTOR = 1.0D0 
!                                                                       
!*****************************************************************      
!     Determine the suitable power of two that shall become the  *      
!     length of the auxiliary vectors F1RE, F1IM  and  GRE, GIM  *      
!*****************************************************************      
!                                                                       
      ITAU = INT (LOG (3.0D0 * DBLE (M) - 2.0D0) / LOG (2.0D0) )        &
      + 1                                                               
      L = 2**ITAU 
      IF (L / 2.GE.3 * M - 2) THEN 
         L = L / 2 
         ITAU = ITAU - 1 
      ENDIF 
!                                                                       
!*****************************************************************      
!     Initialize the vectors  F1RE, F1IM and  GRE, GIM           *      
!*****************************************************************      
!                                                                       
!*****   EW1 = 2*M - th root of unity    *****                          
!*****   EWK = EW1 ** (J**2)             *****                          
!                                                                       
      DO 10 J = 0, L - 1 
         F1RE (J) = 0.0D0 
         F1IM (J) = 0.0D0 
         GRE (J) = 0.0D0 
         GIM (J) = 0.0D0 
   10 END DO 
      F1RE (0) = FRE (0) 
      F1IM (0) = FIM (0) 
      GRE (M - 1) = 1.0D0 
!                                                                       
!     EWK=EW1                                                           
!                                                                       
      EWKR = EW1R 
      EWKI = EW1I 
!                                                                       
!     EW3=EW1                                                           
!                                                                       
      EW3R = EW1R 
      EW3I = EW1I 
      DO 20 J = 1, M - 1 
!                                                                       
!       F1(J)=F(J)*EWK                                                  
!                                                                       
         F1RE (J) = CMLR (FRE (J), FIM (J), EWKR, EWKI) 
         F1IM (J) = CMLI (FRE (J), FIM (J), EWKR, EWKI) 
!                                                                       
!       G(M-1+J)=CONJG(EWK)                                             
!                                                                       
         GRE (M - 1 + J) = EWKR 
         GIM (M - 1 + J) = - EWKI 
!                                                                       
!       G(M-1-J)=G(M-1+J)                                               
!                                                                       
         GRE (M - 1 - J) = GRE (M - 1 + J) 
         GIM (M - 1 - J) = GIM (M - 1 + J) 
!                                                                       
!       EW3=EW3*EW2                                                     
!                                                                       
         H = CMLR (EW3R, EW3I, EW2R, EW2I) 
         EW3I = CMLI (EW3R, EW3I, EW2R, EW2I) 
         EW3R = H 
!                                                                       
!       EWK=EWK*EW3                                                     
!                                                                       
         H = CMLR (EWKR, EWKI, EW3R, EW3I) 
         EWKI = CMLI (EWKR, EWKI, EW3R, EW3I) 
         EWKR = H 
   20 END DO 
!                                                                       
!*****************************************************************      
!     Discrete convolution of the vectors F1RE, F1IM and GRE, GIM*      
!     using the FFT for radix two (via subroutine FFT)           *      
!*****************************************************************      
!                                                                       
      CALL FFT (ITAU, F1RE, F1IM, 0) 
      CALL FFT (ITAU, GRE, GIM, 0) 
      DO 30 K = 0, L - 1 
!                                                                       
!       F1(K)=F1(K)*G(K)                                                
!                                                                       
         H = CMLR (F1RE (K), F1IM (K), GRE (K), GIM (K) ) 
         F1IM (K) = CMLI (F1RE (K), F1IM (K), GRE (K), GIM (K) ) 
         F1RE (K) = H 
   30 END DO 
      CALL FFT (ITAU, F1RE, F1IM, 1) 
!                                                                       
!*****************************************************************      
!     Store needed values in  FRE, FIM                           *      
!*****************************************************************      
!                                                                       
!*****   EW1 = 2*M - th root of unity    *****                          
!*****   EWK = EW1 ** (K**2)             *****                          
!                                                                       
      FAKTL = FAKTOR * DBLE (L) 
      FRE (0) = F1RE (M - 1) * FAKTL 
      FIM (0) = F1IM (M - 1) * FAKTL 
      EWKR = EW1R 
      EWKI = EW1I 
      EW3R = EW1R 
      EW3I = EW1I 
      DO 40 K = 1, M - 1 
!                                                                       
!       F(K)=F1(K+M-1)*EWK*FAKTL                                        
!                                                                       
         FRE (K) = CMLR (F1RE (K + M - 1), F1IM (K + M - 1), EWKR *     &
         FAKTL, EWKI * FAKTL)                                           
         FIM (K) = CMLI (F1RE (K + M - 1), F1IM (K + M - 1), EWKR *     &
         FAKTL, EWKI * FAKTL)                                           
!                                                                       
!       EW3=EW3*EW2                                                     
!                                                                       
         H = CMLR (EW3R, EW3I, EW2R, EW2I) 
         EW3I = CMLI (EW3R, EW3I, EW2R, EW2I) 
         EW3R = H 
!                                                                       
!       EWK=EWK*EW3                                                     
!                                                                       
         H = CMLR (EWKR, EWKI, EW3R, EW3I) 
         EWKI = CMLI (EWKR, EWKI, EW3R, EW3I) 
         EWKR = H 
   40 END DO 
      RETURN 
      END SUBROUTINE FFTB                           
