      SUBROUTINE FFT (ITAU, FRE, FIM, IR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  If IR = 0, this programm uses the fast Fourier transform      *      
!  ( F F T ) in order to determine the discrete Fourier          *      
!  coefficients   F^(-M/2), ... , F^(M/2-1)  for  M = 2**ITAU    *      
!  given real or complex functional values F(0), F(1), ... ,     *      
!  F(M-1).                                                       *      
!  These coefficients define the discrete Fourier series         *      
!                                                                *      
!     (SUM from K=-M/2 to M/2-1)  F^(K)*EXP(I*K*OMEGA*X)         *      
!                                                                *      
!  ( I : imaginary unit with  I**2 = -1 ; OMEGA = 2*PI/L ,       *      
!    L : period ) .                                              *      
!  For  IR = 1  the inverse Fourier transform is performed.      *      
!  The double precision vectors  FRE  and  FIM  contain the real *      
!  and imaginary parts of  F  and  F^, respectively.             *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  ITAU   : the number of data points is  M = 2**ITAU            *      
!  FRE, FIM : DOUBLE PRECISION vectors for the real and imaginary*      
!             parts of M values.                                 *      
!             Depending on IR, FRE, FIM  must contain the        *      
!             following:                                         *      
!             IR=0: FRE, FIM  contain the real and imaginary     *      
!                   parts of the functional values               *      
!             IR=1: FRE, FIM  contain the real and imaginary     *      
!                   parts of the discrete Fourier coefficients   *      
!                   as follows:                                  *      
!                   FRE(K), FIM(K) : Real and imaginary parts of *      
!                   F^(K)  for  K = 0 to M/2-1  and of F^(K-M)   *      
!                   for  K = M/2 to M-1.                         *      
!  IR     : determines the direction of the desired              *      
!           transfornation                                       *      
!           If IR=0, we compute the discrete Fourier             *      
!                    coefficients                                *      
!           If IR=1: we compute the function values              *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  FRE, FIM : DOUBLE PRECISION vectors for the real and          *      
!             imaginary parts of M complex numbers.              *      
!             Depending on IR the following values are returned: *      
!             If IR=0: The vectors FRE and FIM contain the       *      
!                      real and imaginary parts of the discrete  *      
!                      Fourier coefficients as follows:          *      
!                      FRE(K+M), FIM(K+M): The real and imaginary*      
!                           parts of  F^(K) for  K = -M/2 to -1  *      
!                      FRE(K), FIM(K) : The real and imaginary   *      
!                           parts of  F^(K)  for  K = 0 to M/2-1 *      
!             If IR=1: The functional values are returned in FRE *      
!                      and FIM.                                  *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: none                                    *      
!                                                                *      
!  Reference: Niederdrenk, K., see [NIED84]                      *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author    : Klaus Niederdrenk                                 *      
!  Date      : 03.22.1993                                        *      
!  Source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      INTEGER ITAU, IR, SIGMA, M 
      DIMENSION FRE (0:2**ITAU - 1), FIM (0:2**ITAU - 1) 
!                                                                       
!     Define FUNCTION for complex multiplication: real part             
!                                                                       
      CMLR (AR, AI, BR, BI) = AR * BR - AI * BI 
!                                                                       
!     Define FUNCTION for complex multiplication: real part             
!                                                                       
      CMLI (AR, AI, BR, BI) = AR * BI + AI * BR 
!                                                                       
      M = 2**ITAU 
      FAKTOR = 1.0D0 / DBLE (M) 
      PI = 4.0D0 * ATAN (1.0D0) 
      VZ = DBLE (2 * IR - 1) 
      EWPHI = VZ * 2.0D0 * PI * FAKTOR 
      EWR = COS (EWPHI) 
      EWI = SIN (EWPHI) 
      IF (IR.EQ.1) FAKTOR = 1.0D0 
!                                                                       
!*****************************************************************      
!     Restore using bit reversal                                 *      
!     ( normalize, if IR=0 )                                     *      
!*****************************************************************      
!                                                                       
      DO 30 J = 0, M - 1 
         K = J 
         SIGMA = 0 
         DO 20 N = 1, ITAU 
            KD2 = K / 2 
            SIGMA = 2 * SIGMA + K - 2 * KD2 
            K = KD2 
   20    END DO 
         IF (SIGMA.LT.J) GOTO 30 
         UR = FRE (J) 
         FRE (J) = FRE (SIGMA) * FAKTOR 
         FRE (SIGMA) = UR * FAKTOR 
         UI = FIM (J) 
         FIM (J) = FIM (SIGMA) * FAKTOR 
         FIM (SIGMA) = UI * FAKTOR 
   30 END DO 
!                                                                       
!*****************************************************************      
!     Execute the (inverse) transformation                       *      
!*****************************************************************      
!                                                                       
!*****  N min 1 = 2**( N - 1 )         *****                            
!*****  N min 0 = 2**( N )             *****                            
!*****     EW   = root of unity        *****                            
!*****     W    = EW**(2**(ITAU-N))    *****                            
!*****    EPS   = EW**(L*2**(ITAU-N))  *****                            
!                                                                       
      NMIN1 = 1 
      DO 130 N = 1, ITAU 
         WR = EWR 
         WI = EWI 
         DO 100 K = 1, ITAU - N 
            H = CMLR (WR, WI, WR, WI) 
            WI = CMLI (WR, WI, WR, WI) 
            WR = H 
  100    END DO 
         EPSR = 1.0D0 
         EPSI = 0.0D0 
         NMIN0 = NMIN1 + NMIN1 
         DO 120 L = 0, NMIN1 - 1 
            DO 110 J = 0, M - NMIN0, NMIN0 
!                                                                       
!     U=F(J+L+NMIN1)*EPS                                                
!                                                                       
               UR = CMLR (FRE (J + L + NMIN1), FIM (J + L + NMIN1),     &
               EPSR, EPSI)                                              
               UI = CMLI (FRE (J + L + NMIN1), FIM (J + L + NMIN1),     &
               EPSR, EPSI)                                              
               FRE (J + L + NMIN1) = FRE (J + L) - UR 
               FIM (J + L + NMIN1) = FIM (J + L) - UI 
               FRE (J + L) = FRE (J + L) + UR 
               FIM (J + L) = FIM (J + L) + UI 
  110       END DO 
            H = CMLR (EPSR, EPSI, WR, WI) 
            EPSI = CMLI (EPSR, EPSI, WR, WI) 
            EPSR = H 
  120    END DO 
         NMIN1 = NMIN0 
  130 END DO 
      RETURN 
      END SUBROUTINE FFT                            
