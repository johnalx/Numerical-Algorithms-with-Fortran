      SUBROUTINE FFAKO (M, FRE, FIM, HRE, HIM) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This program uses the fast Fourier transform (FFT) to compute *      
!  the discrete values of the convolution                        *      
!                                                                *      
!      Falt(j) = 1/M * (SUM k=0 to M-1) F(j-k)*H(k)              *      
!                                                                *      
!  and of the discrete cyclic correlation                        *      
!                                                                *      
!      Korr(j) = 1/M * (SUM k=0 to M-1) F(j+k)*CONJG(H(k))       *      
!                                                                *      
!  of F and H at the given complex functional values             *      
!  F(0), F(1), ... F(M-1) and  H(0), H(1), ... , H(M-1)  at      *      
!  equidistant nodes in the period interval for                  *      
!  j = 0, 1, ... , M-1 .                                         *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  M        : Number of nodes.                                   *      
!             If M is a power of two ( M = 2**ITAU  for a        *      
!             positive integer ITAU ), we make use of the        *      
!             subroutineso  FFT, which runs for all values of M. *      
!             If M  is not a power of two, then we call on FFTB, *      
!             which limits the available size for its auxiliary  *      
!             vectors F1RE, F1IM, GRE, GIM  to  1 <= M <= 1366.  *      
!             For larger M, the dimension statement inside FFTB  *      
!             must be amended as indicated there.                *      
!  FRE, FIM : DOUBLE PRECISION vectors for the real and          *      
!             imaginary parts for M functional values of F.      *      
!  HRE, HIM : DOUBLE PRECISION vectors for the real and          *      
!             imaginary parts for M functional values of H.      *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  FRE, FIM : DOUBLE PRECISION vectors with real and imaginary   *      
!             parts of the discrete cyclic convolution Falt(j).  *      
!  HRE, HIM : DOUBLE PRECISION vectors with real and imaginary   *      
!             parts of the discrete cyclic correlation Korr(j).  *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Subroutines used: FFT or FFTB                                 *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author      : Klaus Niederdrenk                               *      
!  Date        : 09.08.1994                                      *      
!  Source code : FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      INTEGER M, ITAU 
      DIMENSION FRE (0:M - 1), FIM (0:M - 1), HRE (0:M - 1), HIM (0:M - &
      1)                                                                
      LOGICAL POT2 
!                                                                       
!     complex multiplication: Real part                                 
!                                                                       
      CMLR (AR, AI, BR, BI) = AR * BR - AI * BI 
!                                                                       
!     complex multiplication: Imaginary part                            
!                                                                       
      CMLI (AR, AI, BR, BI) = AR * BI + AI * BR 
!                                                                       
!     Check whether M is a power of two                                 
!                                                                       
      ITAU = NINT (LOG (REAL (M) ) / LOG (2.0) ) 
      IF (2**ITAU.EQ.M) THEN 
         POT2 = .TRUE. 
      ELSE 
         POT2 = .FALSE. 
      ENDIF 
!                                                                       
!*****************************************************************      
!     Compute the discrete Fourier transformation                *      
!*****************************************************************      
!                                                                       
      IF (POT2) THEN 
         CALL FFT (ITAU, FRE, FIM, 0) 
         CALL FFT (ITAU, HRE, HIM, 0) 
      ELSE 
         CALL FFTB (M, FRE, FIM, 0) 
         CALL FFTB (M, HRE, HIM, 0) 
      ENDIF 
!                                                                       
!*****************************************************************      
!     Conpute the transform of the corresponding discrete cyclic *      
!     convolution and of the discrete cyclic correlation         *      
!*****************************************************************      
!                                                                       
      DO 10 K = 0, M - 1 
         H1 = CMLR (FRE (K), FIM (K), HRE (K), HIM (K) ) 
         H2 = CMLI (FRE (K), FIM (K), HRE (K), HIM (K) ) 
         H = CMLR (FRE (K), FIM (K), HRE (K), - HIM (K) ) 
         HIM (K) = CMLI (FRE (K), FIM (K), HRE (K), - HIM (K) ) 
         HRE (K) = H 
         FRE (K) = H1 
         FIM (K) = H2 
   10 END DO 
!                                                                       
!*****************************************************************      
!     Evaluate the functional values                             *      
!*****************************************************************      
!                                                                       
      IF (POT2) THEN 
         CALL FFT (ITAU, FRE, FIM, 1) 
         CALL FFT (ITAU, HRE, HIM, 1) 
      ELSE 
         CALL FFTB (M, FRE, FIM, 1) 
         CALL FFTB (M, HRE, HIM, 1) 
      ENDIF 
      RETURN 
      END SUBROUTINE FFAKO                          
