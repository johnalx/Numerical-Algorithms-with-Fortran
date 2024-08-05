      SUBROUTINE FDICHT (M, FRE, FIM, P, TETA) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This program  computes the values of the trigonometric        *      
!  interpolating polynomial (i.e. the discrete partial Fourier   *      
!  sum) using the Fast Fourier Transform (FFT) for a given set   *      
!  of function values F(0), F(1), ... , F(M-1) of a P-periodic   *      
!  function F for equidistant nodes t(0), t(1), ... , t(M-1),    *      
!  t(j) = j*P/M  at a set of shifted nodes t(j) + TETA  for      *      
!  j = 0, 1, ... , M-1.  ("Increase of number of nodes")         *      
!  The double precision vectors  FRE  and  FIM  contain the real *      
!  and imaginary parts of the values for F.                      *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETER:                                              *      
!  ================                                              *      
!  M        : Number of original nodes.                          *      
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
!  P        : Period interval for F.                             *      
!  TETA     : Shift parameter: the newly computed values         *      
!             correspond to arguments shifted by TETA.           *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  FRE, FIM : DOUBLE PRECISION vectors;                          *      
!             FRE(j) and FIM(j) are the real and imaginary parts *      
!             of  F(j)  of the trigonometric interpolating       *      
!             polynom at the equidistant nodes                   *      
!                 tneu(j) = t(j) + TETA = j*P/M + TETA           *      
!             for   j = 0, 1, ... , M-1 .                        *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Subroutines used: FFT or FFTB                                 *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author      : Klaus Niederdrenk                               *      
!  Date        : 06.30.1994                                      *      
!  Source code : FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      INTEGER M, ITAU 
      DIMENSION FRE (0:M - 1), FIM (0:M - 1) 
      LOGICAL POT2 
!                                                                       
!     Complex multiplication : Real part                                
!                                                                       
      CMLR (AR, AI, BR, BI) = AR * BR - AI * BI 
!                                                                       
!     Complex multiplication : Imaginary part                           
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
!****************************************************************       
!     Determine discrete Fourier transformation                 *       
!****************************************************************       
!                                                                       
      IF (POT2) THEN 
         CALL FFT (ITAU, FRE, FIM, 0) 
      ELSE 
         CALL FFTB (M, FRE, FIM, 0) 
      ENDIF 
!                                                                       
!****************************************************************       
!     Evaluate the discrete Fourier transformation for the      *       
!     shifted nodes                                             *       
!****************************************************************       
!                                                                       
      PI = 4.0D0 * ATAN (1.0D0) 
      FAKTOR = 2.0D0 * PI * TETA / P 
      TETPI = DBLE ( - M / 2) * FAKTOR 
      EKR = COS (TETPI) 
      EKI = SIN (TETPI) 
      HR = CMLR (FRE (M / 2), FIM (M / 2), EKR, EKI) 
      HI = CMLI (FRE (M / 2), FIM (M / 2), EKR, EKI) 
      DO 10 K = - (M - 1) / 2, - 1 
         TETPI = DBLE (K) * FAKTOR 
         EKR = COS (TETPI) 
         EKI = SIN (TETPI) 
         H = CMLR (FRE (K + M), FIM (K + M), EKR, EKI) 
         FIM (K + M) = CMLI (FRE (K + M), FIM (K + M), EKR, EKI) 
         FRE (K + M) = H 
         H = CMLR (FRE ( - K), FIM ( - K), EKR, - EKI) 
         FIM ( - K) = CMLI (FRE ( - K), FIM ( - K), EKR, - EKI) 
         FRE ( - K) = H 
   10 END DO 
      FRE (M / 2) = HR 
      FIM (M / 2) = HI 
!                                                                       
!****************************************************************       
!     Determine functional values                               *       
!****************************************************************       
!                                                                       
      IF (POT2) THEN 
         CALL FFT (ITAU, FRE, FIM, 1) 
      ELSE 
         CALL FFTB (M, FRE, FIM, 1) 
      ENDIF 
      RETURN 
      END SUBROUTINE FDICHT                         
