      SUBROUTINE FOURN (M, FRE, FIM, A, DELTAX) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This program uses the fast Fourier transform (FFT) to compute *      
!  approximate values of the Fourier transform                   *      
!                                                                *      
!       F^(tj) = (INTEGRAL of) F(x) * EXP(-I*tj*x) dx            *      
!                                                                *      
!  ( I : imaginary unit; I**2 = -1 ) for j = -M/2 , ... , M/2 -1.*      
!  Here the nonperiodic function  F is known by its functional   *      
!  values F(0), F(1), ... , F(M-1) at equidistant nodes with     *      
!  uniform distance DELTAX. These nodes lie in the support of F. *      
!  The double precision vectors  FRE  and  FIM  contain the real *      
!  and imaginary parts of the function values of F.              *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  M        : Number of nodes; M must be even.                   *      
!             If M is a power of two ( M = 2**ITAU  for a        *      
!             positive integer ITAU ), we make use of the        *      
!             subroutineso  FFT, which runs for all values of M. *      
!             If M  is not a power of two, then we call on FFTB, *      
!             which limits the available size for its auxiliary  *      
!             vectors F1RE, F1IM, GRE, GIM  to  1 <= M <= 1366.  *      
!             For larger M, the dimension statement inside FFTB  *      
!             must be amended as indicated there.                *      
!  FRE, FIM : DOUBLE PRECISION vectors for the real and          *      
!             imaginary parts for M functional values of F at the*      
!             nodes  x(j) = A + j * DELTAX , j = 0, 1, ... , M-1.*      
!  A        : starting point of the functional values.           *      
!  DELTAX   : uniform distance of the nodes.                     *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  FRE, FIM : DOUBLE PRECISION vectors for real and imaginary    *      
!             parts of M complex numbers.                        *      
!             FRE and FIM contain the values of the Fourier      *      
!             transform  F^  as follows:                         *      
!               FRE(k+M) , FIM(k+M) : Real and imaginary parts of*      
!                               F^(tk)  for  k = -M/2  to -1 ,   *      
!               FRE(k) , FIM(k) : Real and imaginary parts of    *      
!                               F^(tk)  for  k = 0  to  M/2 -1.  *      
!             Here tk = k / (M*DELTAX).                          *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Subroutines used : FFT bzw. FFTB                              *      
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
!     complex multiplication: Real part                                 
!                                                                       
      CMLR (AR, AI, BR, BI) = AR * BR - AI * BI 
!                                                                       
!     complex multiplication: Imaginary part                            
!                                                                       
      CMLI (AR, AI, BR, BI) = AR * BI + AI * BR 
!                                                                       
!     Determine whether M is a power of two                             
!                                                                       
      ITAU = NINT (LOG (REAL (M) ) / LOG (2.0) ) 
      IF (2**ITAU.EQ.M) THEN 
         POT2 = .TRUE. 
      ELSE 
         POT2 = .FALSE. 
      ENDIF 
!                                                                       
!*****************************************************************      
!     Determine discrete Fourier coeffizients                    *      
!*****************************************************************      
!                                                                       
      IF (POT2) THEN 
         CALL FFT (ITAU, FRE, FIM, 0) 
      ELSE 
         CALL FFTB (M, FRE, FIM, 0) 
      ENDIF 
!                                                                       
!*****************************************************************      
!     Transform to the nonperiodic case:                         *      
!     Adjust values to those of the Fourier transform            *      
!*****************************************************************      
!                                                                       
      X = DBLE (M) * DELTAX 
      PI = 4.0D0 * ATAN (1.0D0) 
      FAKTOR = 2.0D0 * PI * A / X 
      DO 10 K = - M / 2 + 1, - 1 
         ARG = DBLE (K) * FAKTOR 
         EKR = X * COS (ARG) 
         EKI = X * SIN (ARG) 
         H = CMLR (FRE (K + M), FIM (K + M), EKR, - EKI) 
         FIM (K + M) = CMLI (FRE (K + M), FIM (K + M), EKR, - EKI) 
         FRE (K + M) = H 
         H = CMLR (FRE ( - K), FIM ( - K), EKR, EKI) 
         FIM ( - K) = CMLI (FRE ( - K), FIM ( - K), EKR, EKI) 
         FRE ( - K) = H 
   10 END DO 
      FRE (0) = X * FRE (0) 
      FIM (0) = X * FIM (0) 
      ARG = DBLE (M / 2) * FAKTOR 
      EKR = X * COS (ARG) 
      EKI = X * SIN (ARG) 
      H = CMLR (FRE (M / 2), FIM (M / 2), EKR, EKI) 
      FIM (M / 2) = CMLI (FRE (M / 2), FIM (M / 2), EKR, EKI) 
      FRE (M / 2) = H 
      RETURN 
      END SUBROUTINE FOURN                          
