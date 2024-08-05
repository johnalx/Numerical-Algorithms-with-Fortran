      SUBROUTINE FFAKON (M, FRE, FIM, N, HRE, HIM, L, DELTAX, IND, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This program uses the fast Fourier transform (FFT) to compute *      
!  approximations for the convolution                            *      
!                                                                *      
!      Falt(j) = (INTEGRAL of) F(xj-t) * H(t) dt                 *      
!                                                                *      
!  or approximations for the correlation                         *      
!                                                                *      
!      Korr(j) = (INTEGRAL of) F(xj+t)*CONJG(H(t)) dt            *      
!                                                                *      
!  for  j = 0, 1, ... , M+N .                                    *      
!  Here F and H are two nonperiodic functions given by their     *      
!  funktional values  F(0), F(1), ... F(M) and  H(0), H(1), ... ,*      
!  H(N)  at equidistant nodes with uniform distance DELTAX.      *      
!  These nodes lie in the support of F and G.                    *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  M        : M+1  functional values of F are given with real    *      
!             and imaginary parts  FRE(j)  and  FIM(j)  for      *      
!             j = 0, 1, ... , M.                                 *      
!  FRE, FIM : DOUBLE PRECISION vectors of length L . The first   *      
!             M+1 entries contain the real and imaginary parts   *      
!             of the functional values of F.                     *      
!  N        : N+1 functional values of H are given by their      *      
!             real and imaginary parts  HRE(j)  and  HIM(j)  for *      
!             j = 0, 1, ... , N.                                 *      
!  HRE, HIM : DOUBLE PRECISION vectors of length L . The first   *      
!             N+1 entries contain the real and imaginary parts   *      
!             of the function values of H .                      *      
!  L        : Length of the vectors  FRE, FIM andnd  HRE, HIM.   *      
!             L must be a power of 2 ( L = 2**ITAU for a positive*      
!             integer ITAU ) and  L >= M+N+1 .                   *      
!  DELTAX   : distance of nodes of the given function values. .  *      
!  IND      : parameter to create approximations of the          *      
!             convolution or the correlation as desired:         *      
!             IND = 0              - Compute approximations for  *      
!                                    the convolution.            *      
!             IND different from 0 - Compute approximations for  *      
!                                    the correlation.            *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  FRE, FIM : DOUBLE PRECISION vectors for the real and imaginary*      
!             parts of  M+N  complex numbers.  The contents of   *      
!             FRE(j) and FIM(j) depends on the setting of IND and*      
!             either the approximate values of the convolution   *      
!             Falt(j) or of the correlation Korr(j) are given    *      
!             for j = 0, 1, ... , M+N .                          *      
!  HRE, HIM : DOUBLE PRECISION vectors for the real and imaginary*      
!             of comlpex numbers. Contain the discrete Fourier   *      
!             coefficients of H as they are needed for           *      
!             convoluting or correlating. (Auxiliary arrays.)    *      
!  IERR     : Error parameter                                    *      
!             = 0 : all o.k.                                     *      
!             = 1 : L is improper.                               *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Subroutines used: FFT                                         *      
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
      INTEGER M, N, L, IND, IERR, ITAU 
      DIMENSION FRE (0:L - 1), FIM (0:L - 1), HRE (0:L - 1), HIM (0:L - &
      1)                                                                
!                                                                       
!     complex multiplication: Real part                                 
!                                                                       
      CMLR (AR, AI, BR, BI) = AR * BR - AI * BI 
!                                                                       
!     complex multiplication: Imaginary part                            
!                                                                       
      CMLI (AR, AI, BR, BI) = AR * BI + AI * BR 
!                                                                       
!     Check input parameter                                             
!                                                                       
      ITAU = NINT (LOG (REAL (L) ) / LOG (2.0) ) 
      IF ( (2**ITAU.NE.L) .OR. (L.LE.M + N) ) THEN 
         IERR = 1 
         RETURN 
      ELSE 
         IERR = 0 
      ENDIF 
!                                                                       
!*****************************************************************      
!     Preassign values to elongated vectors                      *      
!*****************************************************************      
!                                                                       
      IF (IND.EQ.0) THEN 
         DO 10 J = M + 1, L - 1 
            FRE (J) = 0.0D0 
            FIM (J) = 0.0D0 
   10    END DO 
      ELSE 
         DO 20 J = M + N + 1, L - 1 
            FRE (J) = 0.0D0 
            FIM (J) = 0.0D0 
   20    END DO 
         DO 30 J = M + N, N, - 1 
            FRE (J) = FRE (J - N) 
            FIM (J) = FIM (J - N) 
   30    END DO 
         DO 40 J = 0, N - 1 
            FRE (J) = 0.0D0 
            FIM (J) = 0.0D0 
   40    END DO 
      ENDIF 
      DO 50 J = N + 1, L - 1 
         HRE (J) = 0.0D0 
         HIM (J) = 0.0D0 
   50 END DO 
!                                                                       
!*****************************************************************      
!     Compute necessary discrete Fourier coefficients via FFT    *      
!     for powers of two                                          *      
!*****************************************************************      
!                                                                       
      CALL FFT (ITAU, FRE, FIM, 0) 
      CALL FFT (ITAU, HRE, HIM, 0) 
!                                                                       
      IF (IND.EQ.0) THEN 
         DO 100 K = 0, L - 1 
            H = CMLR (FRE (K), FIM (K), HRE (K), HIM (K) ) 
            FIM (K) = CMLI (FRE (K), FIM (K), HRE (K), HIM (K) ) 
            FRE (K) = H 
  100    END DO 
      ELSE 
         DO 110 K = 0, L - 1 
            H = CMLR (FRE (K), FIM (K), HRE (K), - HIM (K) ) 
            FIM (K) = CMLI (FRE (K), FIM (K), HRE (K), - HIM (K) ) 
            FRE (K) = H 
  110    END DO 
      ENDIF 
!                                                                       
!*****************************************************************      
!     Compute corresponding functional values                    *      
!*****************************************************************      
!                                                                       
      CALL FFT (ITAU, FRE, FIM, 1) 
!                                                                       
      X = DBLE (L) * DELTAX 
      DO 200 J = 0, M + N 
         FRE (J) = X * FRE (J) 
         FIM (J) = X * FIM (J) 
  200 END DO 
      DO 210 J = M + N + 1, L - 1 
         FRE (J) = 0.0D0 
         FIM (J) = 0.0D0 
  210 END DO 
      RETURN 
      END SUBROUTINE FFAKON                         
