![KA{P 3}{Roots of Polynomials}                                         
![       {Roots of Polynomials}*)                                       
      SUBROUTINE MULLRP (NPOL, POLYNM, ITERMX, NZ, ZERO, P, PHELPD) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This SUBROUTINE computes all zeros of a polynomial by using   *      
!  the method of Muller.                                         *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  NPOL   : degree of the polynomial                             *      
!  POLYNM : (NPOL+1)-vector POLYNM(0:NPOL) of the polynomial     *      
!           coefficients in ascending order                      *      
!  ITERMX : maximum number of iterations per zero                *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  NZ     : number of zeros found                                *      
!  ZERO   : (2,NZ)-array containing the zeros with real their    *      
!           parts stored first and their imaginary parts stored  *      
!           in the second component.                             *      
!                                                                *      
!                                                                *      
!  AUXILIARY VECTORS:                                            *      
!  ==================                                            *      
!  P      : vector P(0:NPOL)  of type DOUBLE PRECISION           *      
!  PHELPD : vector PHELPD(0:NPOL)  of type DOUBLE PRECISION      *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: MULLER, HORNC, HORNCE, POLDIV,          *      
!                        YEPS, COMPAR                            *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Eberhard Heyne                                   *      
!  date       : 06.14.1992                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER NPOL, NFND, ITERMX, IERR, NZ, NC, N, I, ITMX, NB, NFIRST, &
      NR                                                                
      DOUBLEPRECISION POLYNM (0:NPOL) 
      DOUBLEPRECISION P (0:NPOL) 
      DOUBLEPRECISION ZERO (2, 1:NPOL), X0R, X1R, X2R, X0I, X1I, X2I 
      DOUBLEPRECISION F0R, F1R, F2R, F0I, F1I, F2I, F1BQ, F2BQ, F2BQB 
      DOUBLEPRECISION XR, XI, XEPS, YEPS 
      DOUBLEPRECISION BHELPD (0:2), PHELPD (0:NPOL), R (0:1) 
      COMMON / MULLCI / ITMX 
      COMMON / MULLCD / X0R, X0I, X1R, X1I, X2R, X2I, F0R, F0I, F1R,    &
      F1I, F2R, F2I, F1BQ, F2BQ, F2BQB, XEPS                            
!                                                                       
!     The variable IERR (error indicator) may be added to the parameter 
!     as an additional output parameter. Since the number NPOL of zeros 
!     of the given polynomial is known, the inquiry (NZ .NE. NPOL) might
!     suffice if the presence of an error needs to be detected but      
!     the type of this error does not need to be analyzed.              
!                                                                       
      DATA NFIRST / 0 / 
      IF (NFIRST.NE.1) THEN 
!                                                                       
!        Determine machine constant                                     
!                                                                       
         XEPS = YEPS () 
         NFIRST = 1 
      ENDIF 
!                                                                       
!     Initialize: number NZ of zeros found, error indicator IERR        
!                                                                       
      NZ = 0 
      IERR = 0 
!                                                                       
!     Test for meaningful polynomial degree                             
!                                                                       
      IF (NPOL.LE.0) THEN 
         IERR = 1 
         RETURN 
      ENDIF 
!                                                                       
!     Relabel polynomial coefficients and determine the polynomial degre
!                                                                       
      N = - 1 
      DO 10 I = 0, NPOL 
         P (I) = POLYNM (I) 
         IF (DABS (P (I) ) .NE.0.0D0) N = I 
   10 END DO 
      IF (N.LE.0) THEN 
         IF (N.EQ. - 1) THEN 
!                                                                       
!           The polynomial is identically equal to zero                 
!                                                                       
            IERR = 2 
            RETURN 
         ELSEIF (N.EQ.0) THEN 
!                                                                       
!           The polynomial is identical to a constant different from zer
!                                                                       
            IERR = 3 
            RETURN 
         ENDIF 
      ENDIF 
   12 CONTINUE 
      IF (N.EQ.0) THEN 
         IERR = 4 
         RETURN 
      ENDIF 
      IF (N.EQ.1) THEN 
!                                                                       
!        Solve the linear polynomial if the degree is 1                 
!                                                                       
         NZ = NZ + 1 
         ZERO (1, NZ) = - P (0) / P (1) 
         ZERO (2, NZ) = 0.0D0 
         RETURN 
      ENDIF 
!                                                                       
!     Preset the number ITMX of iteration steps allowed per zero        
!                                                                       
      ITMX = ITERMX 
!                                                                       
!     Automatic start                                                   
!                                                                       
      X0R = - 1.0D0 
      F0R = P (0) - P (1) + P (2) 
      X1R = 1.0D0 
      F1R = P (0) + P (1) + P (2) 
      X2R = 0.0D0 
      F2R = P (0) 
      X0I = 0.0D0 
      F0I = 0.0D0 
      X1I = 0.0D0 
      F1I = 0.0D0 
      X2I = 0.0D0 
      F2I = 0.0D0 
!                                                                       
!     Muller-iteration for one zero                                     
!                                                                       
      CALL MULLER (N, P, NFND, XR, XI, PHELPD) 
      IF (NFND.EQ.0) RETURN 
      IF (NFND.EQ.1) THEN 
!                                                                       
!        A real zero                                                    
!                                                                       
         NZ = NZ + 1 
         ZERO (1, NZ) = XR 
         ZERO (2, NZ) = 0.0D0 
         BHELPD (1) = 1.0D0 
         BHELPD (0) = - XR 
         NB = 1 
      ELSEIF (NFND.EQ.2) THEN 
!                                                                       
!        If a real polynomial has a complex zero,                       
!        then its complex-conjugate is a zero as well                   
!                                                                       
         NZ = NZ + 1 
         ZERO (1, NZ) = XR 
         ZERO (2, NZ) = XI 
         NZ = NZ + 1 
         ZERO (1, NZ) = XR 
         ZERO (2, NZ) = - XI 
         BHELPD (2) = 1.0D0 
         BHELPD (1) = - (XR + XR) 
         BHELPD (0) = + (XR**2 + XI**2) 
         NB = 2 
      ENDIF 
!                                                                       
!     Deflate by the found zero(s)                                      
!                                                                       
      CALL POLDIV (P, N, BHELPD, NB, PHELPD, NC, R, NR) 
      DO 14 I = 0, NC 
         P (I) = PHELPD (I) 
   14 END DO 
      N = NC 
      IF (ITMX.NE.0) GOTO 12 
      RETURN 
      END SUBROUTINE MULLRP                         
!                                                                       
!                                                                       
      SUBROUTINE MULLER (N, P, NFND, XR, XI, PHELPD) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Auxiliary routine for MULLRP                                  *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: HORNC, HORNCE                           *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Eberhard Heyne                                   *      
!  date       : 06.14.1992                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER N, NFND, ITMX 
      DOUBLEPRECISION P (0:N), F1BQ 
      DOUBLEPRECISION F2BQ, F2BQB, XEPS, XEFR, XEFI, XEB1, XEB0, FF 
      DOUBLEPRECISION H1R, H1I, H1BQ, H2R, H2I, HHR, HHI, HXR, HXI, HYR,&
      HYI                                                               
      DOUBLEPRECISION WSIN, WCOS, WBQ, WSINH, WCOSH 
      DOUBLEPRECISION AR, AI, BR, BI, CR, CI, QR, QI, WNR, WNI 
      DOUBLEPRECISION XN1R, XN1I, XN1BQ, XN2R, XN2I, XN2BQ 
      DOUBLEPRECISION XR, XI, FR, FI, XBR, XBI, PHELPD (0:N) 
      DOUBLEPRECISION X0R, X1R, X2R, X0I, X1I, X2I 
      DOUBLEPRECISION F0R, F1R, F2R, F0I, F1I, F2I 
      DOUBLEPRECISION AHR, AHI, BHR, BHI, CMLR, CMLI 
      COMMON / MULLCI / ITMX 
      COMMON / MULLCD / X0R, X0I, X1R, X1I, X2R, X2I, F0R, F0I, F1R,    &
      F1I, F2R, F2I, F1BQ, F2BQ, F2BQB, XEPS                            
!                                                                       
!     Complex multiplication: the real part                             
!                                                                       
      CMLR (AHR, AHI, BHR, BHI) = AHR * BHR - AHI * BHI 
!                                                                       
!     Complex multiplication: the imaginary part                        
!                                                                       
      CMLI (AHR, AHI, BHR, BHI) = AHR * BHI + AHI * BHR 
!                                                                       
!     The number FF helps decide whether one should iterate further     
!     when the computed functional values increase or whether the value 
!     with the smallest absolute functional value so far attained is    
!     acceptable as a zero.                                             
!     FF=10 means that one decimal place of the computed functional valu
!     must be correct. FF=1 means that we only test whether the absolute
!     value of the functional value is less than the rounding error.    
!                                                                       
      DATA FF / 10.0D0 / 
!                                                                       
!     Automatic start                                                   
!                                                                       
!     H1=X1-X0                                                          
!     H2=X2-X1                                                          
!     Q=H2/H1                                                           
!                                                                       
      H1R = X1R - X0R 
      H1I = 0.0D0 
      H2R = X2R - X1R 
      H2I = 0.0D0 
!                                                                       
!     Q=H2R/H1R                                                         
!                                                                       
      H1BQ = H1R**2 + H1I**2 
      QR = CMLR (H2R, H2I, H1R / H1BQ, - H1I / H1BQ) 
      QI = CMLI (H2R, H2I, H1R / H1BQ, - H1I / H1BQ) 
      F1BQ = 1.0D38 
      F2BQ = 1.0D38 
      XBR = X2R 
      XBI = 0.0D0 
      F2BQB = 1.0D38 
   10 CONTINUE 
!                                                                       
!     A=Q*F2-Q*(1+Q)*F1+Q**2*F0                                         
!     B=(2.*Q+1.)*F2-(1+Q)**2*F1+Q**2*F0                                
!     C=(1.+Q)*F2                                                       
!                                                                       
      HHR = CMLR (QR, QI, QR, QI) 
      HHI = CMLI (QR, QI, QR, QI) 
!                                                                       
!     Compute        A=Q*(F2-F1) + Q**2*(F0-F1)                         
!     instead of     A=Q*F2-Q*(1+Q)*F1+Q**2*F0                          
!                                                                       
      HXR = F2R - F1R 
      HXI = F2I - F1I 
      HYR = F0R - F1R 
      HYI = F0I - F1I 
      AR = CMLR (QR, QI, HXR, HXI) + CMLR (HHR, HHI, HYR, HYI) 
      AI = CMLI (QR, QI, HXR, HXI) + CMLI (HHR, HHI, HYR, HYI) 
!                                                                       
!     Compute         B=F2-F1+ 2.*Q*(F2-F1)+Q**2*(F0-F1)                
!     instead of      B=(2.*Q+1.)*F2-(1+Q)**2*F1+Q**2*F0                
!                                                                       
      BR = HXR + 2.0D0 * CMLR (QR, QI, HXR, HXI) + CMLR (HHR, HHI, HYR, &
      HYI)                                                              
      BI = HXI + 2.0D0 * CMLI (QR, QI, HXR, HXI) + CMLI (HHR, HHI, HYR, &
      HYI)                                                              
!                                                                       
!     C=(1.+Q)*F2                                                       
!                                                                       
      CR = F2R + CMLR (QR, QI, F2R, F2I) 
      CI = F2I + CMLI (QR, QI, F2R, F2I) 
!                                                                       
!     The square-root expression appearing in the denominator           
!     of the expression for Q                                           
!                                                                       
!     WN=B**2-4.*A*C  und WN=SQRT(WN)                                   
!                                                                       
      WNR = CMLR (BR, BI, BR, BI) - 4.0D0 * CMLR (AR, AI, CR, CI) 
      WNI = CMLI (BR, BI, BR, BI) - 4.0D0 * CMLI (AR, AI, CR, CI) 
      IF (WNR * WNR + WNI * WNI.GE.1.0D-40) THEN 
         WBQ = DSQRT (WNR * WNR + WNI * WNI) 
         WSIN = WNI / WBQ 
         WCOS = WNR / WBQ 
!                                                                       
!        Half-angle formula                                             
!                                                                       
         WSINH = DSIGN (1.0D0, WSIN) * DSQRT ( (1.0D0 - WCOS) * 0.5D0) 
         WCOSH = DSQRT ( (1.0D0 + WCOS) * 0.5D0) 
         WBQ = DSQRT (WBQ) 
         WNR = WBQ * WCOSH 
         WNI = WBQ * WSINH 
      ELSE 
         WNR = 0.0D0 
         WNI = 0.0D0 
      ENDIF 
!                                                                       
!     The two possible denominators for this root                       
!                                                                       
      XN1R = BR + WNR 
      XN1I = BI + WNI 
      XN1BQ = XN1R**2 + XN1I**2 
      XN2R = BR - WNR 
      XN2I = BI - WNI 
      XN2BQ = XN2R**2 + XN2I**2 
!                                                                       
!     The denominator with the larger absolute value determines the X   
!     closest to X2                                                     
!                                                                       
      IF (XN1BQ.GT.XN2BQ) THEN 
!                                                                       
!        Q=-2.*C/XN1                                                    
!                                                                       
         QR = - (2.0D0 / XN1BQ) * CMLR (CR, CI, XN1R, - XN1I) 
         QI = - (2.0D0 / XN1BQ) * CMLI (CR, CI, XN1R, - XN1I) 
      ELSEIF (XN2BQ.GT.0.0D0) THEN 
!                                                                       
!        Q=-2.*C/XN2                                                    
!                                                                       
         QR = - (2.0D0 / XN2BQ) * CMLR (CR, CI, XN2R, - XN2I) 
         QI = - (2.0D0 / XN2BQ) * CMLI (CR, CI, XN2R, - XN2I) 
      ELSE 
!                                                                       
!        The denominator is zero;                                       
!        Follow the suggestion by Muller: set Q=1 and continue calculati
!                                                                       
!        Q=(1.,0.)                                                      
!                                                                       
         QR = 1.0D0 
         QI = 0.0D0 
      ENDIF 
!                                                                       
!     Prepare the next iteration,                                       
!     in which some instructions become redundant; we will              
!     label them as comment lines: C   ...                              
!                                                                       
!     X0=X1                                                             
!     X1=X2                                                             
!     H1=H2                                                             
!                                                                       
      F0R = F1R 
      F0I = F1I 
      F1R = F2R 
      F1I = F2I 
!     X0R=X1R                                                           
!     X0I=X1I                                                           
!     X1R=X2R                                                           
!     X1I=X2I                                                           
!     H1R=H2R                                                           
!     H1I=H2I                                                           
!                                                                       
!     The new value of H2 is calculated before the new X iterate        
!                                                                       
      F1BQ = F2BQ 
!                                                                       
!     H2=H2*Q                                                           
!                                                                       
      HHR = CMLR (H2R, H2I, QR, QI) 
      H2I = CMLI (H2R, H2I, QR, QI) 
      H2R = HHR 
      X2R = X2R + H2R 
      X2I = X2I + H2I 
   12 CONTINUE 
!                                                                       
!     Determine the functional value                                    
!                                                                       
      XR = X2R 
      XI = X2I 
      CALL HORNC (P, N, XR, XI, FR, FI, PHELPD) 
      F2R = FR 
      F2I = FI 
      F2BQ = FR**2 + FI**2 
!                                                                       
!     Decrease the iteration counter                                    
!                                                                       
      ITMX = ITMX - 1 
      IF (ITMX.EQ.0) THEN 
!                                                                       
!        Maximal number of iterations has been reached                  
!                                                                       
!                                                                       
!        Error estimate for the x-value with the minimal                
!        absolute value of the function                                 
!                                                                       
         CALL HORNCE (XEPS, XEFR, XEFI, XEB1, XEB0, P, N, XBR, XBI, FR, &
         FI, PHELPD)                                                    
!                                                                       
!        Determine whether the x-value with the smallest absolute       
!        functional value so far found can be regarded as a zero        
!                                                                       
         IF (DABS (FR) + DABS (FI) .LE.FF * (DABS (XEFR) + DABS (XEFI) )&
         ) THEN                                                         
            XR = XBR 
            XI = XBI 
            GOTO 14 
         ELSE 
            NFND = 0 
            RETURN 
         ENDIF 
      ENDIF 
!                                                                       
!     Muller-modification for improved convergence                      
!                                                                       
      IF (F2BQ.GT.100.0D0 * F1BQ) THEN 
         QR = QR / 2.0D0 
         QI = QI / 2.0D0 
!                                                                       
!        H2 and X2 already contain the old Q                            
!                                                                       
         H2R = H2R / 2.0D0 
         H2I = H2I / 2.0D0 
         X2R = X2R - H2R 
         X2I = X2I - H2I 
         GOTO 12 
      ENDIF 
!                                                                       
!     As long as the absolute value of the functional values decreases, 
!     we can hope for an improvement                                    
!                                                                       
      IF (F2BQ.LT.F1BQ) THEN 
         IF (F2BQ.LT.F2BQB) THEN 
!                                                                       
!          The absolute value of the new function value is less than    
!          the minimal value so far                                     
!                                                                       
            F2BQB = F2BQ 
            XBR = XR 
            XBI = XI 
         ENDIF 
!                                                                       
!       We do not continue iterating for an exakt zero                  
!                                                                       
         IF (F2BQ.NE.0.0D0) GOTO 10 
      ELSE 
!                                                                       
!        Error estimate the same as for the x-value that, so far, repres
!        the minimal absolute function value                            
!                                                                       
         CALL HORNCE (XEPS, XEFR, XEFI, XEB1, XEB0, P, N, XBR, XBI, FR, &
         FI, PHELPD)                                                    
!                                                                       
!        Check whether iteration is to be continued                     
!                                                                       
         XR = XBR 
         XI = XBI 
         IF (DABS (FR) + DABS (FI) .GE.FF * (DABS (XEFR) + DABS (XEFI) )&
         ) GOTO 10                                                      
      ENDIF 
   14 CONTINUE 
!                                                                       
!     Using an error estimate we determine whether the root             
!     may be real                                                       
!                                                                       
      IF (DABS (PHELPD (1) ) + DABS (PHELPD (0) ) .LE.FF * (XEB1 + XEB0)&
      ) THEN                                                            
!                                                                       
!        Complex-conjugate pair of zeros                                
!                                                                       
         NFND = 2 
      ELSE 
!                                                                       
!        A real solution                                                
!                                                                       
         NFND = 1 
         XI = 0.0D0 
      ENDIF 
      RETURN 
      END SUBROUTINE MULLER                         
!                                                                       
!                                                                       
      SUBROUTINE HORNC (A, NA, XR, XI, FR, FI, B) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  auxiliary routine for MULLRP                                  *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Eberhard Heyne                                   *      
!  date       : 05.09.1988                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
!     Simple Horner-scheme for a COMPLEX argument                       
!                                                                       
      INTEGER NA, I 
      DOUBLEPRECISION A (0:NA) 
      DOUBLEPRECISION XR, XI, FR, FI, P, Q, B (0:NA) 
      P = XR + XR 
      Q = - (XR**2 + XI**2) 
      B (NA) = A (NA) 
      B (NA - 1) = A (NA - 1) + P * B (NA) 
      DO 10 I = NA - 2, 1, - 1 
         B (I) = A (I) + Q * B (I + 2) + P * B (I + 1) 
   10 END DO 
      B (0) = A (0) + Q * B (2) 
      FR = B (1) * XR + B (0) 
      FI = B (1) * XI 
      RETURN 
      END SUBROUTINE HORNC                          
!                                                                       
!                                                                       
      SUBROUTINE HORNCE (XEPS, XEFR, XEFI, XEB1, XEB0, A, NA, XR, XI,   &
      FR, FI, B)                                                        
!                                                                       
!*****************************************************************      
!                                                                *      
!  auxiliary routine for MULLRP                                  *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Eberhard Heyne                                   *      
!  date       : 06.14.1992                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
!     This SUBROUTINE is very similar to HORNC except that              
!     error estimates are carried out for each instruction.             
!                                                                       
!     XEFR   error estimate of the real part of F      ( FR )           
!     XEFI   error estimate of the imaginary part F of ( FI )           
!     XEB1   error estimate of B(1)                                     
!     XEB0   error estimate of B(0)                                     
!                                                                       
      INTEGER NA, I 
      DOUBLEPRECISION A (0:NA) 
      DOUBLEPRECISION XEPS, XEP, XEQ, XEB0, XEB1, XEB2, XEFR, XEFI 
      DOUBLEPRECISION XR, XI, FR, FI, P, Q, B (0:NA) 
      P = XR + XR 
      XEP = DABS (P) * XEPS 
      Q = - (XR**2 + XI**2) 
      XEQ = DABS (Q) * XEPS 
      B (NA) = A (NA) 
      XEB2 = 0.0D0 
      B (NA - 1) = A (NA - 1) + P * B (NA) 
      XEB1 = XEP * DABS (B (NA) ) + XEPS * (DABS (B (NA - 1) ) ) 
      DO 10 I = NA - 2, 1, - 1 
         B (I) = A (I) + Q * B (I + 2) + P * B (I + 1) 
         XEB0 = XEQ * DABS (B (I + 2) ) + XEP * DABS (B (I + 1) )       &
         + XEB2 * DABS (Q) + XEB1 * DABS (P) + XEPS * (DABS (B (I) )    &
         + DABS (A (I) ) )                                              
         XEB2 = XEB1 
         XEB1 = XEB0 
   10 END DO 
      B (0) = A (0) + Q * B (2) 
      XEB0 = XEQ * DABS (B (2) ) + XEB2 * DABS (Q) + XEPS * (DABS (B (0)&
      ) + DABS (A (0) ) )                                               
      FR = B (1) * XR + B (0) 
      XEFR = XEB1 * DABS (XR) + XEB0 
      FI = B (1) * XI 
      XEFI = XEB1 * DABS (XI) 
      RETURN 
      END SUBROUTINE HORNCE                         
!                                                                       
!                                                                       
      SUBROUTINE POLDIV (A, NA, B, NB, C, NC, R, NR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  auxiliary routine for MULLRP                                  *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Eberhard Heyne                                   *      
!  date       : 05.09.1988                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER NA, NB, NC, NR, IC, K, M, IR 
      DOUBLEPRECISION A (0:NA) 
      DOUBLEPRECISION B (0:NB), C (0: * ), SUM, R (0: * ) 
!                                                                       
!     Long division of polynomials:                                     
!     determine polynomials C, R so that the following                  
!     holds:  A = B * C + R, i.e., C = A/B with remainder R             
!                                                                       
      NC = NA - NB 
      DO 10 IC = NC, 0, - 1 
         K = IC + NB 
         SUM = A (K) 
         DO 12 M = IC + 1, MIN (NC, K) 
            SUM = SUM - B (K - M) * C (M) 
   12    END DO 
         C (IC) = SUM / B (NB) 
   10 END DO 
      NR = NB - 1 
      DO 14 IR = 0, NR 
         SUM = A (IR) 
         DO 16 M = 0, IR 
            SUM = SUM - B (IR - M) * C (M) 
   16    END DO 
         R (IR) = SUM 
   14 END DO 
      RETURN 
      END SUBROUTINE POLDIV                         
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION YEPS () 
!                                                                       
!*****************************************************************      
!                                                                *      
!  auxiliary routine for MULLRP                                  *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: COMPAR                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Eberhard Heyne                                   *      
!  date       : 05.09.1988                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION R, S 
      INTEGER M 
!                                                                       
!     Function that determines the machine constant                     
!                                                                       
      S = 1.0D0 
   10 R = S 
      S = S / 2.0D0 
      CALL COMPAR (M, 1.0D0 + S, 1.0D0 + R) 
      IF (M.NE.0) GOTO 10 
      YEPS = R 
      RETURN 
      END FUNCTION YEPS                             
!                                                                       
!                                                                       
      SUBROUTINE COMPAR (M, A, B) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  auxiliary routine for MULLRP                                  *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION A, B 
      INTEGER M 
!                                                                       
!     The reason behind this routine is to avoid internal               
!     compiler optimizations and to force storing the values            
!     for A and B.                                                      
!                                                                       
      M = 1 
      IF (A.EQ.B) M = 0 
      RETURN 
      END SUBROUTINE COMPAR                         
