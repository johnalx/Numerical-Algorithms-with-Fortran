      SUBROUTINE BAUPOL (COEFRE, COEFIM, N, COMPL, ROOTRE, ROOTIM,      &
      NUMIT)                                                            
!                                                                       
!*****************************************************************      
!                                                                *      
!  Without knowing any approximations of the roots, this         *      
!  SUBROUTINE finds N approximate values Z(I), I=1, ..., N for   *      
!  the N zeros of a polynomial PN of degree N with real or       *      
!  complex coefficients. The polynomial is described as follows: *      
!                                                                *      
!     PN(Z)=COEF(0)+COEF(1)*Z+COEF(2)*Z**2+...+COEF(N)*Z**N,     *      
!                                                                *      
!  with COEF(I) = (COEFRE(I),COEFIM(I)) for I=0, ..., N (complex *      
!  coefficients).                                                *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  COEFRE   : (N+1)-vector COEFRE(0:N) containing the real       *      
!             part of each coefficient of the polynomial PN in   *      
!             DOUBLE PRECISION.                                  *      
!  COEFIM   : (N+1)-vector COEFIM(0:N) containing the imaginary  *      
!             part of each coefficient of the polynomial PN in   *      
!             DOUBLE PRECISION.                                  *      
!  N        : degree of the polynomial PN.                       *      
!  COMPL    : boolean variable :                                 *      
!             COMPL=.TRUE.  , if the coefficients are COMPLEX.   *      
!             COMPL=.FALSE. , if the coefficients are REAL.      *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  ROOTRE   : N-vector ROOTRE(1:N) containing the approximate    *      
!             real parts of the computed zeros of PN in          *      
!             DOUBLE PRECISION.                                  *      
!  ROOTIM   : N-vector ROOTIM(1:N) containing the approximate    *      
!             imaginary parts of the computed zeros of PN in     *      
!             DOUBLE PRECISION.                                  *      
!  NUMIT    : maximum number of iteration steps.                 *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: BAUZRO, ABSCOM, CDIV, SCALFC, COMHOR    *      
!                        MCONST.                                 *      
!                                                                *      
!                                                                *      
!  sources: Bauhuber, see [BAUH70].                              *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Guido Dubois                                     *      
!  date       : 11.01.1985                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      COMMON / GLOBAL / A, B, C 
      DOUBLEPRECISION COEFRE (0:100), COEFIM (0:100), E (101), ROOTRE ( &
      100), ROOTIM (100), A (202), B (202), C (200), X0, Y0, XNEW, YNEW,&
      GAMMA, FMACHP, INFINY, SMALNO, BASE, SCALFC, BND, MAX, MIN, X,    &
      ABSCOM                                                            
      LOGICAL COMPL 
!                                                                       
!  Initializing the iteration step counter NUMIT and the error          
!  bound GAMMA.                                                         
!                                                                       
      NUMIT = 0 
      GAMMA = 5.0D-18 
      IF (N.EQ.1) THEN 
!                                                                       
!  If the degree of PN is N=1, then the zero of the polynomial PN is    
!  Z=-COEF(0)/COEF(1), where COEF(I)=(COEFRE(I),COEFIM(I)) for I=0,1.   
!                                                                       
         CALL CDIV ( - COEFRE (0), - COEFIM (0), COEFRE (1), COEFIM (1),&
         ROOTRE (1), ROOTIM (1) )                                       
         RETURN 
      ELSE 
         N1 = N + 1 
!                                                                       
!  Scaling via SCALFC.                                                  
!                                                                       
         DO 10 I = 1, N1 
            E (I) = ABSCOM (COEFRE (N1 - I), COEFIM (N1 - I) ) 
   10    END DO 
         CALL MCONST (FMACHP, INFINY, SMALNO, BASE) 
         BND = SCALFC (N1, E, FMACHP, INFINY, SMALNO, BASE) 
         IF (BND.EQ.1.0D0) THEN 
!                                                                       
!  Normalizing, in case scaling by SCALFC did not normalize the coeffici
!                                                                       
            MAX = 0.0D0 
            MIN = 1.0D+300 
            DO 20 I = N, 0, - 1 
               X = ABSCOM (COEFRE (I), COEFIM (I) ) 
               IF (X.GT.MAX) MAX = X 
               IF (X.LT.MIN.AND.X.NE.0.0D0) MIN = X 
   20       END DO 
            BND = DSQRT (MAX * MIN) 
            BND = 1.0D0 / BND 
         ENDIF 
         DO 30 K = 1, N1 
            B (2 * K - 1) = COEFRE (N1 - K) * BND 
            B (2 * K) = 0.0D0 
            IF (COMPL) B (2 * K) = COEFIM (N1 - K) * BND 
   30    END DO 
         X0 = 0.0D0 
         Y0 = 0.0D0 
         DO 40 I = 1, N 
            L = 2 * (N + 2 - I) 
            DO 50 K = 1, L 
               A (K) = B (K) 
   50       END DO 
!                                                                       
!  Calculating the I-th zero of PN.                                     
!                                                                       
            CALL BAUZRO (X0, Y0, N + 1 - I, GAMMA, XNEW, YNEW, NUMIT) 
            ROOTRE (I) = XNEW 
            ROOTIM (I) = YNEW 
            X0 = XNEW 
            Y0 = - YNEW 
   40    END DO 
      ENDIF 
      RETURN 
      END SUBROUTINE BAUPOL                         
!                                                                       
!                                                                       
      SUBROUTINE BAUZRO (X0, Y0, N, GAMMA, XNEW, YNEW, NUMIT) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This SUBROUTINE calculates a zero of a polynomial PN with     *      
!  complex coefficients. We solve the equation PN(Z)/PN'(Z)=0    *      
!  via Newton's method with spiralization and extrapolation to   *      
!  improve convergence.                                          *      
!  The initial approximation (X0+I*Y0) can be chosen arbitrarily.*      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  X0       : real part of the initial approximation.            *      
!  Y0       : imaginary part of the initial approximation.       *      
!  N        : degree of the polynomial PN.                       *      
!  GAMMA    : error bound.                                       *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  XNEW     : real part of the computed zero of PN.              *      
!  YNEW     : imaginary part of the zero of PN.                  *      
!  NUMIT    : maximum number of iteration steps.                 *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: ABSCOM, CDIV, COMHOR.                   *      
!                                                                *      
!                                                                *      
!  sources: Bauhuber, see [BAUH70].                              *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Guido Dubois                                     *      
!  date       : 11.01.1985                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      COMMON / GLOBAL / A, B, C 
      DOUBLEPRECISION A (202), B (202), C (200), ABSCOM, X0, Y0, XNEW,  &
      YNEW, GAMMA, BETA, RHO, BDZE, QR, QI, UNEW, VNEW, UDNEW, VDNEW,   &
      UDDNEW, VDDNEW, PBNEW, PBOLD, PSBNEW, BD, BDD, DZMAX, DZMIN, DX,  &
      DY, U, V, XOLD, YOLD, H, H1, H2, H3, H4, H5, P1RE, P1IM, Q1RE,    &
      Q1IM, P12RE, P12IM, RA1RE, RA1IM, RARE, RAIM, RABE, RTRE, RTIM    
      LOGICAL ENDIT 
!                                                                       
      IF (N.EQ.2) THEN 
!                                                                       
!  Solving the remaining 2nd degree polynomial exactly.                 
!                                                                       
         CALL CDIV (A (3), A (4), A (1), A (2), P1RE, P1IM) 
         CALL CDIV (A (5), A (6), A (1), A (2), Q1RE, Q1IM) 
         P12RE = - P1RE / 2.0D0 
         P12IM = - P1IM / 2.0D0 
         RA1RE = P12RE * P12RE-P12IM * P12IM 
         RA1IM = 2.0D0 * P12RE * P12IM 
         RARE = RA1RE-Q1RE 
         RAIM = RA1IM - Q1IM 
         IF (RAIM.EQ.0.0D0) THEN 
            IF (RARE.LT.0.0D0) THEN 
!                                                                       
!  Purely imaginary root.                                               
!                                                                       
               RTIM = DSQRT ( - RARE) 
               RTRE = 0.0D0 
               XNEW = P12RE 
               YNEW = P12IM + RTIM 
               RETURN 
            ELSE 
!                                                                       
!  Real root.                                                           
!                                                                       
               RTRE = DSQRT (RARE) 
               RTIM = 0.0D0 
               XNEW = P12RE+RTRE 
               YNEW = P12IM 
               RETURN 
            ENDIF 
         ELSE 
!                                                                       
!  Complex root.                                                        
!                                                                       
            RABE = ABSCOM (RARE, RAIM) 
            IF (RARE.GT.0.0D0) THEN 
               RTRE = DSQRT (0.5D0 * (RABE+RARE) ) 
               IF (RAIM.LT.0.0D0) RTRE = - RTRE 
               RTIM = 0.5D0 * RAIM / RTRE 
            ELSE 
               RTIM = DSQRT (0.5D0 * (RABE-RARE) ) 
               RTRE = 0.5D0 * RAIM / RTIM 
            ENDIF 
            XNEW = P12RE+RTRE 
            YNEW = P12IM + RTIM 
            RETURN 
         ENDIF 
      ELSEIF (N.EQ.1) THEN 
!                                                                       
!  Polynomial of 1st degree.                                            
!                                                                       
         XNEW = P12RE-RTRE 
         YNEW = P12IM - RTIM 
         RETURN 
      ELSE 
         I = 0 
         ENDIT = .FALSE. 
         RHO = DSQRT (GAMMA) 
         BETA = 10.0D0 * GAMMA 
         QR = 0.1D0 
         QI = 0.9D0 
         XNEW = X0 
         YNEW = Y0 
         CALL COMHOR (XNEW, YNEW, N, GAMMA, UNEW, VNEW, UDNEW, VDNEW,   &
         UDDNEW, VDDNEW, BD, BDD)                                       
         NUMIT = NUMIT + 1 
         PBNEW = ABSCOM (UNEW, VNEW) 
         IF (PBNEW.LE.BD) THEN 
            RETURN 
         ELSE 
            PBOLD = 2.0D0 * PBNEW 
            DZMIN = BETA * (RHO + ABSCOM (XNEW, YNEW) ) 
   10       PSBNEW = ABSCOM (UDNEW, VDNEW) 
!                                                                       
!  Spiralization.                                                       
!                                                                       
            IF (PBNEW.LT.PBOLD) THEN 
               DZMAX = 1.0D0 + ABSCOM (XNEW, YNEW) 
               NUMIT = NUMIT + 1 
               H1 = UDNEW * UDNEW - VDNEW * VDNEW - UNEW * UDDNEW +     &
               VNEW * VDDNEW                                            
               H2 = 2.0D0 * UDNEW * VDNEW - UNEW * VDDNEW - VNEW *      &
               UDDNEW                                                   
               IF (PSBNEW.GT.10.0D0 * BDD.AND.ABSCOM (H1, H2)           &
               .GT.100.0D0 * BDD * BDD) THEN                            
!                                                                       
!  Applying Newton's method.                                            
!                                                                       
                  I = I + 1 
                  IF (I.GT.2) I = 2 
                  U = UNEW * UDNEW - VNEW * VDNEW 
                  V = UNEW * VDNEW + VNEW * UDNEW 
                  CALL CDIV ( - U, - V, H1, H2, DX, DY) 
                  IF (ABSCOM (DX, DY) .GT.DZMAX) THEN 
                     H = DZMAX / ABSCOM (DX, DY) 
                     DX = DX * H 
                     DY = DY * H 
                     I = 0 
                  ENDIF 
                  IF (I.EQ.2.AND.ABSCOM (DX, DY) .LT.DZMIN /            &
                  RHO.AND.ABSCOM (DX, DY) .GT.0.0D0) THEN               
!                                                                       
!  Extrapolation.                                                       
!                                                                       
                     I = 0 
                     CALL CDIV (XNEW - XOLD, YNEW - YOLD, DX, DY, H3,   &
                     H4)                                                
                     H3 = 1.0D0 + H3 
                     H1 = H3 * H3 - H4 * H4 
                     H2 = 2.0D0 * H3 * H4 
                     CALL CDIV (DX, DY, H1, H2, H3, H4) 
                     IF (ABSCOM (H3, H4) .LT.50.0D0 * DZMIN) THEN 
                        DX = DX + H3 
                        DY = DY + H4 
                     ENDIF 
                  ENDIF 
                  XOLD = XNEW 
                  YOLD = YNEW 
                  PBOLD = PBNEW 
               ELSE 
!                                                                       
!  In a neighborhood of a saddle point.                                 
!                                                                       
                  I = 0 
                  H = DZMAX / PBNEW 
                  DX = H * UNEW 
                  DY = H * VNEW 
                  XOLD = XNEW 
                  YOLD = YNEW 
                  PBOLD = PBNEW 
   20             CALL COMHOR (XNEW + DX, YNEW + DY, N, GAMMA, U, V, H, &
                  H1, H2, H3, H4, H5)                                   
                  IF (DABS (ABSCOM (U, V) / PBNEW - 1.0D0) .LE.RHO)     &
                  THEN                                                  
                     DX = 2.0D0 * DX 
                     DY = 2.0D0 * DY 
                     GOTO 20 
                  ENDIF 
               ENDIF 
            ELSE 
               I = 0 
               NUMIT = NUMIT + 1 
               H = QR * DX - QI * DY 
               DY = QR * DY + QI * DX 
               DX = H 
            ENDIF 
            IF (ENDIT) THEN 
               IF (ABSCOM (DX, DY) .LT.0.1D0 * BDZE) THEN 
                  XNEW = XNEW + DX 
                  YNEW = YNEW + DY 
               ENDIF 
               CALL COMHOR (XNEW, YNEW, N, GAMMA, UNEW, VNEW, UDNEW,    &
               VDNEW, UDDNEW, VDDNEW, BD, BDD)                          
               RETURN 
            ELSE 
               XNEW = XOLD+DX 
               YNEW = YOLD+DY 
               DZMIN = BETA * (RHO + ABSCOM (XNEW, YNEW) ) 
               CALL COMHOR (XNEW, YNEW, N, GAMMA, UNEW, VNEW, UDNEW,    &
               VDNEW, UDDNEW, VDDNEW, BD, BDD)                          
               PBNEW = ABSCOM (UNEW, VNEW) 
               IF (PBNEW.EQ.0.0D0) THEN 
                  RETURN 
               ELSEIF (ABSCOM (DX, DY) .GT.DZMIN.AND.PBNEW.GT.BD) THEN 
                  GOTO 10 
               ELSE 
                  ENDIT = .TRUE. 
                  BDZE = ABSCOM (DX, DY) 
                  GOTO 10 
               ENDIF 
            ENDIF 
         ENDIF 
      ENDIF 
      END SUBROUTINE BAUZRO                         
!                                                                       
!                                                                       
      SUBROUTINE COMHOR (X, Y, N, GAMMA, U, V, UD, VD, UDD, VDD, BDP,   &
      BDPD)                                                             
!                                                                       
!*****************************************************************      
!                                                                *      
!  This SUBROUTINE calculates the complex functional value       *      
!  PN(X+I*Y)=(U+I*V), the complex valued derivatives             *      
!  PN'(X+I*Y)=(UD+I*VD) and PN''(X+I*Y)=(UDD+I*VDD) of a polymial*      
!  PN(Z) of degree N (N>0) with complex coefficients by using    *      
!  the Horner-scheme.                                            *      
!  Additionally bounds BDP and BDPD are computed for the         *      
!  rounding errors in computing DABS(PN(Z)) and DABS(PN'(Z)).    *      
!  The complex coefficients of PN are stored in a 2-dimensional  *      
!  array A(1:2*(N+1)), arranged in descending order of the       *      
!  powers (they will remain unchanged by this subroutine).       *      
!  The complex coefficients of the polynomial Q(Z) of degree N-1 *      
!  are stored in the array B(1:2*(N+1)), arranged in descending  *      
!  order. Here Q(Z) is defined by PN(Z)=Q(Z)*(Z-Z0)+PN(Z0).      *      
!  The array C(1:2*N) is used as an auxiliary array.             *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  X        : real part of the value for which the functional    *      
!             value and its 1st and 2nd derivatives are to be    *      
!             computed for the polynomial PN.                    *      
!  Y        : imaginary part of the value for which the          *      
!             functional value and its 1st and 2nd derivatives   *      
!             are to be computed for the polynomial PN.          *      
!  GAMMA    : error bound.                                       *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  U        : real part of PN(X+I*Y).                            *      
!  V        : imaginary part of PN(X+I*Y).                       *      
!                                                                *      
!  UD       : real part of PN'(X+I*Y).                           *      
!  VD       : imaginary part of PN'(X+I*Y).                      *      
!                                                                *      
!  UDD      : real part of PN''(X+I*Y).                          *      
!  VDD      : imaginary part of PN''(X+I*Y).                     *      
!                                                                *      
!  BDP      : bound for the rounding error of DABS(PN(Z)).       *      
!  BDPD     : bound for the rounding error of DABS(PN'(Z)).      *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: ABSCOM.                                 *      
!                                                                *      
!                                                                *      
!  sources: Bauhuber, see [BAUH70].                              *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Guido Dubois                                     *      
!  date       : 11.01.1985                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      COMMON / GLOBAL / A, B, C 
      DOUBLEPRECISION A (202), B (202), C (200), ABSCOM, GAMMA, X, Y, U,&
      V, UD, VD, UDD, VDD, BDP, BDPD, H, H1, H2, H3, H4                 
      C (1) = A (1) 
      B (1) = A (1) 
      C (2) = A (2) 
      B (2) = A (2) 
      BDPD = ABSCOM (A (1), A (2) ) 
      BDP = BDPD 
      MS = N - 1 
      M = N 
      J = N 
      NM2P1 = N * 2 + 1 
      DO 10 K = 3, NM2P1, 2 
         J = J - 1 
         H1 = X * B (K - 2) - Y * B (K - 1) 
         H2 = Y * B (K - 2) + X * B (K - 1) 
         B (K) = A (K) + H1 
         B (K + 1) = A (K + 1) + H2 
         H3 = ABSCOM (A (K), A (K + 1) ) 
         H4 = ABSCOM (H1, H2) 
         H = H3 
         IF (H3.LT.H4) H = H4 
         IF (H.GT.BDP) THEN 
            BDP = H 
            M = J 
         ENDIF 
         IF (K.EQ.NM2P1) THEN 
            GOTO 20 
         ELSE 
            H1 = X * C (K - 2) - Y * C (K - 1) 
            H2 = Y * C (K - 2) + X * C (K - 1) 
            C (K) = B (K) + H1 
            C (K + 1) = B (K + 1) + H2 
            H3 = ABSCOM (B (K), B (K + 1) ) 
            H4 = ABSCOM (H1, H2) 
            H = H3 
            IF (H3.LT.H4) H = H4 
            IF (BDPD.LT.H) THEN 
               BDPD = H 
               MS = J - 1 
            ENDIF 
         ENDIF 
   10 END DO 
   20 CONTINUE 
      U = B (2 * N + 1) 
      V = B (2 * N + 2) 
      UD = C (2 * N - 1) 
      VD = C (2 * N) 
      H = ABSCOM (X, Y) 
      IF (H.NE.0.0D0) THEN 
         BDP = BDP * FLOAT (M + 1) * H**M 
         BDPD = BDPD * FLOAT (MS + 1) * H**MS 
      ELSE 
         BDP = ABSCOM (U, V) 
         BDPD = ABSCOM (UD, VD) 
      ENDIF 
      BDP = BDP * GAMMA 
      BDPD = BDPD * GAMMA 
      IF (N.GT.1) THEN 
         H1 = C (1) 
         H2 = C (2) 
         NM2M3 = N * 2 - 3 
         DO 30 K = 3, NM2M3, 2 
            H = C (K) + X * H1 - Y * H2 
            H2 = C (K + 1) + Y * H1 + X * H2 
            H1 = H 
   30    END DO 
         UDD = 2.0D0 * H1 
         VDD = 2.0D0 * H2 
         RETURN 
      ELSE 
         UDD = 0.0D0 
         VDD = 0.0D0 
         RETURN 
      ENDIF 
      END SUBROUTINE COMHOR                         
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION ABSCOM (X, Y) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This FUNCTION-subroutine calculates the absolute value of a   *      
!  complex number (X+I*Y).                                       *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  X        : real part of the complex number.                   *      
!  Y        : imaginary part of the complex number.              *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETER:                                             *      
!  =================                                             *      
!  ABSCOM   : absolute value of the complex number.              *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none.                                   *      
!                                                                *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Guido Gubois                                     *      
!  date       : 11.01.1985                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION X, Y 
      IF (X.NE.0.0D0.OR.Y.NE.0.0D0) THEN 
         IF (DABS (X) .GE.DABS (Y) ) THEN 
            ABSCOM = DABS (X) * DSQRT (Y / X * Y / X + 1.0D0) 
            RETURN 
         ELSE 
            ABSCOM = DABS (Y) * DSQRT (X / Y * X / Y + 1.0D0) 
            RETURN 
         ENDIF 
      ELSE 
         ABSCOM = 0.0D0 
         RETURN 
      ENDIF 
      END FUNCTION ABSCOM                           
!                                                                       
!                                                                       
      SUBROUTINE MCONST (FMACHP, INFINY, SMALNO, BASE) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This subroutine sets up some constants that are machine       *      
!  dependent.                                                    *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  FMACHP   : machine constant for DOUBLE PRECISION.             *      
!  INFINY   : largest representable floating-point number.       *      
!  SMALNO   : smallest representable floating-point number.      *      
!  BASE     : base of the floating-point number system used to   *      
!             represent machine numbers.                         *      
!                                                                *      
!                                                                *      
!  Description of the auxiliary variables:                       *      
!  =======================================                       *      
!  I        : number of digits of the floating-point mantissa    *      
!             of DOUBLE PRECISION numbers.                       *      
!  M        : largest allowed exponent.                          *      
!  N        : smallest allowed exponent.                         *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none.                                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Guido Gubois                                     *      
!  date       : 11.01.1985                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION FMACHP, INFINY, SMALNO, BASE 
      BASE = 8.0D0 
      I = 29 
      M = 322 
      N = - 293 
      FMACHP = 1.0D0 
   10 FMACHP = 0.5D0 * FMACHP 
      IF (1.0D0.LT.1.0D0 + FMACHP) GOTO 10 
      FMACHP = 2.0D0 * FMACHP 
      INFINY = BASE * (1.0D0 - BASE** ( - I) ) * BASE** (M - 1) 
      SMALNO = (BASE** (N + 3) ) / BASE**3 
      RETURN 
      END SUBROUTINE MCONST                         
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION SCALFC (NN, PT, FMACHP, INFINY, SMALNO,  &
      BASE)                                                             
!                                                                       
!*****************************************************************      
!                                                                *      
!  This FUNCTION-subroutine calculates a scaling factor which    *      
!  is used to scale the polynomial coefficients.                 *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  NN       : 1 + the degree of the polynomial.                  *      
!  PT       : nn-vector PT(1:NN) containing the absolute         *      
!             values of the polynomial's coefficients.           *      
!  FMACHP   : machine constant for DOUBLE PRECISION.             *      
!  INFINY   : largest representable floating-point number.       *      
!  SMALNO   : smallest representable floating-point number.      *      
!  BASE     : base for the floating-point number system used by  *      
!             the machine.                                       *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETER:                                             *      
!  =================                                             *      
!  SCALFC   : scaling factor.                                    *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none.                                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Guido Dubois                                     *      
!  date       : 11.01.1985                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION PT (NN), FMACHP, INFINY, SMALNO, BASE, HI, LO,    &
      MAX, MIN, X, SC                                                   
      HI = DSQRT (INFINY) 
      LO = SMALNO / FMACHP 
      MAX = 0.0D0 
      MIN = INFINY 
      DO 10 I = 1, NN 
         X = PT (I) 
         IF (X.GT.MAX) MAX = X 
         IF (X.NE.0.0D0.AND.X.LT.MIN) MIN = X 
   10 END DO 
      SCALFC = 1.0D0 
      IF (MIN.GE.LO.AND.MAX.LE.HI) THEN 
         RETURN 
      ELSE 
         X = LO / MIN 
         IF (X.GT.1.0D0) THEN 
            SC = X 
            IF (INFINY / SC.GT.MAX) SC = 1.0D0 
         ELSE 
            SC = 1.0D0 / (DSQRT (MAX) * DSQRT (MIN) ) 
         ENDIF 
         L = DLOG (SC) / DLOG (BASE) + 0.5D0 
         SCALFC = BASE**L 
      ENDIF 
      RETURN 
      END FUNCTION SCALFC                           
!                                                                       
!                                                                       
      SUBROUTINE CDIV (A, B, C, D, X, Y) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This SUBROUTINE performs a complex division                   *      
!           (X+I*Y) := (A+I*B)/(C+I*D).                          *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  A        : real part of the numerator.                        *      
!  B        : imaginary part of the numerator.                   *      
!                                                                *      
!  C        : real part of the denominator.                      *      
!  D        : imaginary part of the denominator.                 *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  X        : real part of the quotient.                         *      
!  Y        : imaginary part of the quotient.                    *      
!                                                                *      
!                                                                *      
!  NOTE: If the denominator's real and imaginary parts are both  *      
!        equal to zero, the program is aborted with a detailed   *      
!        error message.                                          *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none.                                   *      
!                                                                *      
!                                                                *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author     : Guido Dubois                                     *      
!  date       : 11.01.1985                                       *      
!  source     : FORTRAN 77                                       *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION A, B, C, D, X, Y, U, V, AM, AN, P, Q, F 
      IF (C.NE.0.0D0.OR.D.NE.0.0D0) THEN 
         IF (A.NE.0.0D0.OR.B.NE.0.0D0) THEN 
            IF (DABS (A) .GT.DABS (B) ) THEN 
               U = A 
               AM = 1.0D0 
               AN = B / A 
            ELSE 
               U = B 
               AM = A / B 
               AN = 1.0D0 
            ENDIF 
            IF (DABS (C) .GT.DABS (D) ) THEN 
               V = C 
               P = 1.0D0 
               Q = D / C 
            ELSE 
               V = D 
               P = C / D 
               Q = 1.0D0 
            ENDIF 
            F = U / V 
            V = P * P + Q * Q 
            U = (AM * P + AN * Q) / V 
            X = U * F 
            U = ( - AM * Q + AN * P) / V 
            Y = U * F 
            RETURN 
         ELSE 
            X = 0.0D0 
            Y = 0.0D0 
            RETURN 
         ENDIF 
      ELSE 
         WRITE ( * , * ) 'DIVISION BY ZERO IN SUBROUTINE CDIV' 
         STOP 
      ENDIF 
      END SUBROUTINE CDIV                           
