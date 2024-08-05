      SUBROUTINE LAGUER (A, N, ABSERR, RELERR, MAXIT, XI, NITER, INUM,  &
      WORK, IERR)                                                       
!                                                                       
!*****************************************************************      
!                                                                *      
!  For a real polynomial PN of degree N, this program computes   *      
!  all real roots using Laguerre's method, provided PN has only  *      
!  real roots.                                                   *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  A      : (N+1)-vector A(0:N) containing the coefficients      *      
!           of the real polynomial PN, where                     *      
!              PN(X) = A(0) + A(1)*X + ... + A(N)*X**N           *      
!  N      : degree of PN, N > 2                                  *      
!  ABSERR : ) error bounds, each of which must be nonnegative,   *      
!  RELERR : ) with their sum positive. The following break-off   *      
!             criterion is used:                                 *      
!                ABS(X1-X2) <= ABS(X2)*RELERR+ABSERR.            *      
!             If RELERR=0.0, then we test for the absolute error,*      
!             if ABSERR=0.0, then we test the relative error.    *      
!             The input values for  ABSERR and RELERR  are used  *      
!             without modification only if both exceed four times*      
!             the machine constant or , in case one is zero, the *      
!             other must exceed that constant. Otherwise on or   *      
!             both of the bounds are adjusted internally to that *      
!             value.                                             *      
!  IMAX   :   maximal number of iterations allowed for each zero,*      
!             IMAX >= 1                                          *      
!                                                                *      
!                                                                *      
!  AUXILIARY VARIABLES:                                          *      
!  ====================                                          *      
!  WORK   : (N+1)-vector WORK(0:N)                               *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  XI     : N-vector XI(1:N) containing the N real roots of PN   *      
!  NITER  : INTEGER N-vector NITER(1:N), containing the number of*      
!           iterations performed for finding the root with the   *      
!           same index                                           *      
!  IANZ   : number of roots found                                *      
!  IERR   : error parameter                                      *      
!           =0, invalid input values for ABSERR, RELERR, IMAX    *      
!               or N                                             *      
!           =1, all roots have been found                        *      
!           =2, maximally allowed number of iterations IMAX has  *      
!               been reached; the roots found earlier are stired *      
!               in XI(1:IANZ)                                    *      
!           =3, the intermediate variable S is negative, i.e.,   *      
!               SQRT(S) is not real and there may be complex     *      
!               roots of PN                                      *      
!           =4, error message from subroutine QUADRA when        *      
!               computing the last two roots of PN               *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: HORN1, HORN2, QUADRA, MACHPD            *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author    : Gisela Engeln-MÅllges                             *      
!  Date      : 06.03.1992                                        *      
!  Source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
!  Declarations                                                         
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A (0:N), WORK (0:N), XI (1:N) 
      INTEGER NITER (1:N) 
!                                                                       
!  check input values of  ABSERR, RELERR, IMAX and N                    
!                                                                       
      IF (ABSERR.LT.0.0D0.OR.RELERR.LT.0.0D0.OR.ABSERR +                &
      RELERR.LE.0.0D0.OR.MAXIT.LT.1.OR.N.LE.2) THEN                     
         IERR = 0 
         RETURN 
      ENDIF 
!                                                                       
!  Compute the machine constant and adjust error bounds if needed       
!                                                                       
      FMACHP = 1.0D0 
   10 FMACHP = 0.5D0 * FMACHP 
      IF (MACHPD (1.0D0 + FMACHP) .EQ.1) GOTO 10 
      FMACHP = 2.0D0 * FMACHP 
      EPS = 4.0D0 * FMACHP 
      IF (RELERR.EQ.0.0D0) THEN 
         IF (ABSERR.LT.EPS) ABSERR = EPS 
      ELSEIF (ABSERR.EQ.0.0D0) THEN 
         IF (RELERR.LT.EPS) RELERR = EPS 
      ELSE 
         IF (ABSERR.LT.EPS) ABSERR = EPS 
         IF (RELERR.LT.EPS) RELERR = EPS 
      ENDIF 
!                                                                       
!   Initialize                                                          
!                                                                       
      INUM = 0 
      NACT = N 
      IERR = 1 
!                                                                       
!  Loop to compute N-2 roots                                            
!                                                                       
      DO 20 I = 1, N - 2 
!                                                                       
!  Initialize iteration counter for the I-th root                       
!                                                                       
         ITER = 0 
!                                                                       
!  Initialize starting value for all roots                              
!                                                                       
         X = 0.0D0 
!                                                                       
!  Iteration loop                                                       
!                                                                       
   30    IF (ITER.GE.MAXIT) THEN 
            IERR = 2 
            RETURN 
         ENDIF 
         CALL HORN2 (A, N, X, NACT, PN, PN1, PN2, WORK) 
         ITER = ITER + 1 
!                                                                       
!  Calculate S from the iteration rule and find SQRT(S)                 
!                                                                       
         S = (NACT - 1) * ( (NACT - 1) * PN1**2 - NACT * PN * PN2) 
         IF (DABS (S) .LE.EPS) S = 0.0D0 
         IF (S.LT.0.0D0) THEN 
            IERR = 3 
            RETURN 
         ENDIF 
         SRS = DSQRT (S) 
!                                                                       
!  Choose the sign of  SQRT(S) so that it coincides                     
!  with the sign of PN1                                                 
!                                                                       
         SRS = DSIGN (SRS, PN1) 
!                                                                       
!  Calculate the denominator of the iteration rule.                     
!  If it is less than four times the machine constant,                  
!  we set it equal to  EPS                                              
!                                                                       
         XNENN = PN1 + SRS 
         IF (DABS (XNENN) .LT.EPS) THEN 
            XNENN = DSIGN (EPS, XNENN) 
         ENDIF 
!                                                                       
!  Find the difference of successive iterates and a new iterate         
!                                                                       
         DIFF = NACT * PN / XNENN 
         X = X - DIFF 
!                                                                       
!  check the break-off criterion                                        
!                                                                       
         IF (DABS (DIFF) .GT.DABS (X) * RELERR + ABSERR) GOTO 30 
         XI (I) = X 
         NITER (I) = ITER 
         INUM = INUM + 1 
         CALL HORN1 (A, N, X, NACT) 
         NACT = NACT - 1 
   20 END DO 
!                                                                       
!  Finally compute the remaining two roots for the                      
!  quadratic equation:                                                  
!  A(N)*X**2 + A(N-1)*X + A(N-2) = 0                                    
!                                                                       
      CALL QUADRA (A (N), A (N - 1), A (N - 2), EPS, X1, X2, IERR1) 
      IF (IERR1.NE.1) THEN 
         IERR = 4 
         RETURN 
      ENDIF 
      INUM = INUM + 2 
      XI (N - 1) = X1 
      XI (N) = X2 
      NITER (N - 1) = 0 
      NITER (N) = 0 
      RETURN 
      END SUBROUTINE LAGUER                         
!                                                                       
!                                                                       
      SUBROUTINE HORN1 (A, N, X0, NACT) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  If P has degree NACT, this subroutine finds the remainder     *      
!  polynomial PAB of degree NACT-1 that results from dividing off*      
!  one known root of P. It is neeeded in the SUBROUTINE LAGUER.  *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  A      : (N+1)-vector A(0:N) containing the coefficients of P:*      
!           P(X) = A(N)*X**NACT + A(N-1)*X**(NACT-1) + ...       *      
!                               + A(N-NACT+1)*X + A(N-NACT)      *      
!  N      : dimension of A + 1, as specified in the calling      *      
!           program                                              *      
!  X0     : known root of P; P is divided by (X-X0) to obtain PAB*      
!  NACT   : actual degree of P                                   *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  A      : (N+1)-vector A(0:N), containing the coefficients of  *      
!           the rest polynomial PAB                              *      
!           PAB(X) = A(N)*X**(NACT-1) + A(N-1)*X**(NACT-2) + ... *      
!                     + A(N-NACT+2)*X + A(N-NACT-1)              *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author    : Gisela Engeln-MÅllges                             *      
!  Date      : 06.03.1992                                        *      
!  Source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A (0:N) 
      HELP = 0.0D0 
      DO 10 K = N, N - NACT, - 1 
         HELP = HELP * X0 + A (K) 
         A (K) = HELP 
   10 END DO 
      RETURN 
      END SUBROUTINE HORN1                          
!                                                                       
!                                                                       
      SUBROUTINE HORN2 (A, N, X, NACT, P, P1, P2, WORK) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This subroutine calculates the function value and that of the *      
!  first and second derivative of a given polynom P of degree    *      
!  NACT at X. It is needed in SUBROUTINE LAGUER.                 *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  A      : (N+1)-vector A(0:N) containing the coefficients of P:*      
!           P(X) = A(N)*X**NACT + A(N-1)*X**(NACT-1) + ...       *      
!                               + A(N-NACT+1)*X + A(N-NACT)      *      
!  N      : dimension of A + 1, as specified in the calling      *      
!           program                                              *      
!  X      : value where P(X), P'(X) and P''(X) has to be computed*      
!  NACT   : actual degree of P                                   *      
!                                                                *      
!                                                                *      
!  AUXILIARY PARAMETERS:                                         *      
!  =====================                                         *      
!  WORK   : (N+1)-vector WORK(0:N)                               *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  P      : P(X)                                                 *      
!  P1     : P'(X)                                                *      
!  P2     : P''(X)                                               *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author    : Gisela Engeln-MÅllges                             *      
!  Date      : 02.18.1992                                        *      
!  Source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A (0:N), WORK (0:N) 
      DO 10 I = N - NACT, N 
         WORK (I) = A (I) 
   10 END DO 
      DO 20 K = 0, 2 
         HELP = 0.0D0 
         DO 30 I = N, N - NACT + K, - 1 
            HELP = HELP * X + WORK (I) 
            WORK (I) = HELP 
   30    END DO 
   20 END DO 
      P = WORK (N - NACT) 
      P1 = WORK (N - NACT + 1) 
      P2 = 2.0D0 * WORK (N - NACT + 2) 
      RETURN 
      END SUBROUTINE HORN2                          
!                                                                       
!                                                                       
      SUBROUTINE QUADRA (A, B, C, EPS, X1, X2, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This subroutine solves a real quadratic equation              *      
!        A*X**2 + B*X + C = 0  with A not equal to 0,            *      
!  that has real roots. It serves the subroutine LAGUER.         *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  A, B, C: real coefficients of the quadratic equation          *      
!  EPS    : EPS = 4.0 * machine constant                         *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  X1, X2 : real roots                                           *      
!  IERR   : error parameter                                      *      
!           =1, all is ok                                        *      
!           =2, error: either A=0 or complex roots               *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Authors   : Gisela Engeln-MÅllges                             *      
!  Date      : 09.03.1992                                        *      
!  Source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      IERR = 1 
      IF (A.EQ.0.0D0) THEN 
         IERR = 2 
         RETURN 
      ENDIF 
      D = B * B - 4.0D0 * A * C 
!                                                                       
!  If the magnitude of the discriminant  D = B**2 - 4*A*C  is less than 
!  EPS = 4.0 * machine constant, we set D = 0. If D <  -EPS, we stop:   
!  the roots are complex conjugate.                                     
!                                                                       
      IF (DABS (D) .LE.EPS) D = 0.0D0 
      IF (D.LT.0.0D0) THEN 
         IERR = 2 
         RETURN 
      ENDIF 
      WURZ = DSQRT (D) 
!                                                                       
!  Berechnung der beiden reellen Nullstellen der quadratischen          
!  Gleichung                                                            
!                                                                       
      V1 = B + WURZ 
      V2 = B - WURZ 
      IF (DABS (V1) .GE.DABS (V2) ) THEN 
         X1 = - 2.0D0 * C / V1 
      ELSE 
         X1 = - 2.0D0 * C / V2 
      ENDIF 
      IF (DABS (X1) .LE.EPS) THEN 
         X1 = 0.0D0 
         X2 = - (B / A) 
      ELSE 
         X2 = C / (A * X1) 
      ENDIF 
      RETURN 
      END SUBROUTINE QUADRA                         
