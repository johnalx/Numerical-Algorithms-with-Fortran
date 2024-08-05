      SUBROUTINE JACOBI (FX, X, M, N, DF, LDDF, EPS, WORK) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  The SUBROUTINE JACOBI determines the Jacobi matrix for a      *      
!  vector valued function made up of M+1 real valued functions   *      
!  in N+1 real variables at a point X.                           *      
!  The partial derivatives are approximated by central           *      
!  difference quotients.                                         *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  FX     : denotes a SUBROUTINE that has to be provided by      *      
!           the user. It evaluates the given function. In the    *      
!           calling program FX has to be defined as EXTERNAL.    *      
!           It has to be formatted as follows:                   *      
!               SUBROUTINE FX (X, N, F, M)                       *      
!               INTEGER M, N                                     *      
!               DOUBLE PRECISION X(0:N), F(0:M)                  *      
!               ------------------                               *      
!               F(0) = F0 (X(0), ... , X(N)                      *      
!                .                                               *      
!                .                                               *      
!               F(M) = FM (X(0), ... , X(N)                      *      
!               ------------------                               *      
!               RETURN                                           *      
!               END                                              *      
!  X      : (N+1)-vector X(0:N) containing the coordinates of    *      
!           the point, where the Jacobi matrix of FX is to be    *      
!           determined                                           *      
!  M      : M+1 = number of component functions of FX            *      
!  N      : N+1 = number of variables                            *      
!  DF     : 2-dimensional array DF(0:LDDF, 0:N); storage for the *      
!           Jacobi matrix                                        *      
!  LDDF   : leading dimension of the array DF as defined in the  *      
!           calling program                                      *      
!  EPS    : indicates the precision with which the partial       *      
!           derivatives are to be computed                       *      
!                                                                *      
!                                                                *      
!  AUXILIARY PARAMETER:                                          *      
!  ====================                                          *      
!  WORK   : (N+1)-vector WORK(0:M) for intermediate storage      *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETER:                                             *      
!  =================                                             *      
!  DF     : 2-dimensional array DF(0:LDDF,0:N) containing the    *      
!           approximate Jacobi matrix                            *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Ilona Westermann                                   *      
!  date     : 01.09.1987                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      INTEGER M, N, LDDF 
      DOUBLEPRECISION X (0:N), DF (0:LDDF, 0:N), EPS, WORK (0:M) 
      FACTOR = EPS** (1.0D0 / 3.0D0) 
      DO 20 K = 0, N 
         XK = X (K) 
         IF (XK.EQ.0.0D0) THEN 
            HK = FACTOR 
         ELSE 
            HK = FACTOR * DABS (XK) 
         ENDIF 
         ZHK = 1.0D0 / (2.0D0 * HK) 
         X (K) = XK + HK 
         CALL FX (X, N, DF (0, K), M) 
         X (K) = XK - HK 
         CALL FX (X, N, WORK, M) 
         DO 10 I = 0, M 
            DF (I, K) = (DF (I, K) - WORK (I) ) * ZHK 
   10    END DO 
         X (K) = XK 
   20 END DO 
      RETURN 
      END SUBROUTINE JACOBI                         
