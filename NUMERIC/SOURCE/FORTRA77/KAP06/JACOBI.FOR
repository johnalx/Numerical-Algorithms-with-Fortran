C[BA*)
C[LE*)
      SUBROUTINE JACOBI (FX, X, M, N, DF, LDDF, EPS, WORK)
C[IX{JACOBI}*)
C
C*****************************************************************
C                                                                *
C  The SUBROUTINE JACOBI determines the Jacobi matrix for a      *
C  vector valued function made up of M+1 real valued functions   *
C  in N+1 real variables at a point X.                           *
C  The partial derivatives are approximated by central           *
C  difference quotients.                                         *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  FX     : denotes a SUBROUTINE that has to be provided by      *
C           the user. It evaluates the given function. In the    *
C           calling program FX has to be defined as EXTERNAL.    *
C           It has to be formatted as follows:                   *
C               SUBROUTINE FX (X, N, F, M)                       *
C               INTEGER M, N                                     *
C               DOUBLE PRECISION X(0:N), F(0:M)                  *
C               ------------------                               *
C               F(0) = F0 (X(0), ... , X(N)                      *
C                .                                               *
C                .                                               *
C               F(M) = FM (X(0), ... , X(N)                      *
C               ------------------                               *
C               RETURN                                           *
C               END                                              *
C  X      : (N+1)-vector X(0:N) containing the coordinates of    *
C           the point, where the Jacobi matrix of FX is to be    *
C           determined                                           *
C  M      : M+1 = number of component functions of FX            *
C  N      : N+1 = number of variables                            *
C  DF     : 2-dimensional array DF(0:LDDF, 0:N); storage for the *
C           Jacobi matrix                                        *
C  LDDF   : leading dimension of the array DF as defined in the  *
C           calling program                                      *
C  EPS    : indicates the precision with which the partial       *
C           derivatives are to be computed                       *
C                                                                *
C                                                                *
C  AUXILIARY PARAMETER:                                          *
C  ====================                                          *
C  WORK   : (N+1)-vector WORK(0:M) for intermediate storage      *
C                                                                *
C                                                                *
C  OUTPUT PARAMETER:                                             *
C  =================                                             *
C  DF     : 2-dimensional array DF(0:LDDF,0:N) containing the    *
C           approximate Jacobi matrix                            *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Ilona Westermann                                   *
C  date     : 01.09.1987                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER M, N, LDDF
      DOUBLE PRECISION X(0:N), DF(0:LDDF,0:N), EPS, WORK(0:M)
      FACTOR = EPS ** (1.0D0/3.0D0)
      DO 20 K=0,N
         XK = X(K)
         IF (XK .EQ. 0.0D0) THEN
            HK = FACTOR
         ELSE
            HK = FACTOR * DABS(XK)
         ENDIF
         ZHK = 1.0D0/(2.0D0*HK)
         X(K) = XK + HK
         CALL FX (X, N, DF(0,K), M)
         X(K) = XK - HK
         CALL FX (X, N, WORK, M)
         DO 10 I=0,M
            DF(I,K) = (DF(I,K) - WORK(I)) * ZHK
   10    CONTINUE
         X(K) = XK
   20 CONTINUE
      RETURN
      END
