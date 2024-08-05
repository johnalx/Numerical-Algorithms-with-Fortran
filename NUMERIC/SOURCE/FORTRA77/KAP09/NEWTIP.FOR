C[BA*)
C[KA{F 9}
C[  {Polynomial and Rational Interpolation}
C[  {Polynomial and Rational Interpolation}*)
C[FE{F 9.5.1}
C[  {Newton Formula for Arbitrary Nodes}
C[  {Newton Formula for Arbitrary Nodes}*)
C[LE*)
      SUBROUTINE NEWTIP(N,X,Y,B,IERR)
C[IX{NEWTIP}*)
C
C*****************************************************************
C                                                                *
C     NEWTIP determines the coefficients of the interpolation    *
C     polynomial in Newton's form.                               *
C[BE*)
C     Functional values of this polynomial can then be determined*
C     by FUNCTION NIPFCT.                                        *
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     N       : degree of the interpolation polynomial           *
C               (= index of the last node, if counted from 0 on) *
C     X       : (N+1)-vector X(0:N) ] value pairs (X(I),Y(I)) to *
C     Y       : (N+1)-vector Y(0:N) ] be interpolated by a       *
C                                   ] polynomial                 *
C                                                                *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     B       : (N+1)-vector B(0:N), the coefficients of the     *
C               polynomial in the following form:                *
C                  NP(X) = B(0) + B(1)*(X-X(0)) +                *
C                               + B(2)*(X-X(0))*(X-X(1)) + ...   *
C                           ... + B(N)*(X-X(0))*...*(X-X(N-1))   *
C     IERR    : error code                                       *
C               = 0 : no error                                   *
C               = 1 : N < 0                                      *
C               = 2 : there are two identical nodes              *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Elmar Pohl                                         *
C  date     : 09.28.1985                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION X(0:N),Y(0:N),B(0:N)
      IERR=0
      IF(N .LT. 0) THEN
         IERR=1
         RETURN
      ENDIF
      DO 10 I=0,N
         B(I)=Y(I)
   10 CONTINUE
      DO 20 I=1,N
         DO 30 K=N,I,-1
            H=X(K)-X(K-I)
            IF(H .EQ. 0.0D0) THEN
               IERR=2
               RETURN
            ENDIF
            B(K)=(B(K)-B(K-1))/H
   30    CONTINUE
   20 CONTINUE
      RETURN
      END
C
C
C[BA*)
C[LE*)
      DOUBLE PRECISION FUNCTION NIPFCT(X0,X,B,N)
C[IX{NIPFCT}*)
C
C*****************************************************************
C                                                                *
C     NIPFCT determines the functional value of the Newton       *
C     interpolation polynomial with the coefficients B(I) for    *
C     I=0,1, ...,N at X0 by a generalized Horner scheme.         *
C[BE*)
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     X0      : X-value where the polynomial is to be evaluated  *
C     X       : (N+1)-vector X(0:N) of nodes X(I), I=0,1,...,N,  *
C               for the interpolating polynomial                 *
C     B       : (N+1)-vector B(0:N), the coefficients of the     *
C               Newton interpolation polynomial, in the following*
C               form:
C                  NP(X) = B(0) + B(1)*(X-X(0)) +                *
C                               + B(2)*(X-X(0))*(X-X(1)) + ...   *
C                           ... + B(N)*(X-X(0))*...*(X-X(N-1))   *
C     N       : degree of the interpolating polynomial           *
C                                                                *
C                                                                *
C     OUTPUT PARAMETER:                                          *
C     =================                                          *
C     NIPFCT  : NP(X0), the value of the interpolating polynomial*
C               at X0                                            *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Elmar Pohl                                         *
C  date     : 09.28.1985                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION B(0:N),X(0:N)
      NIPFCT=B(N)
      DO 10 I=N-1,0,-1
         NIPFCT=NIPFCT*(X0-X(I))+B(I)
   10 CONTINUE
      RETURN
      END
