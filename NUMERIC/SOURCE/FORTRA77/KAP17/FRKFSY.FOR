C[BA*)
C[LE*)
      SUBROUTINE FRKFSY(A,DA,N,Y,DES,H,HMX,ABSERR,RELERR,IERR)
C[IX{FRKFSY}*)
C
C*****************************************************************
C                                                                *
C  A system of ordinary differential equations of 1st order is   *
C  integrated by applying the RUNGE-KUTTA-FEHLBERG method        *
C  [ of order O(H**5) ] with estimates for the local error and   *
C  step size control.                                            *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  A      : starting value for the integration interval          *
C  DA     : length of the integration interval;                  *
C           DA may be < 0.0 if we want to integrate to the left. *
C  N      : number of equations; N < 11.                         *
C  Y      : vector Y(1:N); the initial values at A               *
C  DES    : SUBROUTINE, that describes the system of differential*
C           equations, given in the following form:              *
C                      SUBROUTINE  DES(X,Y,YS)                   *
C                      X : independent variable                  *
C                      Y : vector of dependent variables         *
C                      YS: vector YS(I)=DY(I)/DX of derivatives  *
C                          at X, I=1,...,N                       *
C                      Y and YS are dimensioned as DOUBLE        *
C                      PRECISION Y(1), YS(1), however, they may  *
C                      be to be used as a vector of length N.    *
C           example :  SUBROUTINE DES(X,Y,YS)                    *
C                      DOUBLE PRECISION Y(1),YS(1)               *
C                      YS(1)=Y(1)                                *
C                      YS(2)=-Y(2)                               *
C                      RETURN                                    *
C                      END                                       *
C  H      : initial step size; if H is provided unrealistically, *
C           H is modified internally; H may be negative if       *
C           DA < 0.0.                                            *
C  HMX    : upper bound for the step size magnitude used during  *
C           calculation. HMX > 0.0                               *
C  ABSERR :] bounds for the acceptable local error, relative to  *
C  RELERR :] the current step size. If the following holds for   *
C         :] each component of the computed solution Y(I)        *
C                  ABS ( estimate of the local error) .LE.       *
C                      ABS(H)*(RELERR*ABS(Y(I))+ABSERR),         *
C           then the solution is accepted in the current step.   *
C           If ABSERR = 0.0, we test for the relative error;     *
C           If RELERR = 0.0, we test for the absolute error.     *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  A      : last x value for which a solution was successfully   *
C           determined. Normally the following will hold:        *
C           A on output = A on input + DA.                       *
C  Y      : computed solution vector at  A on output             *
C  H      : optimal step size, which was used for the last step. *
C  IERR   : = 1, everything o.k.; solution found at A + DA.      *
C           = 2, after 3000 calls of SUBROUTINE DES we stop with-*
C                out having reached the endpoint A+DA. If com-   *
C                putations are to be continued, call FRKFSY again*
C                with unchanged parameters.                      *
C           = 3, false input data; i.e.                          *
C                ABSERR.LT.0.0     or    RELERR.LT.0.0     or    *
C                ABSERR + RELERR = 0.0  or  HMX.LE.0.0: Return.  *
C           = 4, the optimal step size cannot be achieved for the*
C                computer.     RETURN                            *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required : MACHPD                                 *
C                                                                *
C                                                                *
C  sources : SHAMPINE/ALLEN, see [SHAM73].                       *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Richard Reuter                                     *
C  date     : 02.09.1983                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION YT(10),T(10),R(10),K1(10),K2(10)
      DOUBLE PRECISION K3(10),K4(10),K5(10),K6(10)
      DOUBLE PRECISION Y(N)
C
C     determine machine constant
C
      FMACHP = 1.0D0
    2 FMACHP = 0.5D0 * FMACHP
      IF (MACHPD(1.0D0+FMACHP) .EQ. 1) GOTO 2
      FMACHP = FMACHP * 2.0D0
C
C     check the input data
C
      IERR=3
      IF(RELERR .LT. 0.0D0 .OR. ABSERR .LT. 0.0D0 .OR.
     1   RELERR+ABSERR .EQ. 0.0D0 .OR. HMX .LE. 0.0D0) RETURN
C
      IERR=4
      B=A+DA
      IF(DABS(DA) .LE. 13.0D0*FMACHP*DMAX1(DABS(A),DABS(B))) RETURN
C
      HMAX=DMIN1(HMX,DABS(DA))
      IF(DABS(H) .LE. 13.0D0*FMACHP*DABS(A)) H=HMAX
C
C     Initialize counter for calls of SUBROUTINE DES
C
      LFD=0
      IAD=0
C
C     H is bounded by HMAX and is chosen so that the
C     endpoint B is reached, if possible.
C
   3  H=DSIGN(DMIN1(DABS(H),HMAX),DA)
      IF(DABS(B-A) .LE. 1.25D0*DABS(H)) THEN
         HF=H
C
C        if IAD=1 and H=B-A acceptable, we stop after
C        the next integration step.
C
         IAD=1
         H=B-A
      END IF
C
C     an integration step is executed
C
      CALL DES(A,Y,K1)
      LFD=LFD+1
   5  CONTINUE
      X=0.25D0*H
      DO 6 I=1,N
         YT(I)=Y(I)+X*K1(I)
   6  CONTINUE
      X=A+X
      CALL DES(X,YT,K2)
      DO 7 I=1,N
         YT(I)=Y(I)+H*(K1(I)*(3.0D0/32.0D0)+K2(I)*(9.0D0/32.0D0))
   7  CONTINUE
      X=A+H*(3.0D0/8.0D0)
      CALL DES(X,YT,K3)
      DO 8 I=1,N
         YT(I)=Y(I)+H*(K1(I)*(1932.0D0/2197.0D0)
     1                -K2(I)*(7200.0D0/2197.0D0)
     1                +K3(I)*(7296.0D0/2197.0D0))
   8  CONTINUE
      X=A+H*(12.0D0/13.0D0)
      CALL DES(X,YT,K4)
      DO 9 I=1,N
         YT(I)=Y(I)+H*(K1(I)*(439.0D0/216.0D0)-8.0D0*K2(I)
     1                +K3(I)*(3680.0D0/513.0D0)
     1                -K4(I)*(845.0D0/4104.0D0))
   9  CONTINUE
      X=A+H
      CALL DES(X,YT,K5)
      DO 10 I=1,N
         YT(I)=Y(I)+H*(-K1(I)*(8.0D0/27.0D0)+2.0D0*K2(I)
     1                 -K3(I)*(3544.0D0/2565.0D0)
     1                 +K4(I)*(1859.0D0/4104.0D0)
     1                 -K5(I)*(11.0D0/40.0D0))
  10     CONTINUE
      X=A+0.5D0*H
      CALL DES(X,YT,K6)
      DO 11 I=1,N
         T(I)=K1(I)*(25.0D0/216.0D0)+K3(I)*(1408.0D0/2565.0D0)
     1         +K4(I)*(2197.0D0/4104.0D0)-K5(I)*0.20D0
         YT(I)=Y(I)+H*T(I)
  11  CONTINUE
C
C     YT(I) now represents the latest result of this pass.
C     Determine R(I), the estimate of the local
C     error, relative to the current step size.
C
      DO 12 I=1,N
         R(I)=K1(I)/360.0D0-K3(I)*(128.0D0/4275.0D0)
     1         -K4(I)*(2197.0D0/75240.0D0)+K5(I)/50.0D0
     1         +K6(I)*(2.0D0/55.0D0)
  12  CONTINUE
C
C     Check accuracy
C
      QUOT=0.0D0
      DO 13 I=1,N
         TR=DABS(R(I))/(RELERR*DABS(YT(I))+ABSERR)
         QUOT=DMAX1(QUOT,TR)
  13  CONTINUE
C
C     If  QOUOT.LE.1.0   ==> integration step is accepted
C
      IF(QUOT .LE. 1.0D0) THEN
C
C        result is accepted
C
         DO 14 I=1,N
            Y(I)=YT(I)
  14     CONTINUE
         A=A+H
C
C        if A=B ,  RETURN
C
         IF(IAD .EQ. 1) THEN
            IERR=1
            H=HF
            RETURN
         END IF
C
C        prepare next step
C
         QUOT=DMAX1(QUOT,6.5536D-4)
      END IF
      QUOT=DMIN1(QUOT,4096.0D0)
      H=0.8D0*H/DSQRT(DSQRT(QUOT))
C
C     We just achieved that H was increased by at most a factor of 5,
C     or alternatively, that it was decreased by a factor of 10
C     at most
C
      IF(DABS(H) .LE. 13.0D0*FMACHP*DABS(A)) THEN
         IERR=4
         RETURN
      END IF
      LFD=LFD+5
      IF(LFD .GE. 2995) THEN
         IERR=2
         RETURN
      END IF
      IF(QUOT .LE. 1.0D0) THEN
C
C        the step was successful. Continue with another step.
C
         GOTO 3
      ELSE
C
C        the step is repeated for a smaller H.
C
         IAD=0
         GOTO 5
      END IF
      END
