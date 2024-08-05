C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 14.4}
C[  {Differentiation by the Romberg Method}
C[  {Differentiation by the Romberg Method}*)
C[LE*)
      SUBROUTINE DIFROM (FCT,X0,EPS,N,H,RES,ERREST,NEND,
C[IX{DIFROM}*)
     +                   HEND,IERR,D)
C
C*****************************************************************
C                                                                *
C  The subroutine DIFROM approximately computes the first        *
C  derivative of a given function FCT at X0 by using the         *
C  ROMBERG-method.                                               *
C[BE*)
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  FCT      : name of the function that is to be differentiated. *
C             It must be provided by the user in the form:       *
C                  DOUBLE PRECISION FUNCTION FCT(X)              *
C                    ...                                         *
C                  RETURN                                        *
C                  END                                           *
C             in the calling program. It has to be defined as    *
C             EXTERNAL.                                          *
C  X0       : location where the first derivative is wanted.     *
C  EPS      : desired accuracy for the derivative.               *
C  N        : maximum number of rows or columns of the ROMBERG   *
C             scheme - 1 ; N has to be > 0.                      *
C  H        : initial step size; H has to be >= 4 * machine      *
C             constant.                                          *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  RES      : approximate value for the first derivative of the  *
C             function FCT at X0.                                *
C  ERREST   : error estimate for the approximate value RES.      *
C  NEND     : number of rows or columns of the ROMBERG scheme    *
C             that was actually used.                            *
C  HEND     : terminal step size.                                *
C  EPS      : error bound actually used.                         *
C  IERR     : error parameter.                                   *
C             = 0, no error.                                     *
C             = 1, error in the input data.                      *
C             = 2, desired accuracy was not reached after N steps*
C                  ERREST >= EPS.                                *
C             = 3, step size for the difference quotient became  *
C                  too small ( < 4 * machine constant).          *
C                                                                *
C                                                                *
C  AUXILIARY  PARAMETERS:                                        *
C  ======================                                        *
C  D        : (N+1)-vector D(0:N); storage of the current row    *
C             in the ROMBERG scheme.                             *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: MACHPD                                  *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Elmar Pohl                                         *
C  editor   : Guido Dubois                                       *
C  date     : 04.25.88                                           *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DOUBLE PRECISION D(0:N)
C
C  local storage of the minimal error bound EPSMIN in
C  case that the subroutine is called repeatedly
C
      SAVE EPSMIN,IFLAG
      DATA IFLAG/0/
C
C  determine the machine constant and initialize the
C  minimal error bound
C
      IF(IFLAG .EQ. 0) THEN
         IFLAG=1
         FMACHP=1.0D0
   10    FMACHP=0.5D0*FMACHP
         IF(MACHPD(1.0D0+FMACHP) .EQ. 1) GOTO 10
         EPSMIN=8.0D0*FMACHP
      END IF
C
C  test the input data
C
      IF(EPS .LT. EPSMIN) EPS=EPSMIN
      IERR=1
      IF(N .LT. 1 .OR. H .LT. EPSMIN) RETURN
      IERR=0
C
C  determine the first central difference quotient
C
      H2=2.0D0*H
      H0=H
      D(0)=(FCT(X0+H0)-FCT(X0-H0))/H2
      K=0
   30 K=K+1
      D(K)=0.0D0
      D0=D(0)
      H2=H0
      H0=0.5D0*H0
      IF(H0 .LT. EPSMIN) THEN
         IERR=3
      ELSE
C
C  determine the central difference quotient
C
         D(0)=(FCT(X0+H0)-FCT(X0-H0))/H2
         IK=1
C
C  determine the linear combinations
C
         DO 20 J=1,K
            IK=IK*4
            D1=D(J)
            D(J)=(IK*D(J-1)-D0)/(IK-1)
            D0=D1
   20    CONTINUE
C
C  check the break-off criteria
C
         ERREST=DABS(D(K)-D(K-1))
         IF(ERREST .GE. EPS) THEN
            IF(K .LE. N) GOTO 30
            IERR=2
         END IF
      END IF
      RES=D(K)
      NEND=K+1
      HEND=H0
      RETURN
      END
