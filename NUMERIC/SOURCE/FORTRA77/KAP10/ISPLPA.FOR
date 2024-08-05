C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 10.1.3}
C[  {Parametric Cubic Splines}
C[  {Parametric Cubic Splines}*)
C[LE*)
      SUBROUTINE ISPLPA (N,XN,FN,T,MT,IB,ALPHA,BETA,BX,CX,DX,
C[IX{ISPLPA}*)
     +                   BY,CY,DY,DUMMY,IERR)
C
C*****************************************************************
C                                                                *
C  ISPLPA computes the coefficients BX(I), CX(I), DX(I), BY(I),  *
C  CY(I), DY(I) for I=0,1,...,N-1 of a parametric cubic inter-   *
C  polating spline for various end point conditions.             *
C[BE*)
C  The end point conditions are specified via the parameter IB.  *
C  The parametric two-component function depends on the parameter*
C  T(I) for  I=0,1,...,N and has two components SX and SY of the *
C  the following form:                                           *
C                                                                *
C  SX := SX(T) = XN(I) + BX(I)(T-T(I)) + CX(I)(T-T(I))**2 +      *
C                                      + DX(I)(T-T(I))**3        *
C                                                                *
C  SY := SY(T) = FN(I) + BY(I)(T-T(I)) + CY(I)(T-T(I))**2 +      *
C                                      + DY(I)(T-T(I))**3        *
C                                                                *
C  for T in the interval [T(I),T(I+1)], I=0,1,...N-1.            *
C                                                                *
C  SX and SY each are nonparametric cubic splines.               *
C                                                                *
C                                                                *
C  ASSUMPTIONS:    1.         N > 2                              *
C  ============    2.      T(I) < T(I+1), I=0,1,...,N-1          *
C                  3.     FN(0) = FN(N) , for IB = 4             *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  N  :  Index of the last node                                  *
C  XN :  vector XN(0:N); the nodes XN(I), I = 0,1,..,N           *
C  FN :  vector FN(0:N); the functional values FN(I) = FN(XN(I)) *
C  T  :  vector T(0:N); paramete value for XN(I) and FN(I)       *
C  MT :  indicates origin of the curve parameter T(I):           *
C        MT = 0 :  The user specifies the parameter values T(I), *
C                  I=0,1,...,N                                   *
C        MT = 1 :  Parameter values not prescribed. Their values *
C                  are computed in SUBROUTINE PSPPV from the     *
C                  chordal distances.                            *
C        MT = 2 :  Parameter values not prescribed. Their values *
C                  are computed in SUBROUTINE PSPPV from the     *
C                  arc length.                                   *
C                                                                *
C  IB : describes the end point conditions                       *
C       IB = 1: first derivative w.r.t. parameter is given       *
C       IB = 2: second derivative w.r.t. parameter is given      *
C       IB = 3: first derivative DY/DX is given                  *
C       IB = 4: periodic spline                                  *
C                                                                *
C  ALPHA : ] 2-vectors ..(1:2)                                   *
C  BETA  : ]                                                     *
C              if IB = 1 : ] IBth derivative w.r.t. parameter    *
C                 IB = 2 : ]       ALPHA(1)=SX(IB)(T(0))         *
C                                  ALPHA(2)=SY(IB)(T(0))         *
C                                  BETA(1) =SX(IB)(T(N))         *
C                                  BETA(2) =SY(IB)(T(N))         *
C              if IB = 3 :  first derivative DY/DX               *
C                           ALPHA(1) = DY/DX (XN(0))             *
C                           BETA(1)  = DY/DX (XN(N))             *
C                           ALPHA(2) : not used                  *
C                           BETA(2)  : not used                  *
C               If the magnitude of ALPHA(1) or BETA(1) exceeds  *
C               1.E10, then we compute the corresponding tangent *
C               vector as follows:                               *
C                . .                                             *
C               (X,Y) = (0,DSIGN(1,FN(1)-FN(0))  (left end point)*
C                . .                                             *
C               (X,Y) = (0,DSIGN(1,FN(N)-FN(N-1))                *
C                                               (right end point)*
C                                                                *
C              if IB = 4 : not used                              *
C                                                                *
C  (A natural parametric interpolating spline is obtained for    *
C   IB=2 and ALPHA(1)=ALPHA(2)=BETA(1)=BETA(2)=0.0)              *
C                                                                *
C                                                                *
C  AUXILIARY VBARIABLE:                                          *
C  ====================                                          *
C  DUMMY :  vector DUMMY(1:5*N+1)                                *
C                                                                *
C                                                                *
C  AUSGABEPARAMETER:                                             *
C  =================                                             *
C  XN :  ] N+1-vectors ..(0:N);                                  *
C  BX :  ] their elements 0, ..., N-1 describe the component     *
C  CX :  ] function SX.                                          *
C  DX :  ]                                                       *
C                                                                *
C  FN :  ] N+1-vectors ..(0:N);                                  *
C  BY :  ] their elements 0, ..., N-1 describe the component     *
C  CY :  ] function SY.                                          *
C  DY :  ]                                                       *
C          The elements BX(N), CX(N), DX(N), BY(N), CY(N) and    *
C          DY(N) are used for auxiliary purposes.                *
C  IERR :  error parameter                                       *
C          =  0 :  All is ok                                     *
C          = -1 :  N < 3                                         *
C          = -2 :  IB < 1  or  IB > 4                            *
C          = -3 :  Parameter values T(I) non monotonic:          *
C                  T(I) > or = T(I+1) for some I=0,1,...,N-1     *
C          = -4 :  For IB = 4 we have FN(0) not equal to FN(N)   *
C          =  1 :  crash in TRDSY or CYTSY, system matrix numer- *
C                  ically singular.                              *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Subroutines required: ISPL1D, ISPL2D, ISPLPE, PSPPV           *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : GÅnter Palm                                        *
C  Date     : 04.15.1988                                         *
C  Source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XN(0:N), FN(0:N), T(0:N), BX(0:N), CX(0:N),
     +                 DX(0:N), BY(0:N), CY(0:N), DY(0:N),
     +                 ALPHA(2), BETA(2), DUMMY(1:5*N+1)
C
C-----Check input data
C
      IERR = -1
      IF (N .LT. 3) RETURN
      IF (IB .LT. 1  .OR.  IB .GT. 4) THEN
        IERR = -2
        RETURN
      ENDIF
C
C-----Find and check parameter values
C
      IF (MT .GT. 0) THEN
C
C       Find parameter values in SUBROUTINE PSPPV
C
        CALL PSPPV (N,XN,FN,T,MT,IERR)
        IF (IERR .NE. 0) THEN
          IERR = -3
          RETURN
        ENDIF
      ELSE
C
C       Check the given parameter values
C
        IERR = -3
        DO 20 I=0,N-1,1
          IF (T(I+1) .LE. T(I)) RETURN
   20   CONTINUE
      ENDIF
C
C-----Compute the spline coefficients CX(I) and CY(I). Note that
C     on second or subsequent calls, the system matrix need not
C     be decomposed anew.
C
      IF (IB .EQ. 1) THEN
        CALL ISPL1D (N,T,XN,ALPHA(1),BETA(1),1,BX,CX,DX,DUMMY(1),
     +               DUMMY(N+1),DUMMY(2*N),DUMMY(3*N-1),IERR)
        IF (IERR .NE. 0) RETURN
        CALL ISPL1D (N,T,FN,ALPHA(2),BETA(2),2,BY,CY,DY,DUMMY(1),
     +               DUMMY(N+1),DUMMY(2*N),DUMMY(3*N-1),IERR)
      ELSEIF (IB .EQ. 2) THEN
        CALL ISPL2D (N,T,XN,ALPHA(1),BETA(1),1,BX,CX,DX,DUMMY(1),
     +               DUMMY(N+1),DUMMY(2*N),DUMMY(3*N-1),IERR)
        IF (IERR .NE. 0) RETURN
        CALL ISPL2D (N,T,FN,ALPHA(2),BETA(2),2,BY,CY,DY,DUMMY(1),
     +               DUMMY(N+1),DUMMY(2*N),DUMMY(3*N-1),IERR)
      ELSEIF (IB .EQ. 3) THEN
C
C       Compute the tangent vectors using the first
C       derivatives of SX and SY
C
        UB = 1.0D10
        IF (DABS(ALPHA(1)) .GE. UB) THEN
          ALPHAX = 0.0D0
          ALPHAY = DSIGN(1.0D0,FN(1)-FN(0))
        ELSE
          ROOT   = DSQRT(1.0D0/(1.0D0 + ALPHA(1)*ALPHA(1)))
          ALPHAX = DSIGN(ROOT,XN(1)-XN(0))
          ALPHAY = ALPHAX*ALPHA(1)
        ENDIF
        IF (DABS(BETA(1)) .GE. UB) THEN
          BETAX = 0.0D0
          BETAY = DSIGN(1.0D0,FN(N)-FN(N-1))
        ELSE
          ROOT  = DSQRT(1.0D0/(1.0D0 + BETA(1)*BETA(1)))
          BETAX = DSIGN(ROOT,XN(N)-XN(N-1))
          BETAY = BETAX*BETA(1)
        ENDIF
        CALL ISPL1D (N,T,XN,ALPHAX,BETAX,1,BX,CX,DX,DUMMY(1),
     +               DUMMY(N+1),DUMMY(2*N),DUMMY(3*N-1),IERR)
        IF (IERR .NE. 0) RETURN
        CALL ISPL1D (N,T,FN,ALPHAY,BETAY,2,BY,CY,DY,DUMMY(1),
     +               DUMMY(N+1),DUMMY(2*N),DUMMY(3*N-1),IERR)
      ELSE
        CALL ISPLPE (N,T,XN,1,BX,CX,DX,DUMMY(1),DUMMY(N+2),
     +               DUMMY(2*N+2),DUMMY(3*N+2),DUMMY(4*N+2),IERR)
        IF (IERR .NE. 0) RETURN
        CALL ISPLPE (N,T,FN,2,BY,CY,DY,DUMMY(1),DUMMY(N+2),
     +               DUMMY(2*N+2),DUMMY(3*N+2),DUMMY(4*N+2),IERR)
      ENDIF
      RETURN
      END
