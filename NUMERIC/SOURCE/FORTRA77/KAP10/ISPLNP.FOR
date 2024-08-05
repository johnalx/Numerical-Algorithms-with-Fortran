C[BA*)
C[KA{F 10}
C[  {Interpolating Polynomial Splines}
C[  {Interpolating Polynomial Splines for Constructing Smooth
C[   Curves}*)
C[FE{F 10.1.2}
C[  {Non-Parametric Cubic Splines}
C[  {Non-Parametric Cubic Splines}*)
C[LE*)
      SUBROUTINE ISPLNP (N,XN,FN,IB,ALPHA,BETA,B,C,D,DUMMY,IERR)
C[IX{ISPLNP}*)
C
C*********************************************************************
C                                                                    *
C  ISPLNP computes the coefficients B(I), C(I), D(I) for I=0,1,...,  *
C  N-1 of a nonparametric cubic interpolating spline for various end *
C  point conditions.                                                 *
C[BE*)
C  The end point condition can be specified using the parameter IB   *
C  The spline function S has the form:                               *
C                                                                    *
C  S(X) = FN(I) + B(I)(X-XN(I)) + C(I)(X-XN(I))**2 +                 *
C                               + D(I)(X-XN(I))**3                   *
C                                                                    *
C  for X in the interval [XN(I),XN(I+1)], I=0,1,...,N-1.             *
C                                                                    *
C                                                                    *
C  ASSUMPTIONS:    1.         N > 2                                  *
C  ============    2.     XN(I) < XN(I+1), I=0,1,...,N-1             *
C                  3.     FN(0) = FN(N)  , if IB = 4                 *
C                                                                    *
C                                                                    *
C  INPUT PARAMETERS:                                                 *
C  =================                                                 *
C  N  :  index of the final node                                     *
C  XN :  vector XN(0:N); the nodes XN(I), I = 0,1,...,N              *
C  FN :  vector FN(0:N); the functional values FN(I) = FN( XN(I) )   *
C                                                                    *
C  IB :  specifies the end point condition:                          *
C        IB = 1:  first end point derivative                         *
C        IB = 2:  second end point derivative                        *
C        IB = 3:  third end point derivative                         *
C        IB = 4:  periodic spline                                    *
C        IB = 5:  'not-a-node' - condition                           *
C                                                                    *
C  ALPHA :  IBth derivative at XN(0) ]  used only for IB = 1,2,3;    *
C  BETA  :  IBth derivative at XN(N) ]  not used for IB = 4,5.       *
C                                                                    *
C  (A natural spline will result by setting ALPHA = BETA = 0.0 and   *
C   IB = 2.)                                                         *
C                                                                    *
C                                                                    *
C  AUXILIARY VARIABLES:                                              *
C  ====================                                              *
C  DUMMY :  vector DUMMY(1:5*N+1)                                    *
C                                                                    *
C                                                                    *
C  OUTPUT PARAMETERS:                                                *
C  ==================                                                *
C  FN :  ]  N+1-vectors ..(0:N);                                     *
C  B  :  ]  The first N entries of B, C and D are the spline         *
C  C  :  ]  coefficients for S. B(N), C(N), D(N) are auxiliary       *
C  D  :  ]  variables.                                               *
C  IERR :  error parameter                                           *
C          =  0 :  All is ok                                         *
C          = -1 :  N < 4                                             *
C          = -2 :  IB < 1 or IB > 5                                  *
C          = -3 :  nodes not monotone, XN(I) => XN(I+1) for some     *
C                  I = 0, 1, ..., N-1                                *
C          = -4 :  IB = 4 and FN(0) not equal to FN(N)               *
C          =  1 :  crash in SUBROUTINE TRDSY, TRDIG or CYTSY,        *
C                  system matrix is numerically singular             *
C                                                                    *
C--------------------------------------------------------------------*
C                                                                    *
C  Subroutines required: ISPL1D, ISPL2D, ISPL3D, ISPLPE, ISPLNK,     *
C                        TRDSY, TRDSYS, CYTSY, CYTSYS, TRDIG         *
C                                                                    *
C*********************************************************************
C                                                                    *
C  Author   : GÅnter Palm                                            *
C  Date     : 15.04.1988                                             *
C  Source   : FORTRAN 77                                             *
C                                                                    *
C[BA*)
C*********************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION XN(0:N), FN(0:N), B(0:N), C(0:N), D(0:N),
     +                 DUMMY(1:5*N+1)
C
C-----Check input data
C
      IERR = -1
      IF (N .LT. 3) RETURN
      DO 10 I=0,N-1,1
        IF (XN(I) .GE. XN(I+1)) THEN
          IERR = -3
          RETURN
        ENDIF
   10 CONTINUE
C
C-----Compute spline coefficients for the indicated
C     end point condition
C
      IF (IB .EQ. 1) THEN
        CALL ISPL1D (N,XN,FN,ALPHA,BETA,1,B,C,D,DUMMY(1),DUMMY(N+1),
     +               DUMMY(2*N),DUMMY(3*N-1),IERR)
      ELSEIF (IB .EQ. 2) THEN
        CALL ISPL2D (N,XN,FN,ALPHA,BETA,1,B,C,D,DUMMY(1),DUMMY(N+1),
     +               DUMMY(2*N),DUMMY(3*N-1),IERR)
      ELSEIF (IB .EQ. 3) THEN
        CALL ISPL3D (N,XN,FN,ALPHA,BETA,B,C,D,DUMMY(1),DUMMY(N+1),
     +               IERR)
      ELSEIF (IB .EQ. 4) THEN
        CALL ISPLPE (N,XN,FN,1,B,C,D,DUMMY(1),DUMMY(N+2),
     +               DUMMY(2*N+2),DUMMY(3*N+2),DUMMY(4*N+2),IERR)
      ELSEIF (IB .EQ. 5) THEN
        CALL ISPLNK (N,XN,FN,B,C,D,DUMMY(1),DUMMY(N+1),DUMMY(2*N),
     +               IERR)
      ELSE
        IERR = -2
      ENDIF
      RETURN
      END
