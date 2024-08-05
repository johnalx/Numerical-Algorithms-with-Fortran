C[BA*)
C[KA{F 11}
C[  {Polynomial Fitting Splines}
C[  {Polynomial Cubic Fitting Splines for Constructing
C[   Smooth Curves}*)
C[FE{F 11.3}
C[  {Non-Parametric Cubic Fitting Splines}
C[  {Non-Parametric Cubic Fitting Splines}*)
C[LE*)
      SUBROUTINE CFSPNP (N,XN,FN,W,IB,ALPHA,BETA,A,B,C,D,AUXF,
C[IX{CFSPNP}*)
     +                   IERR)
C
C*****************************************************************
C                                                                *
C  CFSPNP computes the coefficients A(I), B(I), C(I), D(I) for   *
C  I=0, 1, ..., N-1 of a nonparametric cubic fitting spline.     *
C  The end point condition is to be prescribed via the parameter *
C  IB.                                                           *
C[BE*)
C  The splinefunction is represented in the form:                *
C                                                                *
C  S(X) = A(I) + B(I)(X-XN(I)) + C(I)(X-XN(I))**2 +              *
C                              + D(I)(X-XN(I))**3                *
C                                                                *
C  for X in the interval [XN(I),XN(I+1)], I=0, 1, ..., N-1.      *
C                                                                *
C                                                                *
C  ASSUMPTIONS:    1.         N > 4      , for IB = 1, 2 or 3    *
C  ============               N > 5      , for IB =  4           *
C                  2.     XN(I) < XN(I+1), I=0, 1, ..., N-1      *
C                  3.      W(I) > 0.0    , I=0, 1, ..., N        *
C                  4.      W(0) = W(N)   , for IB = 4            *
C                  5.     FN(0) = FN(N)  , for IB = 4            *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  N  :  Index of the last node                                  *
C  XN :  vector XN(0:N); XN(I) is the Ith node, I = 0, ..., N    *
C  FN :  vector FN(0:N); FN(I) is the data at the node XN(I)     *
C  W  :  vector W(0:N);  W(I) is the weight of FN(I)             *
C                                                                *
C  IB :  determines end point condition:                         *
C        IB = 1:  first end point derivative prescribed          *
C        IB = 2:  second end point derivative prescribed         *
C        IB = 3:  third end point derivative prescribed          *
C        IB = 4:  periodic spline                                *
C                                                                *
C  ALPHA :  IB end point derivative at XN(0) ] for IB=1,2,3;     *
C  BETA  :  IB end point derivative at XN(N) ] meaningless for   *
C                                              IB=4              *
C                                                                *
C           (A natural fitting spline will be achieved for IB=2  *
C            and ALPHA = BETA =0.0)                              *
C                                                                *
C                                                                *
C  AUXILIARY VARIABLES:                                          *
C  ====================                                          *
C  AUXF :  vector AUXF(1:14*N-10)                                *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  ==================                                            *
C  A    :  Vector A(0:N) ]  The entries in positions 0 to N-1    *
C  B    :  Vector B(0:N) ]  contain the spline coefficients for  *
C  C    :  Vector C(0:N) ]  S. The entries in A(N), B(N), C(N)   *
C  D    :  Vector D(0:N) ]  and D(N) are auxiliary variables.    *
C                                                                *
C  IERR :  error parameter                                       *
C          =  0 :  All is o.k.                                   *
C          = -1 :  N < 5  if IB = 1, 2 or 3                      *
C                  N < 6  if IB = 4                              *
C          = -2 :  IB < 1  or  IB > 4                            *
C          = -3 :  Inadmissable weight W                         *
C          = -4 :  nodes XN(I) not ordered monotonically:        *
C                  XN(I) >= XN(I+1) for some  I=0, 1, ..., N-1   *
C          = -5 :  IB = 4 and FN(0) not equal to FN(N) or        *
C                  W(0) not equal to W(N)                        *
C          =  1 :  Error in FDISY, FDIAG or NCYFSY (numerically  *
C                  singular system matrix)                       *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  required subroutines: CFSP1D, CFSP2D, CFSP3D, CFSPPE,         *
C                        FDISY, FDISYS, NCYFSY, NCYFSP,          *
C                        NCYFSS, FDIAG                           *
C                                                                *
C                                                                *
C  Reference: Engeln-MÅllges, G.; Reutter, F., [ENGE87].         *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : GÅnter Palm                                        *
C  Date     : 04.18.1988                                         *
C  Source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XN(0:N), FN(0:N), W(0:N), A(0:N), B(0:N),
     +                 C(0:N), D(0:N), AUXF(1:14*N-10)
C
C-----Checking the assumptions
C
      IERR = -1
      IF (N .LT. 5) RETURN
      DO 10 I=0,N-1,1
        IF (XN(I) .GE. XN(I+1)) THEN
          IERR = -4
          RETURN
        ENDIF
   10 CONTINUE
      DO 20 I=0,N,1
        IF (W(I) .LE. 0.0D0) THEN
          IERR = -3
          RETURN
        ENDIF
   20 CONTINUE
C
C-----Compute the spline coefficients
C
      IF (IB .EQ. 1) THEN
        CALL CFSP1D (N,XN,FN,W,ALPHA,BETA,1,A,B,C,D,
     +               AUXF(1),AUXF(N+1),AUXF(2*N+1),AUXF(3*N+1),
     +               AUXF(4*N),AUXF(5*N-1),AUXF(6*N-2),IERR)
      ELSEIF (IB .EQ. 2) THEN
        CALL CFSP2D (N,XN,FN,W,ALPHA,BETA,1,A,B,C,D,
     +               AUXF(1),AUXF(N+1),AUXF(2*N+1),AUXF(3*N+1),
     +               AUXF(4*N),AUXF(5*N-1),AUXF(6*N-2),IERR)
      ELSEIF (IB .EQ. 3) THEN
        CALL CFSP3D (N,XN,FN,W,ALPHA,BETA,A,B,C,D,
     +               AUXF(1),AUXF(N),AUXF(2*N-1),AUXF(3*N-2),
     +               AUXF(4*N-3),AUXF(5*N-4),AUXF(6*N-4),IERR)
      ELSEIF (IB .EQ. 4) THEN
        IF (N .LT. 6) RETURN
        CALL CFSPPE (N,XN,FN,W,1,A,B,C,D,AUXF(1),AUXF(N+2),
     +               AUXF(2*N+3),AUXF(3*N+4),AUXF(4*N+5),
     +               AUXF(5*N+5),IERR)
      ELSE
        IERR = -2
      ENDIF
      RETURN
      END
