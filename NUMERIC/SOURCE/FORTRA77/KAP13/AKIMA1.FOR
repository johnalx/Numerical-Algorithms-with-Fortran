C[BA*)
C[KA{F 13}
C[  {Akima and Renner Subsplines}
C[  {Akima and Renner Subsplines}*)
C[FE{F 13.1}{Akima Subsplines}{Akima Subsplines}*)
C[LE*)
      SUBROUTINE AKIMA1(N,XN,FN,NK,BETA,B,C,D,ISWTCH,HELP,IERR)
C[IX{AKIMA1}*)
C
C*****************************************************************
C                                                                *
C  The program AKIMA1 computes the coefficients  B(I), C(I) and  *
C  D(I) for I=0, ... , N-1 of an interpolating cubic AKIMA       *
C  spline, which can be either periodic or nonperiodic.          *
C[BE*)
C  The subspline has the representation                          *
C                                                                *
C  S(X)=FN(I)+B(I)(X-XN(I))+C(I)(X-XN(I))**2+D(I)(X-XN(I))**3    *
C                                                                *
C  for any point X in the subinterval [XN(I), XN(I+1)] for       *
C  I=0,..., N-1.                                                 *
C                                                                *
C                                                                *
C  ASSUMPTIONS:                                                  *
C  ============                                                  *
C                1. N >= 4 or NK >= 4                            *
C                2. If 0.0 < BETA < 1.0, we must have            *
C                   NK >= N + INT((N+1)/2), otherwise NK = N     *
C                3. The nodes XN(I), I=0, ..., N, must be        *
C                   ordered monotonically, i.e., XN(I) < XN(I+1) *
C                   for I=0, ... , N-1                           *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  N       : Index of the last node                              *
C  XN      : DOUBLE PRECISION (NK+1)-vector XN(0:NK), containing *
C            the nodes XN(I) for I=0, ... , N                    *
C  FN      : DOUBLE PRECISION (NK+1)-vector FN(0:NK), containing *
C            the functional values FN(I) at XN(I) for I=0,...,N  *
C  NK      : NK = N+INT((N+1)/2) maximal number of nodes allowed *
C            when using rounded corners. If corners are not      *
C            rounded, we use: NK = N                             *
C  BETA    : If 0.0 < BETA < 1.0 the corners are rounded, other- *
C            wise corners are kept                               *
C            In the periodic case, we do not round a corner that *
C            may exist at XN(0) = XN(N), even if 0.0 < BETA < 1.0*
C  ISWTCH  : = 0, nonperiodic spline                             *
C            = 1, periodic spline                                *
C            In the periodic case, the interval [XN(0), XN(N)]   *
C            must be an interval of periodicity with FN(0) =     *
C            FN(N).                                              *
C                                                                *
C                                                                *
C  AUXILIARY VARIABLES:                                          *
C  ====================                                          *
C  HELP    : DOUBLE PRECISION array HELP(-2:NK+1,1:4)            *
C                                                                *
C                                                                *
C  OUTPUT PARAMETERS:                                            *
C  =================                                             *
C  N       : Number of the last node. If corners are rounded,    *
C            i.e., if 0.0 < BETA < 1.0 , then the output value   *
C            of N can differ from its input. When corners are    *
C            rounded, the set of nodes can maximally be enlarged *
C            by INT((N+1)/2) new nodes.                          *
C  XN      : DOUBLE PRECISION vector XN(0:NK) containing the     *
C            nodes XN(I), I=0,...,N. If 0.0 < BETA < 1.0, then   *
C            the output nodes can differ from the input nodes.   *
C  FN      : DOUBLE PRECISION vector FN(0:NK), containing the    *
C            functional values FN(I) at XN(I) for I = 0, ..., N. *
C            If 0.0 < BETA < 1.0, the node and functional values *
C            XN(I) and FN(I) can differ from their input values. *
C  B       : DOUBLE PRECISION vector B(0:NK-1) ] B, C and D con- *
C  C       : DOUBLE PRECISION vector C(0:NK-1) ] tain the coeffi-*
C  D       : DOUBLE PRECISION vector D(0:NK-1) ] cients of the   *
C                                     subspline for I=0 to NK-1. *
C                                                                *
C  IERR    : error parameter                                     *
C            = 0,  all is ok                                     *
C            =-1, N < 4 or NK < 4                                *
C            =-2, the XN(I) are not monotonically ordered        *
C            =-3, NK < N+INT((N+1)/2) while 0.0 < BETA < 1.0     *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Required subroutine: MACHPD                                   *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : Gisela Engeln-Mllges                              *
C  Date     : 04.09.1993                                         *
C  Source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C     Declarations
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XN(0:NK), FN(0:NK), B(0:NK-1), C(0:NK-1),
     +                 D(0:NK-1), HELP(-2:NK+1, 1:4)
C
C     Check input parameters
C
      IERR = 0
      IF (N .LT. 4 .OR. NK .LT. 4) THEN
         IERR = -1
         RETURN
      ENDIF
      IF (      BETA .GT. 0.0D0 .AND. BETA .LT. 1.0D0
     F    .AND. NK   .LT. N+INT(0.5D0*(N+1))         ) THEN
            IERR = -3
            RETURN
      ENDIF
      IF ( ISWTCH .NE. 0 ) THEN
         ISWTCH = 1
      ENDIF
C
C     Compute the machine constant
C
      FMACHP = 1.0D0
    5 FMACHP = 0.5D0 * FMACHP
      IF (MACHPD(1.0D0+FMACHP) .EQ. 1) GOTO 5
      FMACHP = 2.0D0 * FMACHP
      EPS = 4.0D0 * FMACHP
C
C     Compute the lengths of the subintervals and store in the first
C     column of the auxiliary array HELP
C
      DO 10, I = 0, N-1
         HELP(I,1) = XN(I+1) - XN(I)
         IF ( HELP(I,1) .LE. EPS ) THEN
            IERR = -2
            RETURN
         ENDIF
   10 CONTINUE
C
C     Compute the secant slopes and store in the second column
C     of the auxiliary array HELP
C
      DO 20, I = 0, N-1
         HELP(I,2) = (FN(I+1)-FN(I))/HELP(I,1)
   20 CONTINUE
C
C     Compute the magnitude of the slope differences and store
C     in the third column of HELP
C
      DO 30, I = 0, N-2
         HELP(I,3) = DABS(HELP(I+1,2)-HELP(I,2))
   30 CONTINUE

      IF ( DABS(FN(0)-FN(N)) .GT. EPS ) THEN
         ISWTCH = 0
      ENDIF

      IF ( BETA .GT. 0.0D0 .AND. BETA .LT. 1.0D0 ) THEN
         IF ( ISWTCH .EQ. 1 ) THEN
C
C     Prepare loop
C
            I = 1
            IMAX = N-1
            HELP(-1,3)  = DABS(HELP(0,2)-HELP(N-1,2))
            HELP(N-1,3) = HELP(-1,3)
         ELSE
            I = 2
            IMAX = N-2
         ENDIF

   35    XL = HELP(I-2,3) + HELP(I,3)
         XR = HELP(I-1,3)
C
C     Eliminate existing corners
C
         IF ( XL.LE.EPS .AND. XR.GT.EPS ) THEN
C
C     Relabel points I to N
C
            DO 40, J = N,I, -1
               XN(J+1) = XN(J)
               FN(J+1) = FN(J)
   40       CONTINUE
C
C     Reassign interval lengths and slopes from I to N-1
C
            DO 50, J = N-1, I, -1
               HELP(J+1,1) = HELP(J,1)
               HELP(J+1,2) = HELP(J,2)
   50       CONTINUE
C
C     Reassign slope differences from I to IMAX
C
            DO 60, J = IMAX, I, -1
               HELP(J+1,3) = HELP(J,3)
   60       CONTINUE
C
C     Generate new points labelled I and I+1
C
            XL = HELP(I-1,1)
            XR = HELP(I+1,1)
            XB = BETA*DMIN1(XL,XR)
            XLAMDA = XB/XL
            XMUE   = XB/XR
            XN(I) = XN(I) - XLAMDA*HELP(I-1,1)
            FN(I) = FN(I) - XLAMDA*(FN(I)-FN(I-1))
            XN(I+1) = XN(I+1) + XMUE*HELP(I+1,1)
            FN(I+1) = FN(I+1) + XMUE*(FN(I+2)-FN(I+1))
C
C     Compute new interval lengths
C
            DO 70, J = I-1, I+1
               HELP(J,1) = XN(J+1) - XN(J)
  70        CONTINUE
C
C     Compute new slopes and slope differences
C
            HELP(I,2) = (FN(I+1)-FN(I))/HELP(I,1)
            DO 80, J = I-1, I
               HELP(J,3) = DABS(HELP(J+1,2)-HELP(J,2))
  80        CONTINUE
C
C     Increase number of nodes
C
            N = N+1
            IMAX = IMAX + 1
         ENDIF
C
C     Set index for next point
C
         I = I + 1
         IF ( I .LE. IMAX ) GOTO 35
      ENDIF

      IF ( ISWTCH .EQ. 1 ) THEN
C
C     For the periodic case, form additional slope data
C
         HELP(-2,2)  = HELP(N-2,2)
         HELP(-1,2)  = HELP(N-1,2)
         HELP(N,2)   = HELP(0,2)
         HELP(N+1,2) = HELP(1,2)
      ELSE
C
C     For the nonperiodic case, provide additional slopes as well
C
         HELP(-2,2)  = 3.0D0*HELP(0,2)   - 2.0D0*HELP(1,2)
         HELP(-1,2)  = 2.0D0*HELP(0,2)   - HELP(1,2)
         HELP(N,2)   = 2.0D0*HELP(N-1,2) - HELP(N-2,2)
         HELP(N+1,2) = 3.0D0*HELP(N-1,2) - 2.0D0*HELP(N-2,2)
      ENDIF
C
C     Compute additional slope differences
C
      HELP(-2,3)  = DABS(HELP(-1,2)  - HELP(-2,2))
      HELP(-1,3)  = DABS(HELP(0,2)   - HELP(-1,2))
      HELP(N-1,3) = DABS(HELP(N,2)   - HELP(N-1,2))
      HELP(N,3)   = DABS(HELP(N+1,2) - HELP(N,2))
C
C     Compute the left and right handed slopes at the points 0 and N-1,
C     and store in column 4 of HELP and in B
C
      DO 90, I = 0, N-1
         XL = HELP(I-2,3)
         XR = HELP(I,3)
         IF ( XL+XR .GT. EPS ) THEN
            ALPHA = XL/(XL+XR)
            HELP(I,4) = HELP(I-1,2) + ALPHA*(HELP(I,2)-HELP(I-1,2))
            B(I) = HELP(I,4)
         ELSE
            HELP(I,4) = HELP(I-1,2)
            B(I) = HELP(I,2)
         ENDIF
   90 CONTINUE
C
C     Compute the left handed slope at the point N
C
      XL = HELP(N-2,3)
      XR = HELP(N,3)
      IF (XL+XR .GT. EPS ) THEN
         ALPHA = XL/(XL+XR)
         HELP(N,4) = HELP(N-1,2) + ALPHA*(HELP(N,2)-HELP(N-1,2))
      ELSE
         HELP(N,4) = HELP(N-1,2)
      ENDIF
C
C     Compute the coefficients C(I) and D(I)
C
      DO 100, I = 0, N-1
         H = 1.0D0/HELP(I,1)
         C(I) = (3.0D0*HELP(I,2) - 2.0D0*B(I) - HELP(I+1,4))*H
         D(I) = (B(I) + HELP(I+1,4) - 2.0D0*HELP(I,2))*H*H
  100 CONTINUE

      RETURN
      END
