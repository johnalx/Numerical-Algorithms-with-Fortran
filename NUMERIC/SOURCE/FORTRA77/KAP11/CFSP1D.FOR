C[BA*)
C[LE*)
      SUBROUTINE CFSP1D (N,XN,FN,W,ALPHA,BETA,MREP,A,B,C,D,
C[IX{CFSP1D}*)
     +                   H,H1,H2,DM,DU1,DU2,RS,IERR)
C
C*****************************************************************
C                                                                *
C  CFSP1D computes the coefficients A(I), B(I), C(I), D(I),      *
C  I=0, 1. ..., N-1, of a cubic fitting spline with prescribed   *
C  first end point derivative.                                   *
C[BE*)
C  The spline is represented in the form:                        *
C                                                                *
C  S(X) = A(I) + B(I)(X-XN(I)) + C(I)(X-XN(I))**2 +              *
C                              + D(I)(X-XN(I))**3                *
C                                                                *
C  for X in the interval [XN(I),XN(I+1)], I=0, 1, ..., N-1.      *
C                                                                *
C                                                                *
C  ASSUMPTIONS:    1.         N > 4                              *
C  ============    2.     XN(I) < XN(I+1), I=0, 1, ..., N-1      *
C                  3.      W(I) > 0.0    , I=0, 1, ..., N        *
C                                                                *
C                                                                *
C  REMARK:  CFSP1D should not be called directly, but rather via *
C  =======  the subroutine CFSPNP, or in case of parametric      *
C           splines via the subroutine CFSPPA. The subroutines   *
C           CFSPNP and CFSPPA also check the above assumptions.  *
C                                                                *
C                                                                *
C  INPUT PARAMETERS:                                             *
C  =================                                             *
C  N  :  Index of the last node                                  *
C  XN :  vector XN(0:N); XN(I) is the Ith node, I = 0, ..., N    *
C  FN :  vector FN(0:N); FN(I) is the data at the node XN(I)     *
C  W  :  vector W(0:N);  W(I) is the weight of FN(I)             *
C                                                                *
C  ALPHA :  first end point derivative at XN(0)                  *
C  BETA  :  first end point derivative at XN(N)                  *
C                                                                *
C  MREP  :  indicator used for repeated call of the subroutine:  *
C           MREP = 1: The system matrix elements must be computed*
C                     This matrix must be factored via subroutine*
C                     FDISY in order to find C(I).               *
C           MREP = 2: The right-hand side only needs to be com-  *
C                     puted. We can use the vectors DM, DU1 and  *
C                     DU2 computed during the first pass of sub- *
C                     routine FDISYS to find the solution.       *
C                     This avoids a repeat factorization in case *
C                     of parametric splines.                     *
C                     The elements in H, H1, H2, DM, DU1 and DU2 *
C                     must not be altered after the first call.  *
C  REMARK: For parametric splines with differing weights WX(I)   *
C          not equal to WY(I) for at least one index I, one      *
C          must work with MREP = 1.                              *
C                                                                *
C                                                                *
C  AUXILIARY VARIABLES:                                          *
C  ====================                                          *
C  H   :]                                                        *
C  H1  :]  N-vectors H(0:N-1), H1(0:N-1), H2(0:N-1)              *
C  H2  :]                                                        *
C                                                                *
C  DM  :]                                                        *
C  DU1 :]  (N-1)-vectors dimensioned as  ..(1:N-1)               *
C  DU2 :]                                                        *
C  RS  :]                                                        *
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
C          = -1 :  N < 5                                         *
C          =  1 :  FDISY did not run correctly (matrix singular) *
C          = -6 :  Wrong input for MREP                          *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  required subroutines: FDISY, FDISYS                           *
C                                                                *
C                                                                *
C  Reference: Engeln-M�llges, G.; Reutter, F., [ENGE87].         *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : G�nter Palm                                        *
C  Date     : 04.18.1988                                         *
C  Source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XN(0:N), FN(0:N), W(0:N), A(0:N), B(0:N),
     +                 C(0:N), D(0:N), H(0:N-1), H1(0:N-1),
     +                 H2(0:N-1), DM(1:N-1), DU1(1:N-1),
     +                 DU2(1:N-1), RS(1:N-1), SF(10)
      SAVE SF
C
C-----Check MREP for repeated calls
C
      IERR = -6
      IF (MREP .NE. 1  .AND.  MREP .NE. 2) RETURN
C
C-----Compute auxiliary values and system matrix elements
C     i.e., its main and two co-diagonals, if in the first pass
C
      IF (MREP .EQ. 1) THEN
C
C       Auxiliary variables
C
        DO 10 I=0,N-1,1
          H(I)  = XN(I+1) - XN(I)
          H1(I) = 1.0D0/H(I)
          C(I)  = H1(I)*H1(I)
          B(I)  = 6.0D0/W(I)
   10   CONTINUE
        B(N) = 6.0D0/W(N)
        DO 20 I=0,N-2,1
          H2(I) = H1(I) + H1(I+1)
   20   CONTINUE
C
C       second co-diagonal
C
        DO 30 I=1,N-3,1
          DU2(I) = B(I+1)*H1(I)*H1(I+1)
   30   CONTINUE
C
C       first co-diagonal
C
        DO 40 I=1,N-2,1
          DU1(I) = H(I) - B(I)*H1(I)*H2(I-1) -B(I+1)*H1(I)*H2(I)
   40   CONTINUE
C
C       main diagonal
C
        DO 50 I=1,N-1,1
          K = I-1
          DM(I) = 2.0D0*(H(K)+H(I)) + B(K)*C(K) +
     +            B(I)*H2(K)*H2(K) +  B(I+1)*C(I)
   50   CONTINUE
C
C       Adjust the outside corner elements of the matrix using the
C       auxiliary variables. Their values will be needed in the
C       second pass for parametric splines, hence we will store them
C       in auxilary varibles secured via SAVE commands.
C
        SF(1)  = H(0) - B(0)*C(0) - B(1)*H2(0)*H1(0)
        SF(2)  = B(1)*H1(0)*H1(1)
        SF(3)  = B(N-1)*H1(N-2)*H1(N-1)
        SF(4)  = H(N-1) - B(N-1)*H2(N-2)*H1(N-1) - B(N)*C(N-1)
        SF(5)  = H1(0) * (B(1)+B(0)) + 2.0D0*H(0)*H(0)
        SF(6)  = H1(N-1) * (B(N)+B(N-1)) + 2.0D0*H(N-1)*H(N-1)
        SF(7)  = B(1)*H2(0) + B(0)*H1(0) - H(0)*H(0)
        SF(8)  = B(N-1)*H2(N-2) + B(N)*H1(N-1) - H(N-1)*H(N-1)
        SF(9)  = B(1)*H1(1)
        SF(10) = B(N-1)*H1(N-2)
C
        DM(1)    = DM(1) + (SF(1)/SF(5)*SF(7))
        DU1(1)   = DU1(1) - (SF(1)/SF(5)*SF(9))
        DM(2)    = DM(2) - (SF(2)/SF(5)*SF(9))
        DM(N-2)  = DM(N-2) - (SF(3)/SF(6)*SF(10))
        DU1(N-2) = DU1(N-2) + (SF(3)/SF(6)*SF(8))
        DM(N-1)  = DM(N-1) + (SF(4)/SF(6)*SF(8))
      ENDIF
C
C-----Computing the right-hand side
C
      DUMMY1 = (FN(1) - FN(0)) * H1(0)
      DO 60 I=1,N-1,1
        DUMMY2 = (FN(I+1) - FN(I)) * H1(I)
        RS(I)  = 3.0D0*(DUMMY2 - DUMMY1)
        DUMMY1 = DUMMY2
   60 CONTINUE
C
      DUMMY3  = 3.0D0*((FN(1)-FN(0)) - ALPHA*H(0))
      DUMMY4  = 3.0D0*(BETA*H(N-1) - (FN(N)-FN(N-1)))
C
      RS(1)   = RS(1) - (SF(1)/SF(5)*DUMMY3)
      RS(2)   = RS(2) - (SF(2)/SF(5)*DUMMY3)
      RS(N-2) = RS(N-2) - (SF(3)/SF(6)*DUMMY4)
      RS(N-1) = RS(N-1) - (SF(4)/SF(6)*DUMMY4)
C
C-----Compute the coefficients C(1) to C(N-1) from
C     the system of equations
C
      IF (MREP .EQ. 1) THEN
C
C       In case we must decompose the system matrix
C
        CALL FDISY (N-1,DM,DU1,DU2,RS,C(1),IFLAG)
        IF (IFLAG .NE. 1) THEN
          IF (IFLAG .EQ. -2) THEN
            IERR = -1
          ELSE
            IERR =  1
          ENDIF
          RETURN
        ENDIF
      ELSE
C
C       When no factorization is necessary
C
        CALL FDISYS (N-1,DM,DU1,DU2,RS,C(1))
      ENDIF
      IERR = 0
C
C-----Compute the remaining spline coefficients
C
      C(0) = (DUMMY3 + C(1)*SF(7) - C(2)*SF(9))/SF(5)
      C(N) = (DUMMY4 + C(N-1)*SF(8) - C(N-2)*SF(10))/SF(6)
C
      A(0) = FN(0) + 2.0D0/W(0)*H1(0)*(C(0)-C(1))
      DO 70 I=1,N-1,1
        A(I) = FN(I) - 2.0D0/W(I)*(C(I-1)*H1(I-1) -
     +         H2(I-1)*C(I) + C(I+1)*H1(I))
   70 CONTINUE
      A(N) = FN(N) - 2.0D0/W(N)*H1(N-1)*(C(N-1)-C(N))
C
      DO 80 I=0,N-1,1
        B(I) = H1(I)*(A(I+1)-A(I)) - H(I)/3.0D0*(C(I+1)+2.0D0*C(I))
        D(I) = H1(I)/3.0D0*(C(I+1)-C(I))
   80 CONTINUE
      RETURN
      END
