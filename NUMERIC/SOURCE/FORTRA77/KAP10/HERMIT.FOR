C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 10.2.2}
C[  {Non-Parametric Hermite Splines}
C[  {Non-Parametric Hermite Splines}*)
C[LE*)
      SUBROUTINE HERMIT (N,MARG,X,Y,Y1,BC1,BCN,IREP,A,B,C,D,E,F,
C[IX{HERMIT}*)
     1                 IERR,H,SUPER,AINFER,PRAEC,DEXT,AR1,AR2,AR3)
C
C*****************************************************************
C                                                                *
C     HERMIT computes the coefficients of a hermite polynomial   *
C     spline of fifth degree.                                    *
C[BE*)
C                                                                *
C     NOTE: If several splines are to be computed for the same   *
C           nodes X(I), but for different function values Y(I)   *
C           or derivaltives Y1(I), then the user must set IREP=1 *
C           for the second, third, ... call of HERMIT. In this   *
C           setting the formation and transformation of the      *
C           matrix is not repeated and thus one saves ca 4N      *
C           operations in each case.                             *
C           This occurs naturally when computing parametric      *
C           splines, see SUBROUTINE PARMIT.                      *
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     N       number of nodes X(I)                               *
C     MARG    index for the endpoint condition:                  *
C               MARG = 1 : Periodic spline                       *
C               MARG = 2 : Natural spline                        *
C               MARG = 3 : User specified second derivatives at  *
C                          the end points                        *
C                          In this case the program expects that *
C                           2    2                               *
C                          D Y/DX ( X(1) )  is in BC1  and       *
C                           2    2                               *
C                          D Y/DX ( X(N) )  in BCN.              *
C               MARG = 4 : The user specifies the curvature radii*
C                          R1, RN of the spline at the end points*
C                          In this case the program expects      *
C                          R1 in BC1    and    RN in BCN.        *
C                          Re. concavity : If the radius is      *
C                          positive, the curvature circle is     *
C                          centered above the spline and the     *
C                          curve is concave up. For a negative   *
C                          radius the curve will be concave down.*
C               MARG = 5 : The user specifies the third end point*
C                          derivatives.                          *
C                          In this case the program expects      *
C                           3    3                               *
C                          D Y/DX  ( X(1) )  in BC1  and         *
C                           3    3                               *
C                          D Y/DX  ( X(N) )  in BCN.             *
C     X       vector X(1:N); the nodes X(1),...,X(N),,,,,        *
C     Y       vector Y(1:N); the given functional values         *
C                            Y(1),...,Y(N)                       *
C     Y1      vector Y1(1:N); the given derivatives at the nodes *
C                             DY/DX( X(1) ), ..., DY/DX( X(N) )  *
C     BC1     )  end point conditions for MARG = 3,4 or 5 as     *
C     BCN     )  stipulated by  MARG. left blank for other values*
C             )  of MARG.                                        *
C     IREP    If IREP is not equal to 0, the program assumes that*
C             it is called repeatedly for the same nodes. The    *
C             system matrix is not recomputed in this case.      *
C                                                                *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     A       )                                                  *
C     B       )  N-vectors ..(1:N);                              *
C     C       )  the coefficients of the spline polynomial       *
C     D       )                                                  *
C     E       )                                                  *
C     F       )                                                  *
C     IERR    error code                                         *
C               IERR = 0 : no error                              *
C               IERR > 0 : input error, no spline was computed.  *
C                          In detail:                            *
C               IERR = 1 : MARG < 1  or  MARG > 5                *
C               IERR = 2 : N < 3                                 *
C               IERR = 3 : X(I) not ordered monotonically        *
C               IERR = 4 : Periodic spline specified for non-    *
C                          periodic input:                       *
C                          Y(1).NE.Y(N)  or  Y1(1).NE.Y1(N)      *
C               MARG = 5 : One of the curvature radii is zero    *
C                                                                *
C                                                                *
C     AUXILIARY PARAMETERS:                                      *
C     =====================                                      *
C     H       ) auxiliary vectors ..(1:N);                       *
C     SUPER   ) for repeated calls with IREP = 1 these vectors   *
C     AINFER  ) (except DEXT) must not be altered between calls. *
C     PRAEC   )                                                  *
C     DEXT    )                                                  *
C     AR1     )                                                  *
C     AR2     )                                                  *
C     AR3     )                                                  *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Subroutines required: CYCTR, CYCTRS, TRIDIG, FDEXT            *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : Elmar Pohl                                         *
C  Date     : 09.28.1985                                         *
C  Source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION A(N),B(N),C(N),D(N),E(N),F(N)
      DOUBLE PRECISION SUPER(N),DEXT(N),PRAEC(N),AINFER(N)
      DOUBLE PRECISION X(N),Y(N),Y1(N),H(N),AR1(N),AR2(N),AR3(N)
      DOUBLE PRECISION FDEXT
      IF(IREP .NE. 0) GOTO 4
      IERR=0
      IF(MARG .LT. 1 .OR. MARG .GT. 5) GOTO 501
      IF(N .LT. 3) GOTO 502
      NM1=N-1
      NM2=N-2
      NM3=N-3
      DO 3 I=1,NM1
           HI=X(I+1)-X(I)
           IF(HI .LE. 0.0D0) GOTO 503
           H(I)=HI
    3 CONTINUE
      H(N)=H(1)
      IF(MARG .NE. 4) GOTO2
      IF(BC1 .EQ. 0.0D0 .OR. BCN .EQ. 0.0D0) GOTO 505
    2 ALPHA=1.0D0
      IF(MARG .EQ. 5) ALPHA=8.0D0/9.0D0
      IF(N .NE. 3) GOTO 18
      PRAEC(1)=3.0D0*ALPHA*(1.0D0/H(1)+1.0D0/H(2))
      GOTO 7
   18 REC1=ALPHA/H(1)
      DO 8 I=1,NM3
           REC2=1.0D0/H(I+1)
           PRAEC(I)=3.0D0*(REC1+REC2)
           SUPER(I)=-REC2
           AINFER(I+1)=-REC2
           REC1=REC2
    8 CONTINUE
      PRAEC(NM2)=3.0D0*(1.0D0/H(NM2)+ALPHA/H(NM1))
    7 IF(MARG .NE. 1) GOTO 9
      REC1=1.0D0/H(NM1)
      REC2=1.0D0/H(1)
      PRAEC(NM1)=3.0D0*(REC1+REC2)
      AR2(1)=-REC2
      AINFER(NM1)=-REC1
      AR1(1)=-REC2
      SUPER(NM2)=-REC1
      IF(N .NE. 3) GOTO 9
      SUPER(1)=SUPER(1)-REC2
      AINFER(2)=AINFER(2)-REC2
      GOTO 9
    4 NM1=N-1
      NM2=N-2
    9 CONTINUE
      DO 6 I=1,N
           A(I)=Y(I)
           B(I)=Y1(I)
    6 CONTINUE
      IF(MARG .EQ. 1 .AND. (A(1) .NE. A(N) .OR.
     &                      B(1) .NE. B(N))) GOTO 504
      GOTO (101,101,103,104,105) MARG
  101 BETA1=0.0D0
      BETA2=0.0D0
      GOTO 200
  103 BETA1=0.5D0*BC1/H(1)
      BETA2=0.5D0*BCN/H(NM1)
      GOTO 200
  104 Y21=(1.0D0+B(1)*B(1))**1.5D0/BC1
      Y2N=(1.0D0+B(N)*B(N))**1.5D0/BCN
      BETA1=0.5D0*Y21/H(1)
      BETA2=0.5D0*Y2N/H(NM1)
      GOTO 200
  105 HSQ=H(1)*H(1)
      BETA1=10.0D0*(A(2)-A(1))/3.0D0/H(1)/HSQ-2.0D0*
     1          (2.0D0*B(2)+3.0D0*B(1))/3.0D0/HSQ-BC1/18.0D0
      HSQ=H(NM1)*H(NM1)
      BETA2=-10.0D0*(A(N)-A(NM1))/3.0D0/H(NM1)/HSQ
      BETA2=BETA2+2.0D0*(3.0D0*B(N)+2.0D0*B(NM1))/3.0D0/HSQ+BCN/18.0D0
  200 REC1=1.0D0/H(1)
      A1=A(1)
      A2=A(2)
      B1=B(1)
      B2=B(2)
      DO 10 I=1,NM2
            REC2=1.0D0/H(I+1)
            A3=A(I+2)
            B3=B(I+2)
            DEXT(I)=FDEXT(A1,A2,A3,B1,B2,B3,REC1,REC2)
            A1=A2
            A2=A3
            B1=B2
            B2=B3
            REC1=REC2
   10 CONTINUE
      DEXT(1)=DEXT(1)+BETA1
      DEXT(NM2)=DEXT(NM2)+BETA2
      IF(MARG .NE. 1) GOTO 17
      DEXT(NM1)=FDEXT(A1,A2,A(2),B1,B2,B(2),REC1,1.0D0/H(N))
      IF(N .NE. 3) GOTO 11
      CALL TRIDIG(2,DEXT,AINFER,SUPER,PRAEC,IREP
     1             ,C,ISING,AR1,AR2,AR3)
      GOTO 19
   11 IF(IREP .NE. 0) GOTO 22
      CALL CYCTR(NM1,AINFER,PRAEC,SUPER,AR1,AR2,DEXT,C,IERR1)
      GOTO 19
   22 CALL CYCTRS(NM1,AINFER,PRAEC,SUPER,AR1,AR2,DEXT,C)
   19 DO 14 I=1,NM1
            NMI=N-I
            C(NMI+1)=C(NMI)
   14 CONTINUE
      C(1)=C(N)
      GOTO 13
   17 CALL TRIDIG(NM2,DEXT,AINFER,SUPER,PRAEC,IREP
     1             ,C,ISING,AR1,AR2,AR3)
      DO 12 I=1,NM2
            NMI=N-I
            C(NMI)=C(NMI-1)
   12 CONTINUE
      GOTO (11,202,203,204,205) MARG
  202 C(1)=0.0D0
      C(N)=0.0D0
      GOTO 13
  203 C(1)=0.5D0*BC1
      C(N)=0.5D0*BCN
      GOTO 13
  204 C(1)=0.5D0*Y21
      C(N)=0.5D0*Y2N
      GOTO 13
  205 CONTINUE
      C(1)=H(1)*BETA1+C(2)/3.0D0
      C(N)=H(NM1)*BETA2+C(NM1)/3.0D0
   13 DO 15 I=1,NM1
            D(I)=10.0D0*(A(I+1)-A(I))/H(I)-2.0D0*(2.0D0*B(I+1)+
     &                                            3.0D0*B(I))
            D(I)=(D(I)/H(I)+C(I+1)-3.0D0*C(I))/H(I)
   15 CONTINUE
      D(N)=D(NM1)-(2.0D0*(B(N)-B(NM1))/H(NM1)-2.0D0*(C(N)+
     &                                               C(NM1)))/H(NM1)
      DO 16 I=1,NM1
            HI=H(I)
            E(I)=(0.5D0*(B(I+1)-B(I))/HI-C(I))/HI
            E(I)=(E(I)-0.25D0*(D(I+1)+5.0D0*D(I)))/HI
            F(I)=(((C(I+1)-C(I))/HI-3.0D0*D(I))/HI-
     &                              6.0D0*E(I))/HI/10.0D0
   16 CONTINUE
      RETURN
  501 IERR=1
      RETURN
  502 IERR=2
      RETURN
  503 IERR=3
      RETURN
  504 IERR=4
      RETURN
  505 IERR=5
      RETURN
      END
C
C
C[BA*)
C[LE*)
      DOUBLE PRECISION FUNCTION FDEXT (A1,A2,A3,B1,B2,B3,REC1,REC2)
C[IX{FDEXT}*)
C
C*****************************************************************
C                                                                *
C     Function routine needed for SUBROUTINE HERMIT.             *
C     FDEXT computes one element of the right hand side of the   *
C     system of equations.                                       *
C                                                                *
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      REC1SQ=REC1*REC1
      REC2SQ=REC2*REC2
      FDEXT=10.0D0*((A3-A2)*REC2*REC2SQ-(A2-A1)*REC1*REC1SQ)
      FDEXT=FDEXT+4.0D0*(B1*REC1SQ-1.5D0*(REC2SQ-REC1SQ)*B2-
     &                   B3*REC2SQ)
      RETURN
      END
C
C
C[BA*)
C[LE*)
      DOUBLE PRECISION FUNCTION HMTVAL (N,X0,A,B,C,D,E,F,X,OUTP)
C[IX{HMTVAL}*)
C
C*****************************************************************
C                                                                *
C     HMTVAL evaluates functional values of a hermite spline of  *
C     fifth degree S(X) together with its first 5 derivatives at *
C     X=X0.                                                      *
C[BE*)
C                                                                *
C     While one could use this program to make an equidistant    *
C     value table for the function for plotting, this is not     *
C     advised since HMTVAL performs an expensive interval search *
C     for each value. Moreover one would not normally want to    *
C     find all derivatives. For a table of values we recommend   *
C     programs like HMTAB instead.                               *
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     N       number of nodes X(I)                               *
C     X0      x-value whose function value is desired            *
C     A       )                                                  *
C     B       )  vectors ..(1:N);                                *
C     C       )  the coefficients of the polynomial spline       *
C     D       )                                                  *
C     E       )                                                  *
C     F       )                                                  *
C     X       vector X(1:N); the nodes  X(1), ..., X(N)          *
C                                                                *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     OUTP    vector OUTP(1:5) for the derivatives:              *
C                        (K)                                     *
C             OUTP(K) = S    (X0)  for K = 1, ..., 5             *
C     HMTVAL functional value of the hermite spline of fifth     *
C            degree S(X).                                        *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : Elmar Pohl                                         *
C  Date     : 09.28.1985                                         *
C  Source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION X(N),A(N),B(N),C(N),D(N),E(N),F(N),OUTP(5)
      I=1
      K=N
    1 M=(I+K)/2
      IF (M .EQ. I) GOTO 5
         IF (X0 .LT. X(M)) GOTO 2
            I=M
            GOTO 1
    2       K=M
            GOTO 1
    5 X1=X0-X(I)
      HMTVAL=(((F(I)*X1+E(I))*X1+D(I))*X1+C(I))*X1+B(I)
      HMTVAL=HMTVAL*X1+A(I)
      OUTP(1)=((5.0D0*F(I)*X1+4.0D0*E(I))*X1+3.0D0*D(I))*X1+
     &                                       2.0D0*C(I)
      OUTP(1)=OUTP(1)*X1+B(I)
      OUTP(2)=((20.0D0*F(I)*X1+12.0D0*E(I))*X1+6.0D0*D(I))*X1+
     &                                         2.0D0*C(I)
      OUTP(3)=(60.0D0*F(I)*X1+24.0D0*E(I))*X1+6.0D0*D(I)
      OUTP(4)=120.0D0*F(I)*X1+24.0D0*E(I)
      OUTP(5)=120.0D0*F(I)
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE TRIDIG (N,A,B,C,D,IREP,X,ISING,GAMMA,ALPHA,G)
C[IX{TRIDIG}*)
C
C*****************************************************************
C                                                                *
C     TRIDIG solves a tridiagonal linear system of equations by  *
C     using a condensed Gauss algorithm.                         *
C[BE*)
C     If one wants to solve several systems for different right  *
C     hand sides, one must set IREP different from 0 for the     *
C     second and subsequent calls of TRIDIG with the same system *
C     matrix. Then the matrix will not be recomputed nor factored*
C     again for a saving of 2N-2 operations.                     *
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     N       size of the matrix                                 *
C     A       vector A(1:N);  right hand side                    *
C     B       vector B(1:N);  lower codiagonal                   *
C     C       vector C(1:N);  upper codiagonal                   *
C     D       vector D(1:N);  main diagonal                      *
C     IREP    If IREP differs from 0, then the program assumes   *
C             that it is in a repeat call with the same system   *
C             matrix. The program does not reconfigure the system*
C             matrix and its factors in this case.               *
C                                                                *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     X       vector X(1:N); the solution                        *
C     ISING   )  1, if the matrix is numerically singular        *
C             )     (no solution computed)                       *
C             )  0, otherwise                                    *
C                                                                *
C                                                                *
C     AUXILIARY PARAMETERS:                                      *
C     =====================                                      *
C     GAMMA   ) N-vectors ..(1:N);                               *
C     ALPHA   ) for repeated calls with IREP.NE.0 the vectors    *
C     G       ) GAMMA and ALPHA must nit be altered.             *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  Subroutines required: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : Elmar Pohl                                         *
C  Date     : 09.28.1985                                         *
C  Source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION A(N),B(N),C(N),D(N),X(N),GAMMA(N),ALPHA(N),
     +                 G(N)
      IF (IREP .NE. 0) GOTO 11
         ISING=0
         NM1=N-1
         ALPHA(1)=D(1)
         IF (ALPHA(1) .EQ. 0.0D0) GOTO 6
         IF(N .EQ. 1) GOTO 2
         DO 1 I=1,NM1
              GAMMA(I)=C(I)/ALPHA(I)
              ALPH=D(I+1)-B(I+1)*GAMMA(I)
              IF(ALPH .EQ. 0.0D0) GOTO 6
              ALPHA(I+1)=ALPH
    1    CONTINUE
   11 IF(N .GT. 1) GOTO 3
    2    X(1)=A(1)/ALPHA(1)
         RETURN
    3 G(1)=A(1)/ALPHA(1)
      DO 4 I=2,N
           G(I)=(A(I)-B(I)*G(I-1))/ALPHA(I)
    4 CONTINUE
      X(N)=G(N)
      NM1=N-1
      DO 5 I=1,NM1
           K=N-I
           X(K)=G(K)-GAMMA(K)*X(K+1)
    5 CONTINUE
      RETURN
    6 ISING=1
      RETURN
      END
