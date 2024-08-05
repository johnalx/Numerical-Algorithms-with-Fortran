C[BA*)
C[LE*)
      SUBROUTINE CEPSPM(AP,N,B,KPVT,RCOND,Z,WK)
C[IX{CEPSPM}*)
C
C*****************************************************************
C                                                                *
C     Condition estimate for a symmetric matrix AP which is given*
C     in factored and condensed form from SUBROUTINE ZSPMMK.     *
C[BE*)
C     If a condition estimate is not required, the SUBROUTINE    *
C     ZSPMOK is more time efficient.                             *
C     In order to solve A*X = B, a subsequent call of SESSPM     *
C     is necessary with WK as input matrix.                      *
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     AP    DOUBLE PRECISION vector AP(1:N*(N+1)/2), that        *
C           contains the symmetric matrix A in condensed form.   *
C           The columns of the upper triangle of A are stored in *
C           sequence in the vector AP as subvectors of length    *
C           N*(N+1)/2                                            *
C     N     dimension of the matrix A                            *
C     B     DOUBLE PRECISION vector B(1:N) containing the        *
C           right-hand side of the system of equations A*X = B.  *
C           The right-hand side is required for the SUBROUTINE   *
C           ZSPMMK, since the decomposition and possible alter-  *
C           ation of the right-hand side is performed in this    *
C           subroutine in order to save arithmetic operations.   *
C                                                                *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     AP    DOUBLE PRECISION vector AP(1:N*(N+1)/2), representing*
C           a block diagonal matrix that contains the factors of *
C           a decomposition in condensed form. The decomposition *
C           is given as  A = U*D*TRANS(U). Here U is the product *
C           of the permutation matrix and a unit upper triangular*
C           matrix, TRANS(U) denotes the transpose of U and D is *
C           a block diagonal matrix composed of 1 x 1 and 2 x 2  *
C           blocks                                               *
C     KPVT  INTEGER vector KPVT(1:N) containing the PIVOT indices*
C     RCOND DOUBLE PRECISION estimate of the reciprocal condi-   *
C           tion number of A.  For the system A*X = B, with      *
C           relative input errors of size EPSILON in A and B,    *
C           the relative error in the solution X will have the   *
C           size EPSILON/RCOND.                                  *
C           If RCOND is smaller than the machine constant, then  *
C           the matrix A is numerically singular.                *
C     Z     DOUBLE PRECISION auxiliary vector Z(1:N). Usually    *
C           the contents of Z has no significance. In case A is  *
C           close to being singular then Z is an approximated    *
C           null-vector for A, i.e.,                             *
C                   NORM(A*Z) = RCOND*NORM(A)*NORM(Z).           *
C     B     DOUBLE PRECISION vector B(1:N) containing the right- *
C           side of the system of equations  A*X = B  for use in *
C           the SUBROUTINE SESSPM in the necessary form.         *
C     WK    DOUBLE PRECISION vector WK(1:N*(N+1)/2)  (input      *
C           parameter for SESSPM), auxiliary vector containing   *
C           entries that are required for solving the system of  *
C           equations  A*X = B. This way the solution of  A*X = B*
C           can be determined in one, instead of two elimination *
C           steps by the SUBROUTINE SESSPM.                      *
C                                                                *
C                                                                *
C     condensed form                                             *
C                                                                *
C           The following code condenses the upper triangular    *
C           part of a symmetric matrix A to a vector AP          *
C                                                                *
C                       K = 0                                    *
C                       DO 20 J = 1, N                           *
C                          DO 10 I= 1, N                         *
C                             K = K + 1                          *
C                             AP(K) = A(I,J)                     *
C                    10    CONTINUE                              *
C                    20 CONTINUE                                 *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C     subroutines required: ZSPMMK, PCOSOL, PCOLTG, VECADD,      *
C                           VECMWC, SCAPRO, ABSSUM, INDMAX,      *
C                           VECXCH                               *
C                                                                *
C                                                                *
C     source : Linpack User's Guide, SIAM, Philadelphia, 1979    *
C                                                                *
C              The source was converted to FORTRAN 77. In some   *
C              instances it had to be modified and adjusted for  *
C              the requirements of our specific calling programs.*
C              This program and the corresponding subprograms    *
C              are not compatible with the originals from the    *
C              Linpack User's Guide.                             *
C                                                                *
C*****************************************************************
C                                                                *
C     authors :  Michael Groenheim, Ina Hinze                    *
C     date    :  10.25.1989                                      *
C     source  :  FORTRAN 77                                      *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      INTEGER N,KPVT(N)
      DOUBLE PRECISION AP(N*(N+1)/2),Z(N),B(N),WK(N*(N+1)/2)
      DOUBLE PRECISION RCOND
      DOUBLE PRECISION EK,ANORM,S,ABSSUM,YNORM
      INTEGER I,IJ,IERR,J,JM1,J1
C
C     determine the norm of A making use of symmetry
C
      J1 = 1
      DO 20 J = 1, N
         Z(J) = ABSSUM(J,AP(J1))
         IJ = J1
         J1 = J1 + J
         JM1 = J - 1
         IF (JM1 .GE. 1) THEN
            DO 10 I = 1, JM1
               Z(I) = Z(I) + DABS(AP(IJ))
               IJ = IJ + 1
   10       CONTINUE
         ENDIF
   20 CONTINUE
      ANORM = 0.0D0
      DO 30 J = 1, N
         ANORM = DMAX1(ANORM,Z(J))
   30 CONTINUE
      IF (ANORM .EQ. 0.0D0) THEN
         RCOND = 0.0D0
         RETURN
      ENDIF
C
C     factor the matrix by using SUBROUTINE ZSPMMK
C
      CALL ZSPMMK(AP,N,B,KPVT,WK,IERR)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF THE NORM(INVERSE(A)))).
C     ESTIMATE = NORM(Z)/NORM(Y). Here A*Z = Y and A*Y = E.
C     The elements of E are chosen in such a way that the elements
C     of W become as large as possible, where U*D*W = E.
C
      EK = 1.0D0
      DO 40 J = 1, N
         Z(J) = 0.0D0
   40 CONTINUE
C
C     solve U*D*W = E
C
      CALL PCOSOL (1,AP,N,Z,KPVT,S,EK,YNORM)
C
C     solve TRANS(U) * Y = W
C
      CALL PCOLTG(AP,N,Z,KPVT,S)
      S = 1.0D0/ABSSUM(N,Z)
      CALL VECMWC(N,S,Z)
C
      YNORM = 1.0D0
C
C     solve U*D*V = Y
C
      CALL PCOSOL (2,AP,N,Z,KPVT,S,EK,YNORM)
C
C     solve TRANS(U) * Z = V
C
      CALL PCOLTG(AP,N,Z,KPVT,S)
      YNORM = S * YNORM
C
C     set ZNORM = 1.0
C
      S = 1.0D0/ABSSUM(N,Z)
      CALL VECMWC(N,S,Z)
      YNORM = S*YNORM
C
      RCOND = YNORM/ANORM
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE ZSPMMK (AP,N,B,KPVT,WK,IERR)
C[IX{ZSPMMK}*)
C
C*****************************************************************
C                                                                *
C     Factorization of a symmetric matrix, given in condensed    *
C     format, by elimination that relies on symmetric pivoting.  *
C[BE*)
C     If a condition estimate is not required, the SUBROUTINE    *
C     ZSPMOK will perform the task faster.                       *
C     In order to solve A*X = B the subsequent call of SESSPM is *
C     necessary.                                                 *
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     AP    DOUBLE PRECISION vector AP(1:(N*(N+1)/2)), contain-  *
C           ing the symmetric matrix A in condensed form.        *
C           The columns of its upper triangle are stored sequen- *
C           tially in the vector AP of length N*(N+1)/2          *
C     N     dimension of the matrix A                            *
C     B     DOUBLE PRECISION vector B(1:N), the right-hand side. *
C           The right-hand side is required in SUBROUTINE ZSPMMK *
C           in order to save arithmetic operations.              *
C                                                                *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     AP    DOUBLE PRECISION vector AP(1:N*(N+1)/2), representing*
C           a block diagonal matrix that contains the factors of *
C           a decomposition in condensed form. The decomposition *
C           is given as  A = U*D*TRANS(U). Here U is the product *
C           of the permutation matrix and a unit upper triangular*
C           matrix, TRANS(U) denotes the transpose of U and D is *
C           a block diagonal matrix composed of 1 x 1 and 2 x 2  *
C           blocks                                               *
C     KPVT  INTEGER vector KPVT(1:N) containing the PIVOT indices*
C     B     DOUBLE PRECISION vector B(1:N) containing the right- *
C           side of the system of equations  A*X = B  for use in *
C           the SUBROUTINE SESSPM in the necessary form.         *
C     WK    DOUBLE PRECISION vector WK(1:N*(N+1)/2)  (input      *
C           parameter for SESSPM), auxiliary vector containing   *
C           entries that are required for solving the system of  *
C           equations  A*X = B. This way the solution of  A*X = B*
C           can be determined in one, instead of two elimination *
C           steps by the SUBROUTINE SESSPM.                      *
C     IERR  error parameter                                      *
C           = 0, everything is o.k                               *
C           = K, the K-th PIVOT block is singular.               *
C                For the current subroutine, this does not denote*
C                an error; however, it indicates that the        *
C                SUBROUTINE SESSPM may encounter division by zero*
C                                                                *
C                                                                *
C     condensed form                                             *
C                                                                *
C           The following code condenses the upper triangular    *
C           part of a symmetric matrix A to a vector AP          *
C                                                                *
C                       K = 0                                    *
C                       DO 20 J = 1, N                           *
C                          DO 10 I= 1, N                         *
C                             K = K + 1                          *
C                             AP(K) = A(I,J)                     *
C                    10    CONTINUE                              *
C                    20 CONTINUE                                 *
C                    20 CONTINUE                                 *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C     subroutines required: VECADD, VECXCH, INDMAX               *
C                                                                *
C                                                                *
C     source : Linpack User's Guide, SIAM, Philadelphia, 1979    *
C                                                                *
C*****************************************************************
C                                                                *
C     authors :  Michael Groenheim, Ina Hinze                    *
C     date    :  10.25.1989                                      *
C     source  :  FORTRAN 77                                      *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      INTEGER N,KPVT(N),IERR
      DOUBLE PRECISION AP(N*(N+1)/2),B(N),WK(N*(N+1)/2)
      DOUBLE PRECISION D,D1,D2,T
      DOUBLE PRECISION ABSAKK,ALPHA,COLMAX,ROWMAX
      INTEGER INDMAX,IJ,IJJ,IK,IKM1,IM,IMAX,IMAXP1,IMIM,IMJ,IMK
      INTEGER J,JJ,JK,JKM1,JMAX,JMIM,K,KK,KM1,KM1K,KM1KM1,KM2,KSTEP
      LOGICAL SWAP
C
C     store condensed vector AP in auxiliary vector WK
C
      DO 5 I = 1,N*(N+1)/2
         WK(I) = AP(I)
5     CONTINUE
C
C     initialization
C
C     ALPHA is used to determine the PIVOT block size
C
      ALPHA = (1.0D0 + SQRT(17.0D0))/8.0D0
C
      IERR = 0
C
C     main loop over K; K runs from N back to 1
C
      K = N
      IK = (N*(N - 1))/2
10    CONTINUE
C
C     leave the loop if K=0 or K=1
C
      IF (K .EQ. 0) RETURN
      IF (K .LE. 1) THEN
         KPVT(1) = 1
         IF (AP(1) .EQ. 0.0D0) IERR = 1
         RETURN
      ENDIF
C
C     this part of the program determines the elimination method to be
C     used. After this part has been executed, KSTEP is set to be the size
C     of the PIVOT block and SWAP is set to .TRUE., if swapping is
C     necessary.
C
      KM1 = K - 1
      KK = IK + K
      ABSAKK = DABS(AP(KK))
C
C     determine the largest off-diagonal element in
C     magnitude in column K
C
      IMAX = INDMAX(K-1,AP(IK+1))
      IMK = IK + IMAX
      COLMAX = DABS(AP(IMK))
      IF (ABSAKK .LT. ALPHA*COLMAX) THEN
C
C        determine the largest off-diagonal element
C        in magnitude in row IMAX
C
         ROWMAX = 0.0D0
         IMAXP1 = IMAX + 1
         IM = IMAX*(IMAX - 1)/2
         IMJ = IM + 2*IMAX
         DO 20 J = IMAXP1, K
            ROWMAX = DMAX1(ROWMAX,DABS(AP(IMJ)))
            IMJ = IMJ + J
20       CONTINUE
         IF (IMAX .NE. 1) THEN
            JMAX = INDMAX(IMAX-1,AP(IM+1))
            JMIM = JMAX + IM
            ROWMAX = DMAX1(ROWMAX,DABS(AP(JMIM)))
         ENDIF
         IMIM = IMAX + IM
         IF (DABS(AP(IMIM)) .LT. ALPHA*ROWMAX) THEN
            IF (ABSAKK .LT. ALPHA*COLMAX*(COLMAX/ROWMAX)) THEN
               KSTEP = 2
               SWAP = IMAX .NE. KM1
            ELSE
               KSTEP = 1
               SWAP = .FALSE.
            ENDIF
         ELSE
            KSTEP = 1
            SWAP = .TRUE.
         ENDIF
      ELSE
         KSTEP = 1
         SWAP = .FALSE.
      ENDIF
      IF (DMAX1(ABSAKK,COLMAX) .EQ. 0.0D0) THEN
C
C        column K is the zero column. Modify IERR and reiterate loop
C
         KPVT(K) = K
         IERR = K
         IK = IK - (K - 1)
         IF (KSTEP .EQ. 2) IK = IK - (K - 2)
         K = K - KSTEP
         GOTO 10
      ENDIF
      IF (KSTEP .EQ. 2) THEN
C
C        2 x 2 PIVOT block
C
         KM1K = IK + K - 1
         IKM1 = IK - (K - 1)
         IF (SWAP) THEN
C
C           perform swap
C
            CALL VECXCH(IMAX,AP(IM+1),AP(IKM1+1))
            CALL VECXCH(IMAX,WK(IM+1),WK(IKM1+1))
            IMJ = IKM1 + IMAX
            DO 30 JJ = IMAX, KM1
               J = KM1 + IMAX - JJ
               JKM1 = IKM1 + J
               T = AP(JKM1)
               T1 = WK(JKM1)
               AP(JKM1) = AP(IMJ)
               WK(JKM1) = WK(IMJ)
               AP(IMJ) = T
               WK(IMJ) = T1
               IMJ = IMJ - (J - 1)
30          CONTINUE
            T = AP(KM1K)
            AP(KM1K) = AP(IMK)
            AP(IMK) = T
            T = WK(KM1K)
            WK(KM1K) = WK(IMK)
            WK(IMK) = T
            T = B(K-1)
            B(K-1) = B(IMAX)
            B(IMAX) = T
         ENDIF

C
C        perform elimination
C
         KM2 = K - 2
         IF (KM2 .NE. 0) THEN
            KM1KM1 = IKM1 + K - 1
            D = AP(KM1K) * AP(KM1K) - AP(KK) * AP (KM1KM1)
            IJ = IK - (K - 1) - (K - 2)
            DO 40 JJ = 1, KM2
               J = KM1 - JJ
               JK = IK + J
               JKM1 = IKM1 + J
               D1 = (AP(KM1KM1) * AP(JK) - AP(JKM1) * AP(KM1K)) / D
               D2 = (AP(KK) * AP(JKM1) - AP(JK) * AP(KM1K)) / D
               CALL VECADD(J,D1,AP(IK+1),AP(IJ+1))
               CALL VECADD(1,D1,B(K),B(J))
               CALL VECADD(J,D2,AP(IKM1+1),AP(IJ+1))
               CALL VECADD(1,D2,B(K-1),B(J))
               IF (IKM1 .EQ. 1) WK(IKM1) = AP(IKM1)
               DO 45 ID = 1 , J
                  WK(IJ+ID) = AP(IJ+ID)
45             CONTINUE
               AP(JK) = D1
               AP(JKM1) = D2
               IJJ = IJ + J
               IJ = IJ - (J - 1)
40          CONTINUE
         ENDIF
C
C        set up PIVOT vector
C
         KPVT(K) = 1 - K
         IF (SWAP) KPVT(K) = -IMAX
         KPVT(K-1) = KPVT(K)
      ELSE
C
C        1 x 1 PIVOT block
C
         IF (SWAP) THEN
C
C           perform swap
C
            CALL VECXCH(IMAX,AP(IM+1),AP(IK+1))
            CALL VECXCH(IMAX,WK(IM+1),WK(IK+1))
            IMJ = IK + IMAX
            DO 50 JJ = IMAX, K
               J = K + IMAX - JJ
               JK = IK + J
               T = AP(JK)
               AP(JK) = AP(IMJ)
               AP(IMJ) = T
               T = WK(JK)
               WK(JK) = WK(IMJ)
               WK(IMJ) = T
               IMJ = IMJ - (J - 1)
50          CONTINUE
            T = B(K)
            B(K) = B(IMAX)
            B(IMAX) = T
         ENDIF
C
C        perform elimination
C
         IJ = IK - (K - 1)
         DO 60 JJ = 1, KM1
            J = K - JJ
            JK = IK + J
            D1 = -AP(JK) / AP(KK)
            CALL VECADD(J,D1,AP(IK+1),AP(IJ+1))
            CALL VECADD(1,D1,B(K),B(J))
            DO 70 ID = 1 , J
               WK(IJ+ID) = AP(IJ+ID)
70          CONTINUE
            IJJ = IJ + J
            AP(JK) = D1
            IJ = IJ - (J - 1)
60       CONTINUE
C
C        set up PIVOT vector
C
         KPVT(K) = K
         IF (SWAP) KPVT(K) = IMAX
      ENDIF
      IK = IK - (K - 1)
      IF (KSTEP .EQ. 2) IK = IK - (K - 2)
      K = K - KSTEP
      GOTO 10
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE VECADD (N,SA,SX,SY)
C[IX{VECADD}*)
C
C*****************************************************************
C                                                                *
C     Multiplies the vector SX by the constant SA and then adds  *
C     the vector SX to the vector SY.                            *
C                                                                *
C----------------------------------------------------------------*
C[BE*)
C                                                                *
C     subroutines required: none                                 *
C                                                                *
C*****************************************************************
C                                                                *
C     authors :  Michael Groenheim, Ina Hinze                    *
C     date    :  10.25.1989                                      *
C     source  :  FORTRAN 77                                      *
C                                                                *
C*****************************************************************
C
      DOUBLE PRECISION SX(N),SY(N),SA
      INTEGER I,N
C
      IF (N .LE. 0) RETURN
C
      IF (SA .EQ. 0.0D0) RETURN
C
      DO 10 I = 1,N
         SY(I) = SY(I) + SA * SX(I)
10    CONTINUE
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE VECXCH (N,SX,SY)
C[IX{VECXCH}*)
C
C*****************************************************************
C                                                                *
C     Swaps the vectors SX and SY.                               *
C                                                                *
C----------------------------------------------------------------*
C[BE*)
C                                                                *
C     subroutines required: none                                 *
C                                                                *
C*****************************************************************
C                                                                *
C     authors :  Michael Groenheim, Ina Hinze                    *
C     date    :  10.25.1989                                      *
C     source  :  FORTRAN 77                                      *
C                                                                *
C*****************************************************************
C
      DOUBLE PRECISION SX(N),SY(N),DUMMY
      INTEGER  I,N
C
      DO 10 I = 1,N
         DUMMY  = SX(I)
         SX(I) = SY(I)
         SY(I) = DUMMY
   10 CONTINUE
      RETURN
      END
C
C
C[BA*)
C[LE*)
      INTEGER FUNCTION INDMAX (N,SX)
C[IX{INDMAX}*)
C
C*****************************************************************
C                                                                *
C     Determines the index of the element of the vector SX that  *
C     has the largest entry in magnitude.                        *
C                                                                *
C----------------------------------------------------------------*
C[BE*)
C                                                                *
C     subroutines required: none                                 *
C                                                                *
C*****************************************************************
C                                                                *
C     authors :  Michael Groenheim, Ina Hinze                    *
C     date    :  10.25.1989                                      *
C     source  :  FORTRAN 77                                      *
C                                                                *
C*****************************************************************
C
      DOUBLE PRECISION  SX(N),SMAX
      INTEGER  I,N
C
      INDMAX = 0
C
      IF (N .LT. 1) RETURN
C
      INDMAX = 1
C
      IF (N .EQ. 1) RETURN
C
      SMAX = DABS (SX(1))
C
      DO 10 I = 2,N
         IF (DABS (SX(I)) .GT. SMAX ) THEN
            INDMAX = I
            SMAX = DABS (SX(I))
         ENDIF
   10 CONTINUE
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE PCOSOL (IPOINT,AP,N,Z,KPVT,S,EK,YNORM)
C[IX{PCOSOL}*)
C
C*****************************************************************
C                                                                *
C     Subroutine of CEPSPM.                                      *
C     This subroutine helps to solve a system of equations       *
C                   U*D*W  = E  or  U*D*V = Y.                   *
C     It is called twice by ZSPM.. .                             *
C[BE*)
C     The method of solving the                                  *
C     equations above is identical except for small variations.  *
C     In order to combine both tasks in one subroutine, we       *
C     introduce the parameter IPOINT, which, depending on the    *
C     kind of call, initiates branches at the relevant places in *
C     this subroutine.                                           *
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     IPOINT  flag                                               *
C             = 1  , in case PCOSOL is called for solving        *
C                    U*D*W = E                                   *
C             = 2  , in case PCOSOL is called for solving        *
C                    U*D*V = Y                                   *
C     N       dimension of the matrix A                          *
C     KPVT    INTEGER vector KVPT(1:N) containing the PIVOT      *
C             indices                                            *
C     Z       DOUBLE PRECISION auxiliary vector Z(1:N), needed   *
C             for solving  U*D*W = E  or  U*D*V = Y.             *
C     AP      vector containing the factors of the decomposition *
C             of the symmetric matrix A in condensed form.       *
C     S       DOUBLE PRECISION auxiliary variable                *
C     EK      DOUBLE PRECISION auxiliary variable                *
C                                                                *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     S       DOUBLE PRECISION auxiliary variable                *
C     YNORM   DOUBLE PRECISION norm of Y                         *
C     EK      DOUBLE PRECISION auxiliary variable                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C     subroutines required: VECADD, VECMWC                       *
C                                                                *
C                                                                *
C     source : Linpack User's Guide, SIAM, Philadelphia, 1979    *
C                                                                *
C*****************************************************************
C                                                                *
C     authors :  Michael Groenheim, Ina Hinze                    *
C     date    :  10.25.1989                                      *
C     source  :  FORTRAN 77                                      *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      INTEGER N,IPOINT,KPVT(N)
      DOUBLE PRECISION AP(N*(N+1)/2),Z(N)
      DOUBLE PRECISION D,D1,D2,EK,T
      DOUBLE PRECISION S,YNORM
      INTEGER IK,IKM1,K,KK,KM1K,KM1KM1,KP,KPS,KSTEP
      K = N
      IK = N*(N - 1)/2
   10 IF (K .EQ. 0) RETURN
      KK = IK + K
      IKM1 = IK - (K - 1)
      KSTEP = 1
      IF (KPVT(K) .LT. 0) KSTEP = 2
      IF (IPOINT .EQ. 1) THEN
C
C        call was for solving U*D*W = E
C
         KP = IABS(KPVT(K))
         KPS = K + 1 - KSTEP
         IF (KP .NE. KPS) THEN
            T = Z(KPS)
            Z(KPS) = Z(KP)
            Z(KP) = T
         ENDIF
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,Z(K))
         Z(K) = Z(K) + EK
         CALL VECADD(K-KSTEP,Z(K),AP(IK+1),Z(1))
      ENDIF
      IF (IPOINT .EQ. 2) THEN
C
C        call was for solving U*D*V = Y
C
         IF (K .NE. KSTEP) THEN
            KP = IABS(KPVT(K))
            KPS = K + 1 - KSTEP
            IF (KP .NE. KPS) THEN
               T = Z(KPS)
               Z(KPS) = Z(KP)
               Z(KP) = T
            ENDIF
            CALL VECADD(K-KSTEP,Z(K),AP(IK+1),Z(1))
            IF (KSTEP .EQ. 2) THEN
               CALL VECADD(K-KSTEP,Z(K-1),AP(IKM1+1),Z(1))
            ENDIF
         ENDIF
      ENDIF
C
C     1 x 1 PIVOT block
C
      IF (KSTEP .EQ. 1) THEN
         IF (DABS(Z(K)) .GT. DABS(AP(KK))) THEN
            S = DABS(AP(KK))/DABS(Z(K))
            CALL VECMWC (N,S,Z)
            IF (IPOINT .EQ. 1) EK = S * EK
            IF (IPOINT .EQ. 2) YNORM = S * YNORM
         ENDIF
         IF (AP(KK) .NE. 0.0D0) Z(K) = Z(K)/AP(KK)
         IF (AP(KK) .EQ. 0.0D0) Z(K) = 1.0D0
         K = K - KSTEP
         IK = IK - K
         GOTO 10
      ENDIF
C
C     2 x 2 PIVOT block
C
      IF (KSTEP .EQ. 2) THEN
         IF (IPOINT .EQ. 1) THEN
            IF (Z(K-1) .NE. 0.0D0) EK = DSIGN(EK,Z(K-1))
            Z(K-1) = Z(K-1) + EK
            CALL VECADD(K-KSTEP,Z(K-1),AP(IKM1+1),Z(1))
         ENDIF
         KM1K = IK + K - 1
         KM1KM1 = IKM1 + K - 1
         D1 = AP(KM1KM1) * Z(K) - AP(KM1K) * Z(K-1)
         D2 = AP(KK) * Z(K-1) - AP(KM1K) * Z(K)
         D  = AP(KK) * AP(KM1KM1) - AP(KM1K) * AP(KM1K)
         Z(K)   = D1 / D
         Z(K-1) = D2 / D
         K = K - KSTEP
         IK = IK - K - (K + 1)
         GOTO 10
      ENDIF
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE VECMWC (N,SA,SX)
C[IX{VECMWC}*)
C
C*****************************************************************
C                                                                *
C     Multiplies the vector SX by the constant SA.               *
C                                                                *
C----------------------------------------------------------------*
C[BE*)
C                                                                *
C     subroutines required: none                                 *
C                                                                *
C*****************************************************************
C                                                                *
C     authors :  Michael Groenheim, Ina Hinze                    *
C     date    :  10.25.1989                                      *
C     source  :  FORTRAN 77                                      *
C                                                                *
C*****************************************************************
C
      DOUBLE PRECISION SA,SX(N)
      INTEGER I,N
C
      IF (N. LE. 0) RETURN
C
      DO 10 I = 1,N
         SX(I) = SA * SX(I)
   10 CONTINUE
      RETURN
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE PCOLTG(AP,N,Z,KPVT,S)
C[IX{PCOLTG}*)
C
C*****************************************************************
C                                                                *
C     Subroutine of CEPSPM.                                      *
C     This subroutine helps to solve the system the equations    *
C     TRANS(U) * Y = W   or  TRANS(U) * Z = V.                   *
C[BE*)
C                                                                *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     N     dimension of the matrix A                            *
C     Z     DOUBLE PRECISION auxiliary vector Z(1:N), needed     *
C           to solve TRANS(U) * Y = W  or  TRANS(U) * Z = V      *
C     KPVT  INTEGER vector KPVT(1:N) containing the PIVOT        *
C           indices                                              *
C     AP    vector containing the factors of the decomposition   *
C           of the symmetric matrix A in condensed form          *
C                                                                *
C                                                                *
C     OUTPUT PARAMETER:                                          *
C     =================                                          *
C     S     DOUBLE PRECISION auxiliary variable                  *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C     subroutines required: SCAPRO, VECMWC, ABSSUM               *
C                                                                *
C                                                                *
C     source : Linpack User's Guide, SIAM, Philadelphia, 1979    *
C                                                                *
C*****************************************************************
C                                                                *
C     authors :  Michael Groenheim, Ina Hinze                    *
C     date    :  10.25.1989                                      *
C     source  :  FORTRAN 77                                      *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      INTEGER N,KPVT(N)
      DOUBLE PRECISION AP(N*(N+1)/2),Z(N)
      DOUBLE PRECISION SCAPRO,T,S,ABSSUM
      INTEGER IK,IKP1,K,KP,KSTEP
      S = 1.0D0 / ABSSUM(N,Z)
      CALL VECMWC(N,S,Z)
      K = 1
      IK = 0
   10 IF (K .LE. N) THEN
         KSTEP = 1
         IF (KPVT(K) .LT. 0) KSTEP = 2
         IF (K .NE. 1) THEN
            Z(K) = Z(K) + SCAPRO(K-1,AP(IK+1),Z(1))
            IKP1 = IK + K
            IF (KSTEP .EQ. 2)
     *              Z(K+1) = Z(K+1) + SCAPRO(K-1,AP(IKP1+1),Z(1))
            KP = IABS(KPVT(K))
            IF (KP .NE. K) THEN
               T = Z(K)
               Z(K) = Z(KP)
               Z(KP) = T
            ENDIF
         ENDIF
         IK = IK + K
         IF (KSTEP .EQ. 2) IK = IK + (K + 1)
         K = K + KSTEP
         GOTO 10
      ENDIF
      RETURN
      END
C
C
C[BA*)
C[LE*)
      DOUBLE PRECISION FUNCTION SCAPRO (N,SX,SY)
C[IX{SCAPRO}*)
C
C*****************************************************************
C                                                                *
C     Determines the scalar product of two vectors SX and SY.    *
C                                                                *
C----------------------------------------------------------------*
C[BE*)
C                                                                *
C     subroutines required: none                                 *
C                                                                *
C*****************************************************************
C                                                                *
C     authors :  Michael Groenheim, Ina Hinze                    *
C     date    :  10.25.1989                                      *
C     source  :  FORTRAN 77                                      *
C                                                                *
C*****************************************************************
C
      DOUBLE PRECISION SX(N),SY(N)
      INTEGER I,N
C
      SCAPRO = 0.0D0
C
      DO 10 I = 1,N
         SCAPRO = SCAPRO + SX(I) * SY(I)
   10 CONTINUE
      RETURN
      END
C
C
C[BA*)
C[LE*)
      DOUBLE PRECISION FUNCTION ABSSUM (N,SX)
C[IX{ABSSUM}*)
C
C*****************************************************************
C                                                                *
C     Forms the sum of the absolute values of the entries in SX. *
C                                                                *
C----------------------------------------------------------------*
C[BE*)
C                                                                *
C     subroutine required: none                                  *
C                                                                *
C*****************************************************************
C                                                                *
C     authors :  Michael Groenheim, Ina Hinze                    *
C     date    :  10.25.1989                                      *
C     source  :  FORTRAN 77                                      *
C                                                                *
C*****************************************************************
C
      DOUBLE PRECISION SX(N)
      INTEGER I,N
C
      ABSSUM = 0.0D0
C
      IF (N .LE. 0) RETURN
C
      DO 10 I=1,N
         ABSSUM = ABSSUM + DABS(SX(I))
   10 CONTINUE
      RETURN
      END
