C[BA*)
C[LE*)
C[LE*)
C[LE*)
C[FE{F 7.8}
C[  {QR Algorithm}
C[  {Eigenvalues and Eigenvectors of a Matrix via the
C[   QR Algorithm}*)
C[LE*)
      INTEGER FUNCTION EIGEN (VEC,ORTHO,EVNORM,BASIS,LD,N,MAT,
C[IX{EIGEN}*)
     *                        SCAL,D,EIVEC,VALR,VALI,CNT,LOW,HIGH)
C
C*****************************************************************
C                                                                *
C     PROGRAM OBJECTIVE:                                         *
C     ==================                                         *
C     This FUNCTION-subroutine of type INTEGER determines all    *
C     eigenvalues and eigenvectors of a real square matrix MAT.  *
C[BE*)
C     The eigenvalues are stored in the vectors VALR(1:N)        *
C     (real parts) and VALI(1:N) (imaginary parts).              *
C     The eigenvectors are stored in the array EIVEC(1:N,1:N)    *
C     provided the flag VEC is set.                              *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     VEC:         parameter designating eigenvector computation:*
C                    =  .TRUE. : compute eigenvectors            *
C                    =  .FALSE.: compute eigenvalues only        *
C     ORTHO:       flag indicating desired reduction of MAT to   *
C                  Hessenberg form : if .TRUE. use orthogonal    *
C                  transformations in ORTHES; if .FALSE. use     *
C                  elementary Gauss eliminations via ELMHES which*
C                  is slightly faster. For symmetric matrices MAT*
C                  only ORTHES will preserve symmetry.           *
C     EVNORM:      flag that governs potential normalizing of    *
C                  eigenvectors if .TRUE.                        *
C     N:           order of the square matrix MAT                *
C     LD:          leading dimension of the matrix MAT and the   *
C                  vector EIVEC as defined in the calling program*
C     MAT:         2-dimensional array MAT(1:N,1:N) of type      *
C                  DOUBLE PRECISION, containing the real input   *
C                  matrix                                        *
C     BASIS:       basis of the floating-point representation    *
C                  used by the machine (in most cases 2 or 16)   *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     MAT:         the upper triangle of this (N,N) matrix       *
C                  contains the eigenvectors of the quasi-trian- *
C                  gular matrix, that is produced by the         *
C                  QR-method.                                    *
C     D:           N-vector with transform info from ORTHES      *
C     VALR,VALI:   two N-vectors VALR(1:N), VALI(1:N) of type    *
C                  DOUBLE PRECISION, that contain the real and   *
C                  imaginary parts of the eigenvalues of MAT     *
C     EIVEC:       2-dimensional array EIVEC(1:N,1:N) of type    *
C                  DOUBLE PRECISION, that contains the normalized*
C                  eigenvectors of the input matrix MAT as       *
C                  columns in case VEC = .TRUE. .                *
C                  If the I-th eigenvalue of MAT is real, the    *
C                  I-th column of EIVEC contains the correspond- *
C                  ing real eigenvector. If the I-th and (I+1)-st*
C                  eigenvalues are a complex conjugate pair,     *
C                  the I-th and (I+1)-st columns of EIVEC contain*
C                  the real and imaginary parts of the associated*
C                  eigenvector for the eigenvalue with positive  *
C                  imaginary part.                               *
C     CNT:         N-vector CNT(1:N) of type INTEGER containing  *
C                  the number of iteration steps for each eigen- *
C                  value. If two eigenvalues are found simul-    *
C                  taneously as a complex conjugate pair, the    *
C                  number of iteration steps used for both is    *
C                  recorded with a positive sign for the first   *
C                  and with a negative sign for the second       *
C                  eigenvalue.                                   *
C     SCAL:        N-vector SCAL(1:N) of type DOUBLE PRECISION   *
C                  containing information about the permutations *
C                  and the scaling factors used.                 *
C     LOW,HIGH:    The rows numbered 1 to LOW-1 or HIGH+1 to N   *
C                  contain isolated eigenvalues, i.e., eigen-    *
C                  values for the eigenvectors e_i, the unit     *
C                  vectors                                       *
C                                                                *
C     RETURN VALUES of FUNCTION EIGEN:                           *
C     ================================                           *
C     0:      no error, eigenvalue-eigenvector problem solved.   *
C     401:    order N of the input matrix MAT is less than 1.    *
C     402:    MAT is the zero matrix.                            *
C     403:    the maximum number of steps for the QR-method      *
C             has been reached. However, not all of the eigen-   *
C             values have be determined.                         *
C                                                                *
C     LOCAL VARIABLES:                                           *
C     ================                                           *
C     ONE,TWO,HALF: floating-point constants                     *
C     EPS:          machine constant                             *
C     TEMP:         auxiliary variable                           *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: BALAN, ELMHES, ELMTRA, HQR2, BALBAK     *
C                        NORMAL, SWAP, COMDIV, COMABS, ORTHES,   *
C                        ORTHTRA                                 *
C                                                                *
C                                                                *
C  sources : 1. Martin, R. S. and Wilkinson, J. H.,              *
C               see [MART70].                                    *
C            2. Parlett, B. N. and Reinsch, C., see [PARL69].    *
C            3. Peters, G. and Wilkinson, J. H., see [PETE70].   *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Juergen Dietel                                     *
C  date     : 7.15.1993                                          *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      INTEGER          BASIS,LD,N,CNT(1:N),LOW,HIGH
      LOGICAL          VEC,ORTHO,EVNORM
      DOUBLE PRECISION MAT(1:LD,1:N),SCAL(1:N),EIVEC(1:LD,1:N)
      DOUBLE PRECISION VALR(1:N),VALI(1:N),D(1:N)
      DOUBLE PRECISION ONE,TWO,HALF
      PARAMETER        (ONE = 1.0D0,TWO = 2.0D0,HALF = 0.5D0)
      INTEGER          RES,BALAN,ELMHES,ORTHES,ELMTRA,ORTTRA,HQR2
      INTEGER          BALBAK, NORMAL
      DOUBLE PRECISION EPS,TEMP
C
C     determine the machine constant EPS, i.e., the smallest
C     positive number for which the  1 + EPS > 1)  is true
C
      TEMP = TWO
      EPS = ONE
   10 IF (ONE .LT. TEMP) THEN
         EPS = EPS * HALF
         TEMP = ONE + EPS
         GOTO 10
         ENDIF
      EPS = TWO * EPS
      RES = BALAN(LD,N,MAT,SCAL,LOW,HIGH,BASIS)
      IF (RES .NE. 0) THEN
         EIGEN = RES + 100
         RETURN
      ENDIF
      IF (ORTHO) THEN
         RES = ORTHES(LD,N,LOW,HIGH,MAT,D,EPS)
      ELSE
         RES = ELMHES(LD,N,LOW,HIGH,MAT,CNT)
      ENDIF
      IF (RES .NE. 0) THEN
         EIGEN = RES + 200
         RETURN
      ENDIF
      IF (VEC) THEN
         IF (ORTHO) THEN
            RES = ORTTRA(LD,N,LOW,HIGH,MAT,D,EIVEC)
         ELSE
            RES = ELMTRA(LD,N,LOW,HIGH,MAT,CNT,EIVEC)
         ENDIF
         IF (RES .NE. 0) THEN
            EIGEN = RES + 300
            RETURN
         ENDIF
      ENDIF
      RES = HQR2(VEC,LD,N,LOW,HIGH,MAT,VALR,VALI,EIVEC,CNT,EPS)
      IF (RES .NE. 0) THEN
         EIGEN = RES + 400
         RETURN
      ENDIF
      IF (VEC) THEN
         RES = BALBAK(LD,N,LOW,HIGH,SCAL,EIVEC)
         IF (RES .NE. 0) THEN
            EIGEN = RES + 500
            RETURN
         ENDIF
         IF (EVNORM) THEN
            RES = NORMAL(LD,N,EIVEC,VALI)
            IF (RES .NE. 0) THEN
               EIGEN = RES + 600
               RETURN
            ENDIF
         ENDIF
      ENDIF
      EIGEN = 0
      END
C
C
C[BA*)
C[LE*)
      INTEGER FUNCTION BALAN (LD,N,MAT,SCAL,LOW,HIGH,BASIS)
C[IX{BALAN}*)
C
C*****************************************************************
C                                                                *
C     The procedure BALAN balances a given real matrix with      *
C     respect to the column sum norm.                            *
C[BE*)
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     LD:       leading dimension of the matrix as defined in    *
C               the calling program                              *
C     N:        the order of the given square matrix             *
C     MAT:      2-dimensional array MAT(1:N,1:N) containing the  *
C               input matrix                                     *
C     BASIS:    basis for the floating-point representation of   *
C               the machine                                      *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     MAT:      the balanced matrix                              *
C     LOW,HIGH: two INTEGER numbers, for which MAT(I,J) = 0      *
C               if the following holds:                          *
C               1. I > J, and                                    *
C               2. J = 1, ..., LOW-1 or I = HIGH+1, ..., N       *
C     SCAL:     N-vector SCAL(1:N) containing information about  *
C               permutations and scaling factors used.           *
C                                                                *
C     RETURN VALUE of SUBROUTINE BALAN:                          *
C     =================================                          *
C     0:             no error                                    *
C                                                                *
C     LOCAL VARIABLES:                                           *
C     ============== =                                           *
C     ZERO,ONE,PT95: floating-point constants                    *
C     I,J,K,L:       counting variables                          *
C     B2:            square of the machine floating point basis  *
C     R,C,F,G,S:     auxiliary variables for determining row     *
C                    norms, reciprocals etc.                     *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: SWAP                                    *
C                                                                *
C                                                                *
C  sources : Parlett, B. N. and Reinsch, C., see [PARL69].       *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Juergen Dietel                                     *
C  date     : 04.10.1987                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      INTEGER          LD,N,LOW,HIGH,BASIS
      DOUBLE PRECISION SCAL(1:N),MAT(1:LD,1:N)
      DOUBLE PRECISION ZERO,ONE,PT95
      PARAMETER        (ZERO = 0.0D0,ONE = 1.0D0,PT95 = 0.95D0)
      INTEGER          I,J,K,L,B2
      DOUBLE PRECISION R,C,F,G,S
C
C     reduce the norm of MAT(1:N,1:N) by exact diagonal similarity
C     transformations and store in SCAL(1:N)
C
      B2 = BASIS*BASIS
      L = 1
      K = N
C
C     search for rows isolating an eigenvalue and
C     move them to the bottom
C
   10 DO 50 J=K,1,-1
         R = ZERO
         DO 20 I=1,K
   20       IF (I .NE. J) R = R+ABS(MAT(J,I))
         IF (R .EQ. ZERO) THEN
            SCAL(K) = J
            IF (J .NE. K) THEN
               DO 30 I=1,K
   30             CALL SWAP(MAT(I,J),MAT(I,K))
               DO 40 I=L,N
   40             CALL SWAP(MAT(J,I),MAT(K,I))
            ENDIF
            K = K-1
            GOTO 10
         ENDIF
   50    CONTINUE
C
C     search for columns isolating an eigenvalue and
C     move them to the left
C
   60 DO 100 J=L,K
         C = ZERO
         DO 70 I=L,K
   70       IF (I .NE. J) C = C+ABS(MAT(I,J))
         IF (C .EQ. ZERO) THEN
            SCAL(L) = J
            IF (J .NE. L) THEN
               DO 80 I=1,K
   80             CALL SWAP(MAT(I,J),MAT(I,L))
               DO 90 I=L,N
   90             CALL SWAP(MAT(J,I),MAT(L,I))
            ENDIF
            L = L+1
            GOTO 60
         ENDIF
  100    CONTINUE
C
C     balance the matrix in rows L to K
C
      LOW = L
      HIGH = K
      DO 110 I=L,K
  110    SCAL(I) = ONE
  120 DO 180 I=L,K
         C = ZERO
         R = ZERO
         DO 130 J=L,K
            IF (J .NE. I) THEN
               C = C+ABS(MAT(J,I))
               R = R+ABS(MAT(I,J))
            ENDIF
  130       CONTINUE
         G = R/BASIS
         F = ONE
         S = C+R
  140    IF (C .LT. G) THEN
            F = F*BASIS
            C = C*B2
            GOTO 140
         ENDIF
         G = R*BASIS
  150    IF (C .GE. G) THEN
            F = F/BASIS
            C = C/B2
            GOTO 150
         ENDIF
         IF ((C+R)/F .LT. PT95*S) THEN
            G = ONE/F
            SCAL(I) = SCAL(I)*F
            DO 160 J=L,N
  160          MAT(I,J) = MAT(I,J)*G
            DO 170 J=1,K
  170          MAT(J,I) = MAT(J,I)*F
            GOTO 120
         ENDIF
  180    CONTINUE
      BALAN = 0
      END
C
C
C[BA*)
C[LE*)
      INTEGER FUNCTION BALBAK (LD,N,LOW,HIGH,SCAL,EIVEC)
C[IX{BALBAK}*)
C
C*****************************************************************
C                                                                *
C     BALBAK transforms all right eigenvectors of a balanced     *
C     matrix back into the eigenvectors of the original matrix.  *
C     The balanced matrix was generated by calling the procedure *
C     BALAN.                                                     *
C[BE*)
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     LD:       leading dimension of the array EIVEC as defined  *
C               in the calling program                           *
C     N:        order of the eigenvectors (number of components) *
C     LOW,HIGH: two INTEGER numbers that were created in the     *
C               procedure BALAN                                  *
C     SCAL:     output vector of the procedure BALAN             *
C     EIVEC:    2-dimensional array EIVEC(1:N,1:N), each column  *
C               of EIVEC contains an eigenvector (or its real    *
C               or imaginary parts) for the balanced matrix      *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     EIVEC:    the eigenvectors (or their real or imaginary     *
C               parts) of the original matrix                    *
C                                                                *
C     RETURN VALUE of SUBROUTINE BALBAK:                         *
C     ==================================                         *
C     0:        no error                                         *
C                                                                *
C     LOCAL VARIABLES:                                           *
C     ================                                           *
C     I,J,K:    auxiliary variables used for indexing            *
C     S:        scaling constant                                 *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: SWAP                                    *
C                                                                *
C                                                                *
C  sources : Parlett, B. N. and Reinsch, C., see [PARL69].       *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Juergen Dietel                                     *
C  date     : 04.10.1987                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      INTEGER          LD,N,LOW,HIGH
      DOUBLE PRECISION SCAL(1:N),EIVEC(1:LD,1:N)
      INTEGER          I,J,K
      DOUBLE PRECISION S
C
      DO 20 I=LOW,HIGH
         S = SCAL(I)
C
C        left eigenvectors are back transformed, if the last
C        statement is replaced by:  S = 1.0D0/SCAL(I)
C
         DO 10 J=1,N
   10       EIVEC(I,J) = EIVEC(I,J)*S
   20    CONTINUE
      DO 40 I=LOW-1,1,-1
         K=SCAL(I)
         DO 30 J=1,N
   30       CALL SWAP(EIVEC(I,J),EIVEC(K,J))
   40    CONTINUE
      DO 60 I=HIGH+1,N
         K=SCAL(I)
         DO 50 J=1,N
   50       CALL SWAP(EIVEC(I,J),EIVEC(K,J))
   60    CONTINUE
      BALBAK = 0
      END
C
C
C[BA*)
C[LE*)
      INTEGER FUNCTION ELMHES (LD,N,LOW,HIGH,MAT,PERM)
C[IX{ELMHES}*)
C
C*****************************************************************
C     Reduction of a general square matrix a(1:N,1:N) to         *
C     Hessenberg form by similarity with non-orthogonal Gaussian *
C     elimination matrices.                                      *
C*****************************************************************
C[BE*)
C                                                                *
C     PROGRAM OBJECTIVE:                                         *
C     ==================                                         *
C     For a general square matrix A(1:N,1:N) this program reduces*
C     the principal submatrix of A of order HIGH-LOW+1, defined  *
C     to lie between the elements A(LOW,LOW) and A(HIGH,HIGH),   *
C     to Hessenberg form H by similarity with non-orthogonal     *
C     Gaussian elimination matrices. The principal submatrix is  *
C     overwriten with H. The transforming matrices are stored in *
C     the triangle below H and in the vector PERM.               *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     LD:       leading dimension of the matrix as defined in    *
C               the calling program                              *
C     N:        order of the square matrix A                     *
C     LOW,HIGH: output parameter of the balancing procedure      *
C               BALAN. If A has not been balanced, set LOW:=1,   *
C               and HIGH:=N.                                     *
C     MAT:      2-dimensional array MAT(1:N,1:N), the matrix in  *
C               balanced form                                    *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     MAT:      2-dimensional array MAT(1:N,1:n), which consists *
C               in part of the upper Hessenberg matrix and the   *
C               transforming matrices.                           *
C               The number N(I,R+1), that is needed for the      *
C               reduction, is stored in the (I,R) position.      *
C     PERM:     INTEGER vector that stores the row- and column   *
C               permutations performed during the reduction      *
C                                                                *
C     RETURN VALUE:                                              *
C     =============                                              *
C     0:        no error                                         *
C                                                                *
C     LOCAL VARIABLES:                                           *
C     ================                                           *
C     ZERO,ONE: floating-point constants                         *
C     I,J,M:    counters                                         *
C     X,Y:      auxiliary variables used for storing matrix      *
C               elements and intermediate results                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: SWAP                                    *
C                                                                *
C                                                                *
C  sources : Martin, R. S. and Wilkinson, J. H., see [MART70].   *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Juergen Dietel                                     *
C  date     : 04.10.1987                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C*****************************************************************
C
      INTEGER          N,LOW,HIGH,PERM(1:N)
      DOUBLE PRECISION MAT(1:LD,1:N)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER        (ZERO = 0.0D0,ONE = 1.0D0)
      INTEGER          I,J,M
      DOUBLE PRECISION X,Y
C
      DO 70 M=LOW+1,HIGH-1
         I = M
         X = ZERO
         DO 10 J=M,HIGH
            IF (ABS(MAT(J,M-1)) .GT. ABS(X)) THEN
               X = MAT(J,M-1)
               I=J
            ENDIF
   10       CONTINUE
         PERM(M) = I
         IF (I .NE. M) THEN
C
C           Swap rows and columns of MAT
C
            DO 20 J=M-1,N
   20          CALL SWAP(MAT(I,J),MAT(M,J))
            DO 30 J=1,HIGH
   30          CALL SWAP(MAT(J,I),MAT(J,M))
         ENDIF
         IF (X .NE. ZERO) THEN
            DO 60 I=M+1,HIGH
               Y = MAT(I,M-1)
               IF (Y .NE. ZERO) THEN
                  Y = Y/X
                  MAT(I,M-1) = Y
                  DO 40 J=M,N
   40                MAT(I,J) = MAT(I,J)-Y*MAT(M,J)
                  DO 50 J=1,HIGH
   50                MAT(J,M) = MAT(J,M)+Y*MAT(J,I)
               ENDIF
   60          CONTINUE
         ENDIF
   70    CONTINUE
      ELMHES = 0
      END
C
C
C[BA*)
C[LE*)
      INTEGER FUNCTION ELMTRA (LD,N,LOW,HIGH,MAT,PERM,H)
C[IX{ELMTRA}*)
C
C*****************************************************************
C                                                                *
C     Form the matrix of accumulated transformations from the    *
C     information left by procedure ELMHES in the lower triangle *
C     of the Hessenberg matrix H - in MAT(1:N,1:N) and in the    *
C     INTEGER vector PERM(1:N). Store in the array H(1:N,1:N).   *
C[BE*)
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     LD:       leading dimension of the matrix as defined in    *
C               the calling program                              *
C     N:        order of the Hessenberg matrix H                 *
C     LOW,HIGH: INTEGER numbers, that were produced by procedure *
C               BALAN (if it was used; otherwise set LOW:=1,     *
C               HIGH:=N.)                                        *
C     PERM:     INTEGER N-vector produced by ELMHES              *
C     MAT:      (N,N)-matrix, that was produced by ELMHES and    *
C               contains the Hessenberg matrix H and the LR      *
C               multipliers used in order to create H from the   *
C               given matrix A                                   *
C                                                                *
C     OUTPUT PARAMETER:                                          *
C     =================                                          *
C     H:        (N,N)-matrix that describes the similarity       *
C               transformation of A to Hessenberg form H         *
C                                                                *
C     RETURN VALUE:                                              *
C     =============                                              *
C     0:             no error                                    *
C                                                                *
C     LOCAL FACTORS:                                             *
C     ==============                                             *
C     ZERO,ONE: floating-point constants                         *
C     I,J,K:    index variables                                  *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C                                                                *
C  sources : Peters, G. and Wilkinson, J. H., see [PETE70].      *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Juergen Dietel                                     *
C  date     : 04.10.1987                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      INTEGER          LD,N,LOW,HIGH,PERM(1:N)
      DOUBLE PRECISION MAT(1:LD,1:N),H(1:LD,1:N)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER        (ZERO = 0.0D0,ONE = 1.0D0)
      INTEGER          I,J,K
C
      DO 20 I=1,N
         DO 10 J=1,N
   10       H(I,J) = ZERO
         H(I,I) = ONE
   20    CONTINUE
C
      DO 50 I=HIGH-1,LOW+1,-1
         J=PERM(I)
         DO 30 K=I+1,HIGH
   30       H(K,I) = MAT(K,I-1)
         IF (I .NE. J) THEN
            DO 40 K=I,HIGH
               H(I,K)=H(J,K)
               H(J,K)=ZERO
   40          CONTINUE
            H(J,I) = ONE
         ENDIF
   50    CONTINUE
      ELMTRA = 0
      END
C
C
C[BA*)
C[LE*)
      INTEGER FUNCTION ORTHES (LD,N,LOW,HIGH,MAT,D,EPSM)
C[IX{ORTHES}*)
C
C*****************************************************************
C                                                                *
C  This function transform the matrix MAT to upper Hessenberg    *
C  using Householder transforms.                                 *
C[BE*)
C  The essential transform information is stored in the otherwise*
C  unused lower triangle of MAT and in the vector D.             *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C LD       leading dimension of the matrix MAT, as defined in the*
C          calling program                                       *
C N        order of the matrix MAT                               *
C LOW  \   the rows 1 to LOW-1 and HIGH+1 to N contain the       *
C HIGH  >  isolated eigenvalues, i. e., those eigenvalues that   *
C      /   have unit vectors e_i as eigenvectors.                *
C MAT      the original matrix MAT(1:N,1:N)                      *
C EPSM     machine constant                                      *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C MAT      the desired Hessenberg matrix with part of the        *
C          transform information below the subdiagonal           *
C D        N-vector with the remainder of the transform info     *
C                                                                *
C RETURN VALUE :                                                 *
C ==============                                                 *
C Error code. Only  0 (no error) is possible here.               *
C                                                                *
C LOCAL VARIABLES:                                               *
C ================                                               *
C I,J,M  Loop variables                                          *
C S      Euclidean norm sigma of a column below the subdiagonal  *
C        v of MAT, which must be transformed into a multiple of  *
C        e1 = (1,0,...,0); (v = (v1,...,v(HIGH-M+1))             *
C X      initially leading element of v, then summation value    *
C        inside the Householder transformation                   *
C Y      initially  sigma^2, then ||u||^2, where                 *
C        u := v +- sigma * e1                                    *
C EPS    accuracy bound to check transformation                  *
C ZERO   double precision zero 0.0D0                             *
C                                                                *
C*****************************************************************
C Reference: R.S. Martin and J.H. Wilkinson, Num. Math. 12 (1968)*
C            pp. 359, 360, see [WILK71], contrib. II/13          *
C Author:    Juergen Dietel, Computer Center, RWTH Aachen        *
C Date:      7. 15. 1993                                         *
C[BA*)
C*****************************************************************
C[BE*)
C
      INTEGER          LD,N,LOW,HIGH
      DOUBLE PRECISION MAT(1:LD,1:N),D(1:N),EPSM
C
      INTEGER          I,J,M
      DOUBLE PRECISION S
      DOUBLE PRECISION X
      DOUBLE PRECISION Y
      DOUBLE PRECISION EPS
C
      DOUBLE PRECISION ZERO
      PARAMETER        (ZERO = 0.0D0)
C
      EPS = 128 * EPSM
C
      DO 80 M = LOW+1,HIGH-1
C
         Y = ZERO
         DO 10 I=HIGH,M,-1
            X    = MAT(I,M-1)
            D(I) = X
            Y    = Y + X * X
   10       CONTINUE
         IF (Y .LE. EPS) THEN
            S = ZERO
         ELSE
C
            IF (X .GE. ZERO) THEN
               S = -SQRT(Y)
            ELSE
               S = SQRT(Y)
            ENDIF
            Y    = Y - X * S
            D(M) = X - S
C
C           Multiply MAT on the left by  (E-(u * uT)/y)
C
            DO 40 J=M,N
               X = ZERO
               DO 20 I=HIGH,M,-1
                  X = X + D(I) * MAT(I,J)
   20             CONTINUE
               X = X / Y
               DO 30 I=M,HIGH
                  MAT(I,J) = MAT(I,J) - X * D(I)
   30             CONTINUE
   40          CONTINUE
C
C           Multiply MAT on the right by (E-(u * uT)/y)
C
            DO 70 I=1,HIGH
               X = ZERO
               DO 50 J=HIGH,M,-1
                  X = X + D(J) * MAT(I,J)
   50             CONTINUE
               X = X / Y
               DO 60 J=M,HIGH
                  MAT(I,J) = MAT(I,J) - X * D(J)
   60             CONTINUE
   70          CONTINUE
C
         ENDIF
C
        MAT(M,M-1) = S
   80   CONTINUE
C
      ORTHES = 0
C
      END
C
C
C[BA*)
C[LE*)
      INTEGER FUNCTION ORTTRA (LD,N,LOW,HIGH,MAT,D,V)
C[IX{ORTTRA}*)
C[BE*)
C
C*****************************************************************
C                                                                *
C  Reconstruct the transformation matrix V of the Householder    *
C  reductions of MAT to Hessenberg form from the information     *
C  stored in the lower triangle of MAT and in D.                 *
C  The contents of D is lost.                                    *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C LD       leading dimension of the matrices  MAT and V, as      *
C          defined in the calling routine                        *
C N        size of the matrix MAT                                *
C LOW  \   the rows  1 to LOW-1 and the rows  HIGH+1  to N       *
C HIGH  >  contain the isolated eigenvalues, i.e. those with     *
C      /   unit vectors e_i as eigenvectors.                     *
C MAT      (1:N,1:N) matrix, that ORTHES has reduced to upper    *
C          Hessenberg form; partially filled with transformation *
C          information                                           *
C D        N-vector with the remaining transform information     *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C D        modified input vector                                 *
C V        (1:N,1:N)-matrix, giving the similarity transform     *
C          from A to the upper Hessenberg matrix in MAT          *
C                                                                *
C RETURN VALUE :                                                 *
C ==============                                                 *
C Error code. Only 0 possible here (no error).                   *
C                                                                *
C LOCAL VARIABLES:                                               *
C ================                                               *
C I,J,M  loop variables                                          *
C X      summation variable during Householder transformation    *
C Y      sigma  or  sigma * (v1 +- sigma)                        *
C ZERO   double precision  0.0D0                                 *
C ONE    double precision  1.0D0                                 *
C                                                                *
C*****************************************************************
C Reference: G. Peters and J.H. Wilkinson, Num. Math. 16 (1970)  *
C            p. 191, see [WILK71], contrib. II/15                *
C Author:    Juergen Dietel, Computer Center, RWTH Aachen        *
C Date:      7. 15. 1993                                         *
C*****************************************************************
C
      INTEGER          LD,N,LOW,HIGH
      DOUBLE PRECISION MAT(1:LD,1:N),D(1:N),V(1:LD,1:N)
C
      INTEGER          I,J,M
      DOUBLE PRECISION X
      DOUBLE PRECISION Y
C
      DOUBLE PRECISION ZERO,ONE
      PARAMETER        (ZERO = 0.0D0,ONE = 1.0D0)
C
C     start with identity matrix in V
C
      DO 20 I=1,N
         DO 10 J=1,N
           V(I,J) = ZERO
   10      CONTINUE
         V(I,I) = ONE
   20    CONTINUE
C
C     Perform the transformations that have reduced MAT to
C     Hessenberg form on the unit matrix stored in V in order
C     to generate the transformation matrix.
C
      DO 70 M=HIGH-1,LOW+1,-1
         Y = MAT(M,M-1)
C
         IF (Y .NE. ZERO) THEN
            Y = Y * D(M)
            DO 30 I=M+1,HIGH
               D(I) = MAT(I,M-1)
   30          CONTINUE
            DO 60 J=M,HIGH
            X = ZERO
            DO 40 I=M,HIGH
               X = X + D(I) * V(I,J)
   40          CONTINUE
            X = X / Y
            DO 50 I=M,HIGH
               V(I,J) = V(I,J) + X * D(I)
   50          CONTINUE
   60       CONTINUE
         ENDIF
   70    CONTINUE
C
      ORTTRA = 0
C
      END
C
C
C[BA*)
C[LE*)
      INTEGER FUNCTION HQR2 (VEC,LD,N,LOW,HIGH,H,VALR,VALI,
C[IX{HQR2}*)
     *                       EIVEC,CNT,EPS)
C
C*****************************************************************
C     Finds the eigenvalues and eigenvectors of a real matrix,   *
C     which has been reduced to upper Hessenberg form.           *
C*****************************************************************
C[BE*)
C                                                                *
C     Finds the eigenvalues and eigenvectors (if VEC = .TRUE.) of*
C     a real matrix, which has been reduced to upper Hessenberg  *
C     form and is stored in the array H(1:N,1:N) with the accu-  *
C     mulated transformations stored in the array EIVEC(1:N,1:N).*
C     The real and imaginary parts of the eigenvalues are placed *
C     in the two vectors VALR(1:N), VALI(1:N) while the eigen-   *
C     vectors are stored in the array EIVEC(1:N,1:N), where only *
C     one complex eigenvector corresponding to the eigenvalue    *
C     with positive imaginary part is stored for a complex con-  *
C     jugate eigenvalue pair. LOW and HIGH are two INTEGER       *
C     numbers produced during balancing MAT, so that eigenvalues *
C     are isolated in positions 1 to LOW-1 and HIGH+1 to N. If   *
C     no initial balancing was performed, set LOW:=1, HIGH:=N.   *
C     The subroutine is aborted with an error message if any one *
C     of the eigenvalues requires more than MAXSTP iterations.   *
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     VEC:      parameter indicating eigenvector computation:    *
C                 =  .TRUE.  : compute eigenvectors              *
C                 =  .FALSE. : compute eigenvalues only          *
C     LD:       leading dimension of H and EIVEC, as defined in  *
C               the calling program                              *
C     N:        order of the Hessenberg matrix H                 *
C     LOW,HIGH: INTEGER numbers produced by BALAN, if it was     *
C               used. Otherwise set LOW:=1, HIGH:=N.             *
C     EPS:      machine constant                                 *
C     H:        (N,N) matrix containing H with its relevant      *
C               components                                       *
C     EIVEC:    (N,N) matrix containing the array that describes *
C               the similarity transformation of A to H..        *
C               (This was procuded by ELMTRA or ORTTRA.)         *
C               If H is the original matrix, define EIVEC := I,  *
C               the identity matrix.                             *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     H:           the upper triangle of this (N,N) matrix       *
C                  contains the eigenvectors of the quasi-trian- *
C                  gular matrix, that is produced by the QR steps*
C     VALR,VALI:   two N-vectors, with the real and imaginary    *
C                  parts of the eigenvalues                      *
C     CNT:         INTEGER N-vector, with the count of iterations*
C                  for each eigenvalue. If two eigenvalues are   *
C                  found as a complex conjugate pair, the number *
C                  of iterations is entered with a positive sign *
C                  for the first and with a negative sign for the*
C                  second eigenvalue.                            *
C     EIVEC:       (N,N) matrix, where (for VEC = .TRUE.) the    *
C                  non-normalized eigenvectors of the original   *
C                  matrix are stored in case H is not the        *
C                  original matrix.                              *
C                  If the I-th eigenvalue is real, the I-th      *
C                  column of EIVEC contains the corresponding    *
C                  real eigenvector. If the eigenvalues numbered *
C                  I and I+1 form a complex conjugate eigenvalue *
C                  pair for MAT, then the I-th and (I+1)-st      *
C                  columns of EIVEC contain the real and         *
C                  imaginary parts of the eigenvector for the    *
C                  corresponding eigenvalue with positive real   *
C                  part.                                         *
C                                                                *
C                                                                *
C     RETURN VALUES:                                             *
C     ==============                                             *
C     0:           no error                                      *
C     1:           parameters N, LOW or HIGH are unacceptable    *
C     2:           all eigenvectors are the zero vector.         *
C     3:           the maximum number of QR-steps has been       *
C                  exceeded.                                     *
C                                                                *
C     LOCAL VARIABLES:                                           *
C     ================                                           *
C     ZERO,ONE,TWO,PT75,PT4375: floating-point constants         *
C     MAXSTP:                   maximum number of allowed steps  *
C     I,J,K,L,M,N,NA,EN:        index variables                  *
C     ITER:                     step counter                     *
C     P,Q,R,S,T,W,X,Y,Z,NORM,                                    *
C       RA,SA,VR,VI:            auxiliary variables              *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: COMDIV                                  *
C                                                                *
C                                                                *
C  sources : Peters, G. and Wilkinson, J. H., see [PETE70].      *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Juergen Dietel                                     *
C  date     : 07.14.1993                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C*****************************************************************
C
      LOGICAL          VEC
      INTEGER          LD,N,LOW,HIGH,CNT(1:N)
      DOUBLE PRECISION H(1:LD,1:N),EIVEC(1:LD,1:N),VALR(1:N)
      DOUBLE PRECISION VALI(1:N),EPS
      DOUBLE PRECISION ZERO,ONE,TWO,PT75,PT4375
      INTEGER          MAXSTP
      PARAMETER        (ZERO = 0.0D0,ONE = 1.0D0,TWO = 2.0D0,
     *                 PT75 = 0.75D0,PT4375 = 0.4375D0,
     *                 MAXSTP = 100)
      INTEGER          I,J,K,L,M,NA,EN,ITER
      DOUBLE PRECISION P,Q,R,S,T,W,X,Y,Z,NORM,RA,SA,VR,VI
C
C     error 1: one of the parameters N, LOW or HIGH has an
C              unacceptable value
C
      IF (N .LT. 1 .OR. LOW .LT. 1 .OR. HIGH .GT. N) THEN
         HQR2 = 1
         RETURN
      ENDIF
C
C     initialize the isolated eigenvalues found during balancing:
C
      DO 10 I=1,N
         IF (I .LT. LOW .OR. I .GT. HIGH) THEN
            VALR(I) = H(I,I)
            VALI(I) = ZERO
            CNT(I) = 0
         ELSE
            CNT(I) = -1
         ENDIF
   10    CONTINUE
C
      EN = HIGH
      T = ZERO
   15    IF (EN .LT. LOW) GOTO 333
         ITER = 0
         NA = EN-1
C
C           search for a single small subdiagonal element:
C
   20       DO 30 L=EN,LOW+1,-1
   30          IF (ABS(H(L,L-1)) .LE. EPS*
     *             (ABS(H(L-1,L-1))+ABS(H(L,L)))) GOTO 40
   40       X = H(EN,EN)
            IF (L .EQ. EN) THEN
C
C              found one root:
C
               VALR(EN) = X + T
               H(EN,EN) = VALR(EN)
               VALI(EN) = ZERO
               CNT(EN) = ITER
               EN = NA
               GOTO 15
            ENDIF
C
            Y = H(NA,NA)
            W = H(EN,NA)*H(NA,EN)
            IF (L .EQ. NA) THEN
C
C              found two roots:
C
               P = (Y-X)/TWO
               Q = P*P+W
               Z = SQRT(ABS(Q))
               H(EN,EN) = X+T
               X = H(EN,EN)
               H(NA,NA) = Y+T
               CNT(EN) = -ITER
               CNT(NA) = ITER
               IF (Q .GE. ZERO) THEN
C
C                 found a real pair:
C
                  IF (P .LT. ZERO) Z = -Z
                  Z = P+Z
                  VALR(NA) = X+Z
                  VALR(EN) = X-W/Z
                  VALI(NA) = ZERO
                  VALI(EN) = ZERO
                  X = H(EN,NA)
                  R = SQRT(X*X+Z*Z)

                  IF (VEC) THEN
                     P = X/R
                     Q = Z/R
C
C                    row modification:
C
                     DO 50 J=NA,N
                        Z = H(NA,J)
                        H(NA,J) = Q*Z+P*H(EN,J)
                        H(EN,J) = Q*H(EN,J)-P*Z
   50                   CONTINUE
C
C                    column modification:
C
                     DO 60 I=1,EN
                        Z = H(I,NA)
                        H(I,NA) = Q*Z+P*H(I,EN)
                        H(I,EN) = Q*H(I,EN)-P*Z
   60                   CONTINUE
C
C                    accumulate transformations:
C
                     DO 70 I=LOW,HIGH
                        Z = EIVEC(I,NA)
                        EIVEC(I,NA) = Q*Z+P*EIVEC(I,EN)
                        EIVEC(I,EN) = Q*EIVEC(I,EN)-P*Z
   70                   CONTINUE
                  ENDIF
               ELSE
C
C                 complex conjugate pair:
C
                  VALR(NA) = X+P
                  VALR(EN) = VALR(NA)
                  VALI(NA) = Z
                  VALI(EN) = -Z
               ENDIF
               EN = EN-2
               GOTO 15
            ENDIF
C
            IF (ITER .EQ. MAXSTP) THEN
C
C           error 3: maximal number of iterations exceeded:
C
               CNT(EN) = MAXSTP+1
               HQR2 = 3
               RETURN
            ENDIF
            IF (MOD(ITER,10) .EQ. 0 .AND. ITER .NE. 0) THEN
C
C              use exceptional shift:
C
               T = T+X
               DO 80 I=LOW,EN
   80             H(I,I) = H(I,I)-X
               S = ABS(H(EN,NA))+ABS(H(NA,EN-2))
               X = PT75*S
               Y = X
               W = -PT4375*S*S
            ENDIF
            ITER = ITER+1
C
C        search for two consecutive small subdiagonal elements:
C
            DO 90 M=EN-2,L,-1
               Z = H(M,M)
               R = X-Z
               S = Y-Z
               P = (R*S-W)/H(M+1,M)+H(M,M+1)
               Q = H(M+1,M+1)-Z-R-S
               R = H(M+2,M+1)
               S = ABS(P)+ABS(Q)+ABS(R)
               P = P/S
               Q = Q/S
               R = R/S
               IF (M .EQ. L) GOTO 100
               IF (ABS(H(M,M-1))*(ABS(Q)+ABS(R)) .LE. EPS*ABS(P)*
     *             (ABS(H(M-1,M-1))+ABS(Z)+ABS(H(M+1,M+1))))
     *             GOTO 100
   90          CONTINUE
  100       DO 110 I=M+2,EN
  110          H(I,I-2) = ZERO
            DO 120 I=M+3,EN
  120          H(I,I-3) = ZERO
C
C           double QR-step involving rows L to EN and
C           columns M to EN of the complete array:
C
            DO 200 K=M,NA
               IF (K .NE. M) THEN
                  P = H(K,K-1)
                  Q = H(K+1,K-1)
                  IF (K .NE. NA) THEN
                     R = H(K+2,K-1)
                  ELSE
                     R = ZERO
                  ENDIF
                  X = ABS(P)+ABS(Q)+ABS(R)
                  IF (X .EQ. ZERO) GOTO 200
                  P = P/X
                  Q = Q/X
                  R = R/X
               ENDIF
               S = SQRT(P*P+Q*Q+R*R)
               IF (P .LT. ZERO) S = -S
               IF (K .NE. M) THEN
                  H(K,K-1) = -S*X
               ELSEIF (L .NE. M) THEN
                  H(K,K-1) = -H(K,K-1)
               ENDIF
               P = P+S
               X = P/S
               Y = Q/S
               Z = R/S
               Q = Q/P
               R = R/P
C
C              row modification:
C
               DO 130 J=K,N
                  P = H(K,J)+Q*H(K+1,J)
                  IF (K .NE. NA) THEN
                     P = P+R*H(K+2,J)
                     H(K+2,J) = H(K+2,J)-P*Z
                  ENDIF
                  H(K+1,J) = H(K+1,J)-P*Y
                  H(K,J) = H(K,J)-P*X
  130             CONTINUE
               J = MIN(K+3,EN)
C
C              column modification:
C
               DO 140 I=1,J
                  P = X*H(I,K)+Y*H(I,K+1)
                  IF (K .NE. NA) THEN
                     P = P+Z*H(I,K+2)
                     H(I,K+2) = H(I,K+2)-P*R
                  ENDIF
                  H(I,K+1) = H(I,K+1)-P*Q
                  H(I,K) = H(I,K)-P
  140             CONTINUE
C
               IF (VEC) THEN
C
C                 accumulate transformations:
C
                  DO 150 I=LOW,HIGH
                     P = X*EIVEC(I,K)+Y*EIVEC(I,K+1)
                     IF (K .NE. NA) THEN
                        P = P+Z*EIVEC(I,K+2)
                        EIVEC(I,K+2) = EIVEC(I,K+2)-P*R
                     ENDIF
                     EIVEC(I,K+1) = EIVEC(I,K+1)-P*Q
                     EIVEC(I,K) = EIVEC(I,K)-P
  150                CONTINUE
               ENDIF
  200          CONTINUE
            GOTO 20
C
C
  333 IF (.NOT. VEC) THEN
         HQR2 = 0
         RETURN
      ENDIF
C
C
C     all eigenvalues have been found; now transform back:
C
C     find the 1-norm of H :
C
      NORM = ZERO
      K=1
      DO 201 I=1,N
         DO 101 J=K,N
  101       NORM = NORM+ABS(H(I,J))
  201    K = I
      IF (NORM .EQ. ZERO) THEN
C        Fehler 2: 1-Norm von H ist gleich 0:
         HQR2 = 2
         RETURN
      ENDIF
C
C     back transformation:
C
      DO 207 EN=N,1,-1
         P = VALR(EN)
         Q = VALI(EN)
         NA = EN - 1
         IF (Q .EQ. ZERO) THEN
C
C           real vector:
C
            M = EN
            H(EN,EN) = ONE
            DO 63 I=NA,1,-1
               W = H(I,I)-P
               R = H(I,EN)
               DO 38 J=M,NA
   38             R = R+H(I,J)*H(J,EN)
               IF (VALI(I) .LT. ZERO) THEN
                  Z = W
                  S = R
               ELSE
                  M = I
                  IF (VALI(I) .EQ. ZERO) THEN
                     IF (W .NE. ZERO) THEN
                        H(I,EN) = -R/W
                     ELSE
                        H(I,EN) = -R/(EPS*NORM)
                     ENDIF
                  ELSE
C
C                    solve the linear system:
C                    [ W   X ] [ H(I,EN)   ]   [ -R ]
C                    [       ] [           ] = [    ]
C                    [ Y   Z ] [ H(I+1,EN) ]   [ -S ]
C
                     X = H(I,I+1)
                     Y = H(I+1,I)
                     Q = (VALR(I)-P)*(VALR(I) - P)+
     *                   VALI(I)*VALI(I)
                     T = (X*S-Z*R)/Q
                     H(I,EN) = T
                     IF (ABS(X) .GT. ABS(Z)) THEN
                        H(I+1,EN) = (-R-W*T)/X
                     ELSE
                        H(I+1,EN) = (-S-Y*T)/Z
                     ENDIF
                  ENDIF
               ENDIF
   63          CONTINUE
         ELSEIF (Q .LT. ZERO) THEN
C
C           complexer eigenvector for LAMBDA = P - I * Q :
C
            M = NA
            IF (ABS(H(EN,NA)) .GT. ABS(H(NA,EN))) THEN
               H(NA,NA) = -(H(EN,EN)-P)/H(EN,NA)
               H(NA,EN) = -Q/H(EN,NA)
            ELSE
               CALL COMDIV(-H(NA,EN),ZERO,H(NA,NA)-P,Q,
     *                     H(NA,NA),H(NA,EN))
            ENDIF
            H(EN,NA) = ONE
            H(EN,EN) = ZERO
            DO 190 I=NA-1,1,-1
               W = H(I,I)-P
               RA = H(I,EN)
               SA = ZERO
               DO 75 J=M,NA
                  RA = RA+H(I,J)*H(J,NA)
                  SA = SA+H(I,J)*H(J,EN)
   75             CONTINUE
               IF (VALI(I) .LT. ZERO) THEN
                  Z = W
                  R = RA
                  S = SA
               ELSE
                  M = I
                  IF (VALI(I) .EQ. ZERO) THEN
                      CALL COMDIV(-RA,-SA,W,Q,H(I,NA),H(I,EN))
                  ELSE
C
C                    solve the complex linear system:
C           [ W+Q*I   X   ] [H(I,NA)+H(I,EN)*I    ]   [-RA-SA*I]
C           [             ] [                     ] = [        ]
C           [   Y   Z+Q*I ] [H(I+1,NA)+H(I+1,EN)*I]   [-R-S*I  ]
C
                     X = H(I,I+1)
                     Y = H(I+1,I)
                     VR = (VALR(I)-P)*(VALR(I)-P)+
     *                    VALI(I)*VALI(I)-Q*Q
                     VI = TWO*Q*(VALR(I)-P)
                     IF (VR .EQ. ZERO .AND. VI .EQ. ZERO) VR =
     *                  EPS*NORM*
     *                  (ABS(W)+ABS(Q)+ABS(X)+ABS(Y)+ABS(Z))
                     CALL COMDIV(X*R-Z*RA+Q*SA,X*S-Z*SA-Q*RA,
     *                           VR,VI,H(I,NA),H(I,EN))
                     IF (ABS(X) .GT. ABS(Z)+ABS(Q)) THEN
                        H(I+1,NA) = (-RA-W*H(I,NA)+Q*H(I,EN))/X
                        H(I+1,EN) = (-SA-W*H(I,EN)-Q*H(I,NA))/X
                     ELSE
                       CALL COMDIV(-R-Y*H(I,NA),-S-Y*H(I,EN),Z,Q,
     *                             H(I+1,NA),H(I+1,EN))
                     ENDIF
                  ENDIF
               ENDIF
  190          CONTINUE
         ENDIF
  207    CONTINUE
C
C     find eigenvectors for isolated eigenvalues:
C
      DO 230 I=1,N
         IF (I .LT. LOW .OR. I .GT. HIGH) THEN
            DO 220 J=I+1,N
  220          EIVEC(I,J) = H(I,J)
         ENDIF
  230    CONTINUE
C
C     multiply by transformation matrix in order to
C     obtain eigenvectors of the original matrix MAT:
C
      DO 300 J=N,LOW,-1
         IF (J .LE. HIGH) THEN
            M = J
         ELSE
            M = HIGH
         ENDIF
         L = J-1
         IF (VALI(J) .LT. ZERO) THEN
            DO 330 I=LOW,HIGH
               Y = ZERO
               Z = ZERO
               DO 320 K=LOW,M
                  Y = Y+EIVEC(I,K)*H(K,L)
                  Z = Z+EIVEC(I,K)*H(K,J)
  320             CONTINUE
               EIVEC(I,L) = Y
               EIVEC(I,J) = Z
  330          CONTINUE
         ELSE
            IF (VALI(J) .EQ. ZERO) THEN
               DO 350 I=LOW,HIGH
                  Z = ZERO
                  DO 340 K=LOW,M
  340                Z =Z+EIVEC(I,K)*H(K,J)
  350             EIVEC(I,J) = Z
            ENDIF
         ENDIF
  300    CONTINUE
C
C     Return  0: no error
C
      HQR2 = 0
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE COMDIV (AR,AI,BR,BI,RESR,RESI)
C[IX{COMDIV}*)
C
C*****************************************************************
C                                                                *
C     complex division: RESR+I*RESI := (AR+I*AI)/(BR+I*BI).      *
C     (this procedure should not be called if BR=BI=0.)          *
C[BE*)
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     AR,AI:     real and imaginary parts of the numerator       *
C     BR,BI:     real and imaginary parts of the denominator     *
C                                                                *
C     OUTPUT PARAMETERS:                                         *
C     ==================                                         *
C     RESR,RESI: real and imaginary parts of the quotient        *
C                                                                *
C     LOCAL VARIABLES:                                           *
C     ================                                           *
C     ZERO:              floating-point constant 0               *
C     TEMP1,TEMP2,TEMP3: auxiliary variables for intermediate    *
C                        values                                  *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: none                                    *
C                                                                *
C                                                                *
C  sources : Martin, R. S. and Wilkinson, J. H., see [MART68].   *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Juergen Dietel                                     *
C  date     : 04.10.1987                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      DOUBLE PRECISION AR,AI,BR,BI,RESR,RESI
      DOUBLE PRECISION ZERO
      PARAMETER        (ZERO = 0.0D0)
      DOUBLE PRECISION TEMP1,TEMP2,TEMP3
C
      IF (BR .EQ. ZERO .AND. BI .EQ. ZERO) THEN
         RESR = ZERO
         RESI = ZERO
         RETURN
      ENDIF
      IF (ABS(BR) .GT. ABS(BI)) THEN
         TEMP1 = BI/BR
         TEMP2 = TEMP1*BI+BR
         TEMP3 = (AR+TEMP1*AI)/TEMP2
         RESI = (AI-TEMP1*AR)/TEMP2
         RESR = TEMP3
      ELSE
         TEMP1 = BR/BI
         TEMP2 = TEMP1*BR+BI
         TEMP3 = (TEMP1*AR+AI)/TEMP2
         RESI = (TEMP1*AI-AR)/TEMP2
         RESR = TEMP3
      ENDIF
      END
C
C
C[BA*)
C[LE*)
      DOUBLE PRECISION FUNCTION COMABS (AR,AI)
C[IX{COMABS}*)
C
C*****************************************************************
C                                                                *
C     Determine the absolute value of the complex number AR+I*AI:*
C     COMABS:=SQRT(AR*AR+AI*AI)                                  *
C[BE*)
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     AR,AI: the real and imaginary parts of the complex number  *
C            whose absolute value is to be determined            *
C                                                                *
C     OUTPUT PARAMETER:                                          *
C     =================                                          *
C     none                                                       *
C                                                                *
C     RETURN VALUE:                                              *
C     =============                                              *
C     value of the complex parameter                             *
C                                                                *
C     LOCAL VARIABLES:                                           *
C     ================                                           *
C     ZERO,ONE:    constants                                     *
C     TEMP1,TEMP2: auxiliary variables                           *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: SWAP                                    *
C                                                                *
C                                                                *
C  sources : Martin, R. S. and Wilkinson, J. H., see [MART68].   *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Juergen Dietel                                     *
C  date     : 04.10.1987                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      DOUBLE PRECISION AR,AI
      DOUBLE PRECISION ZERO,ONE
      PARAMETER        (ZERO = 0.0D0,ONE = 1.0D0)
      DOUBLE PRECISION TEMP1,TEMP2
C
      TEMP1 = ABS(AR)
      TEMP2 = ABS(AI)
      IF (AR .EQ. ZERO .AND. AI .EQ. ZERO) THEN
         COMABS = ZERO
         RETURN
      ENDIF
      IF (TEMP2 .GT. TEMP1) CALL SWAP(TEMP1,TEMP2)
      IF (TEMP2 .EQ. ZERO) THEN
         COMABS = TEMP1
      ELSE
         COMABS = TEMP1*SQRT(ONE+(TEMP2/TEMP1)**2)
      ENDIF
      END
C
C
C[BA*)
C[LE*)
      INTEGER FUNCTION NORMAL (LD,N,V,WI)
C[IX{NORMAL}*)
C
C*****************************************************************
C                                                                *
C     NORMAL normalizes a set of eigenvectors in the maximum norm*
C[BE*)
C                                                                *
C     INPUT PARAMETERS:                                          *
C     =================                                          *
C     LD:     leading dimension of array V as defined in the     *
C             calling program                                    *
C     N:      the order of the square array V                    *
C     V:      (N,N) matrix of type DOUBLE PRECISION that contains*
C             the eigenvectors column-wise (see its description  *
C             in EIGEN under EIVEC)                              *
C     WI:     N-vector WI(1:N) of type DOUBLE PRECISION. Its     *
C             components are the imaginary parts of the          *
C             corresponding eigenvalues                          *
C                                                                *
C     OUTPUT PARAMETER:                                          *
C     =================                                          *
C     V:      array containing the normalized eigenvectors       *
C                                                                *
C     LOCAL VARIABLES:                                           *
C     ================                                           *
C     ZERO,ONE: floating-point constants 0 and 1                 *
C     I,J:      indexing variables                               *
C     MAXI:     auxiliary variable for determining the norm of   *
C               a real vector                                    *
C     TR,TI:    auxiliary variables for determining the norm of  *
C               a complex vector                                 *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  subroutines required: COMABS, COMDIV                          *
C                                                                *
C*****************************************************************
C                                                                *
C  author   : Juergen Dietel                                     *
C  date     : 04.10.1987                                         *
C  source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
      INTEGER          LD,N
      DOUBLE PRECISION V(1:LD,1:N),WI(1:N)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER        (ZERO = 0.0D0,ONE = 1.0D0)
      INTEGER          I,J
      DOUBLE PRECISION MAXI,TR,TI,COMABS
C
      J = 1
   10 IF (J .GT. N) GOTO 80
         IF (WI(J) .EQ. ZERO) THEN
             MAXI = V(1,J)
             DO 15 I=2,N
   15           IF (ABS(V(I,J)) .GT. ABS(MAXI)) MAXI = V(I,J)
             IF (MAXI .NE. ZERO) THEN
                MAXI = ONE/MAXI
                DO 20 I=1,N
   20              V(I,J) = V(I,J)*MAXI
             ENDIF
             J = J+1
         ELSE
            TR = V(1,J)
            TI = V(1,J+1)
            DO 30 I=2,N
               IF (COMABS(V(I,J),V(I,J+1)) .GT. COMABS(TR,TI))
     *         THEN
                  TR = V(I,J)
                  TI = V(I,J+1)
               ENDIF
   30          CONTINUE
            IF (TR .NE. ZERO .OR. TI .NE. ZERO) THEN
               DO 40 I=1,N
   40             CALL COMDIV(V(I,J),V(I,J+1),TR,TI,
     *                        V(I,J),V(I,J+1))
            ENDIF
            J = J+2
         ENDIF
         GOTO 10
   80 NORMAL = 0
      END
C
C
C[BA*)
C[LE*)
      SUBROUTINE SWAP (X,Y)
C[IX{SWAP}*)
C
C*****************************************************************
C                                                                *
C     PROGRAM OBJECTIVE:                                         *
C     ==================                                         *
C     SWAP interchanges the values of the two DOUBLE PRECISION   *
C     variables X and Y.                                         *
C                                                                *
C*****************************************************************
C[BE*)
C
      DOUBLE PRECISION X,Y
      DOUBLE PRECISION TEMP
      TEMP = X
      X    = Y
      Y    = TEMP
      END
