![  {QR Algorithm}                                                      
![  {Eigenvalues and Eigenvectors of a Matrix via the                   
![   QR Algorithm}*)                                                    
      INTEGER FUNCTION EIGEN (VEC, ORTHO, EVNORM, BASIS, LD, N, MAT,    &
      SCAL, D, EIVEC, VALR, VALI, CNT, LOW, HIGH)                       
!                                                                       
!*****************************************************************      
!                                                                *      
!     PROGRAM OBJECTIVE:                                         *      
!     ==================                                         *      
!     This FUNCTION-subroutine of type INTEGER determines all    *      
!     eigenvalues and eigenvectors of a real square matrix MAT.  *      
!     The eigenvalues are stored in the vectors VALR(1:N)        *      
!     (real parts) and VALI(1:N) (imaginary parts).              *      
!     The eigenvectors are stored in the array EIVEC(1:N,1:N)    *      
!     provided the flag VEC is set.                              *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     VEC:         parameter designating eigenvector computation:*      
!                    =  .TRUE. : compute eigenvectors            *      
!                    =  .FALSE.: compute eigenvalues only        *      
!     ORTHO:       flag indicating desired reduction of MAT to   *      
!                  Hessenberg form : if .TRUE. use orthogonal    *      
!                  transformations in ORTHES; if .FALSE. use     *      
!                  elementary Gauss eliminations via ELMHES which*      
!                  is slightly faster. For symmetric matrices MAT*      
!                  only ORTHES will preserve symmetry.           *      
!     EVNORM:      flag that governs potential normalizing of    *      
!                  eigenvectors if .TRUE.                        *      
!     N:           order of the square matrix MAT                *      
!     LD:          leading dimension of the matrix MAT and the   *      
!                  vector EIVEC as defined in the calling program*      
!     MAT:         2-dimensional array MAT(1:N,1:N) of type      *      
!                  DOUBLE PRECISION, containing the real input   *      
!                  matrix                                        *      
!     BASIS:       basis of the floating-point representation    *      
!                  used by the machine (in most cases 2 or 16)   *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     MAT:         the upper triangle of this (N,N) matrix       *      
!                  contains the eigenvectors of the quasi-trian- *      
!                  gular matrix, that is produced by the         *      
!                  QR-method.                                    *      
!     D:           N-vector with transform info from ORTHES      *      
!     VALR,VALI:   two N-vectors VALR(1:N), VALI(1:N) of type    *      
!                  DOUBLE PRECISION, that contain the real and   *      
!                  imaginary parts of the eigenvalues of MAT     *      
!     EIVEC:       2-dimensional array EIVEC(1:N,1:N) of type    *      
!                  DOUBLE PRECISION, that contains the normalized*      
!                  eigenvectors of the input matrix MAT as       *      
!                  columns in case VEC = .TRUE. .                *      
!                  If the I-th eigenvalue of MAT is real, the    *      
!                  I-th column of EIVEC contains the correspond- *      
!                  ing real eigenvector. If the I-th and (I+1)-st*      
!                  eigenvalues are a complex conjugate pair,     *      
!                  the I-th and (I+1)-st columns of EIVEC contain*      
!                  the real and imaginary parts of the associated*      
!                  eigenvector for the eigenvalue with positive  *      
!                  imaginary part.                               *      
!     CNT:         N-vector CNT(1:N) of type INTEGER containing  *      
!                  the number of iteration steps for each eigen- *      
!                  value. If two eigenvalues are found simul-    *      
!                  taneously as a complex conjugate pair, the    *      
!                  number of iteration steps used for both is    *      
!                  recorded with a positive sign for the first   *      
!                  and with a negative sign for the second       *      
!                  eigenvalue.                                   *      
!     SCAL:        N-vector SCAL(1:N) of type DOUBLE PRECISION   *      
!                  containing information about the permutations *      
!                  and the scaling factors used.                 *      
!     LOW,HIGH:    The rows numbered 1 to LOW-1 or HIGH+1 to N   *      
!                  contain isolated eigenvalues, i.e., eigen-    *      
!                  values for the eigenvectors e_i, the unit     *      
!                  vectors                                       *      
!                                                                *      
!     RETURN VALUES of FUNCTION EIGEN:                           *      
!     ================================                           *      
!     0:      no error, eigenvalue-eigenvector problem solved.   *      
!     401:    order N of the input matrix MAT is less than 1.    *      
!     402:    MAT is the zero matrix.                            *      
!     403:    the maximum number of steps for the QR-method      *      
!             has been reached. However, not all of the eigen-   *      
!             values have be determined.                         *      
!                                                                *      
!     LOCAL VARIABLES:                                           *      
!     ================                                           *      
!     ONE,TWO,HALF: floating-point constants                     *      
!     EPS:          machine constant                             *      
!     TEMP:         auxiliary variable                           *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: BALAN, ELMHES, ELMTRA, HQR2, BALBAK     *      
!                        NORMAL, SWAP, COMDIV, COMABS, ORTHES,   *      
!                        ORTHTRA                                 *      
!                                                                *      
!                                                                *      
!  sources : 1. Martin, R. S. and Wilkinson, J. H.,              *      
!               see [MART70].                                    *      
!            2. Parlett, B. N. and Reinsch, C., see [PARL69].    *      
!            3. Peters, G. and Wilkinson, J. H., see [PETE70].   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Juergen Dietel                                     *      
!  date     : 7.15.1993                                          *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER BASIS, LD, N, CNT (1:N), LOW, HIGH 
      LOGICAL VEC, ORTHO, EVNORM 
      DOUBLEPRECISION MAT (1:LD, 1:N), SCAL (1:N), EIVEC (1:LD, 1:N) 
      DOUBLEPRECISION VALR (1:N), VALI (1:N), D (1:N) 
      DOUBLEPRECISION ONE, TWO, HALF 
      PARAMETER (ONE = 1.0D0, TWO = 2.0D0, HALF = 0.5D0) 
      INTEGER RES, BALAN, ELMHES, ORTHES, ELMTRA, ORTTRA, HQR2 
      INTEGER BALBAK, NORMAL 
      DOUBLEPRECISION EPS, TEMP 
!                                                                       
!     determine the machine constant EPS, i.e., the smallest            
!     positive number for which the  1 + EPS > 1)  is true              
!                                                                       
      TEMP = TWO 
      EPS = ONE 
   10 IF (ONE.LT.TEMP) THEN 
         EPS = EPS * HALF 
         TEMP = ONE+EPS 
         GOTO 10 
      ENDIF 
      EPS = TWO * EPS 
      RES = BALAN (LD, N, MAT, SCAL, LOW, HIGH, BASIS) 
      IF (RES.NE.0) THEN 
         EIGEN = RES + 100 
         RETURN 
      ENDIF 
      IF (ORTHO) THEN 
         RES = ORTHES (LD, N, LOW, HIGH, MAT, D, EPS) 
      ELSE 
         RES = ELMHES (LD, N, LOW, HIGH, MAT, CNT) 
      ENDIF 
      IF (RES.NE.0) THEN 
         EIGEN = RES + 200 
         RETURN 
      ENDIF 
      IF (VEC) THEN 
         IF (ORTHO) THEN 
            RES = ORTTRA (LD, N, LOW, HIGH, MAT, D, EIVEC) 
         ELSE 
            RES = ELMTRA (LD, N, LOW, HIGH, MAT, CNT, EIVEC) 
         ENDIF 
         IF (RES.NE.0) THEN 
            EIGEN = RES + 300 
            RETURN 
         ENDIF 
      ENDIF 
      RES = HQR2 (VEC, LD, N, LOW, HIGH, MAT, VALR, VALI, EIVEC, CNT,   &
      EPS)                                                              
      IF (RES.NE.0) THEN 
         EIGEN = RES + 400 
         RETURN 
      ENDIF 
      IF (VEC) THEN 
         RES = BALBAK (LD, N, LOW, HIGH, SCAL, EIVEC) 
         IF (RES.NE.0) THEN 
            EIGEN = RES + 500 
            RETURN 
         ENDIF 
         IF (EVNORM) THEN 
            RES = NORMAL (LD, N, EIVEC, VALI) 
            IF (RES.NE.0) THEN 
               EIGEN = RES + 600 
               RETURN 
            ENDIF 
         ENDIF 
      ENDIF 
      EIGEN = 0 
      END FUNCTION EIGEN                            
!                                                                       
!                                                                       
      INTEGER FUNCTION BALAN (LD, N, MAT, SCAL, LOW, HIGH, BASIS) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     PROGRAM OBJECTIVE:                                         *      
!     ==================                                         *      
!     The procedure BALAN balances a given real matrix with      *      
!     respect to the column sum norm.                            *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     LD:       leading dimension of the matrix as defined in    *      
!               the calling program                              *      
!     N:        the order of the given square matrix             *      
!     MAT:      2-dimensional array MAT(1:N,1:N) containing the  *      
!               input matrix                                     *      
!     BASIS:    basis for the floating-point representation of   *      
!               the machine                                      *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     MAT:      the balanced matrix                              *      
!     LOW,HIGH: two INTEGER numbers, for which MAT(I,J) = 0      *      
!               if the following holds:                          *      
!               1. I > J, and                                    *      
!               2. J = 1, ..., LOW-1 or I = HIGH+1, ..., N       *      
!     SCAL:     N-vector SCAL(1:N) containing information about  *      
!               permutations and scaling factors used.           *      
!                                                                *      
!     RETURN VALUE of SUBROUTINE BALAN:                          *      
!     =================================                          *      
!     0:             no error                                    *      
!                                                                *      
!     LOCAL VARIABLES:                                           *      
!     ============== =                                           *      
!     ZERO,ONE,PT95: floating-point constants                    *      
!     I,J,K,L:       counting variables                          *      
!     B2:            square of the machine floating point basis  *      
!     R,C,F,G,S:     auxiliary variables for determining row     *      
!                    norms, reciprocals etc.                     *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: SWAP                                    *      
!                                                                *      
!                                                                *      
!  sources : Parlett, B. N. and Reinsch, C., see [PARL69].       *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Juergen Dietel                                     *      
!  date     : 04.10.1987                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER LD, N, LOW, HIGH, BASIS 
      DOUBLEPRECISION SCAL (1:N), MAT (1:LD, 1:N) 
      DOUBLEPRECISION ZERO, ONE, PT95 
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, PT95 = 0.95D0) 
      INTEGER I, J, K, L, B2 
      DOUBLEPRECISION R, C, F, G, S 
!                                                                       
!     reduce the norm of MAT(1:N,1:N) by exact diagonal similarity      
!     transformations and store in SCAL(1:N)                            
!                                                                       
      B2 = BASIS * BASIS 
      L = 1 
      K = N 
!                                                                       
!     search for rows isolating an eigenvalue and                       
!     move them to the bottom                                           
!                                                                       
   10 DO 50 J = K, 1, - 1 
         R = ZERO 
         DO 20 I = 1, K 
   20    IF (I.NE.J) R = R + ABS (MAT (J, I) ) 
         IF (R.EQ.ZERO) THEN 
            SCAL (K) = J 
            IF (J.NE.K) THEN 
               DO 30 I = 1, K 
   30          CALL SWAP (MAT (I, J), MAT (I, K) ) 
               DO 40 I = L, N 
   40          CALL SWAP (MAT (J, I), MAT (K, I) ) 
            ENDIF 
            K = K - 1 
            GOTO 10 
         ENDIF 
   50 END DO 
!                                                                       
!     search for columns isolating an eigenvalue and                    
!     move them to the left                                             
!                                                                       
   60 DO 100 J = L, K 
         C = ZERO 
         DO 70 I = L, K 
   70    IF (I.NE.J) C = C + ABS (MAT (I, J) ) 
         IF (C.EQ.ZERO) THEN 
            SCAL (L) = J 
            IF (J.NE.L) THEN 
               DO 80 I = 1, K 
   80          CALL SWAP (MAT (I, J), MAT (I, L) ) 
               DO 90 I = L, N 
   90          CALL SWAP (MAT (J, I), MAT (L, I) ) 
            ENDIF 
            L = L + 1 
            GOTO 60 
         ENDIF 
  100 END DO 
!                                                                       
!     balance the matrix in rows L to K                                 
!                                                                       
      LOW = L 
      HIGH = K 
      DO 110 I = L, K 
  110 SCAL (I) = ONE 
  120 DO 180 I = L, K 
         C = ZERO 
         R = ZERO 
         DO 130 J = L, K 
            IF (J.NE.I) THEN 
               C = C + ABS (MAT (J, I) ) 
               R = R + ABS (MAT (I, J) ) 
            ENDIF 
  130    END DO 
         G = R / BASIS 
         F = ONE 
         S = C + R 
  140    IF (C.LT.G) THEN 
            F = F * BASIS 
            C = C * B2 
            GOTO 140 
         ENDIF 
         G = R * BASIS 
  150    IF (C.GE.G) THEN 
            F = F / BASIS 
            C = C / B2 
            GOTO 150 
         ENDIF 
         IF ( (C + R) / F.LT.PT95 * S) THEN 
            G = ONE / F 
            SCAL (I) = SCAL (I) * F 
            DO 160 J = L, N 
  160       MAT (I, J) = MAT (I, J) * G 
            DO 170 J = 1, K 
  170       MAT (J, I) = MAT (J, I) * F 
            GOTO 120 
         ENDIF 
  180 END DO 
      BALAN = 0 
      END FUNCTION BALAN                            
!                                                                       
!                                                                       
      INTEGER FUNCTION BALBAK (LD, N, LOW, HIGH, SCAL, EIVEC) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     PROGRAM OBJECTIVE:                                         *      
!     ==================                                         *      
!     BALBAK transforms all right eigenvectors of a balanced     *      
!     matrix back into the eigenvectors of the original matrix.  *      
!     The balanced matrix was generated by calling the procedure *      
!     BALAN.                                                     *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     LD:       leading dimension of the array EIVEC as defined  *      
!               in the calling program                           *      
!     N:        order of the eigenvectors (number of components) *      
!     LOW,HIGH: two INTEGER numbers that were created in the     *      
!               procedure BALAN                                  *      
!     SCAL:     output vector of the procedure BALAN             *      
!     EIVEC:    2-dimensional array EIVEC(1:N,1:N), each column  *      
!               of EIVEC contains an eigenvector (or its real    *      
!               or imaginary parts) for the balanced matrix      *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     EIVEC:    the eigenvectors (or their real or imaginary     *      
!               parts) of the original matrix                    *      
!                                                                *      
!     RETURN VALUE of SUBROUTINE BALBAK:                         *      
!     ==================================                         *      
!     0:        no error                                         *      
!                                                                *      
!     LOCAL VARIABLES:                                           *      
!     ================                                           *      
!     I,J,K:    auxiliary variables used for indexing            *      
!     S:        scaling constant                                 *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: SWAP                                    *      
!                                                                *      
!                                                                *      
!  sources : Parlett, B. N. and Reinsch, C., see [PARL69].       *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Juergen Dietel                                     *      
!  date     : 04.10.1987                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER LD, N, LOW, HIGH 
      DOUBLEPRECISION SCAL (1:N), EIVEC (1:LD, 1:N) 
      INTEGER I, J, K 
      DOUBLEPRECISION S 
!                                                                       
      DO 20 I = LOW, HIGH 
         S = SCAL (I) 
!                                                                       
!        left eigenvectors are back transformed, if the last            
!        statement is replaced by:  S = 1.0D0/SCAL(I)                   
!                                                                       
         DO 10 J = 1, N 
   10    EIVEC (I, J) = EIVEC (I, J) * S 
   20 END DO 
      DO 40 I = LOW - 1, 1, - 1 
         K = SCAL (I) 
         DO 30 J = 1, N 
   30    CALL SWAP (EIVEC (I, J), EIVEC (K, J) ) 
   40 END DO 
      DO 60 I = HIGH + 1, N 
         K = SCAL (I) 
         DO 50 J = 1, N 
   50    CALL SWAP (EIVEC (I, J), EIVEC (K, J) ) 
   60 END DO 
      BALBAK = 0 
      END FUNCTION BALBAK                           
!                                                                       
!                                                                       
      INTEGER FUNCTION ELMHES (LD, N, LOW, HIGH, MAT, PERM) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     PROGRAM OBJECTIVE:                                         *      
!     ==================                                         *      
!     For a general square matrix A(1:N,1:N) this program reduces*      
!     the principal submatrix of A of order HIGH-LOW+1, defined  *      
!     to lie between the elements A(LOW,LOW) and A(HIGH,HIGH),   *      
!     to Hessenberg form H by similarity with non-orthogonal     *      
!     Gaussian elimination matrices. The principal submatrix is  *      
!     overwriten with H. The transforming matrices are stored in *      
!     the triangle below H and in the vector PERM.               *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     LD:       leading dimension of the matrix as defined in    *      
!               the calling program                              *      
!     N:        order of the square matrix A                     *      
!     LOW,HIGH: output parameter of the balancing procedure      *      
!               BALAN. If A has not been balanced, set LOW:=1,   *      
!               and HIGH:=N.                                     *      
!     MAT:      2-dimensional array MAT(1:N,1:N), the matrix in  *      
!               balanced form                                    *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     MAT:      2-dimensional array MAT(1:N,1:n), which consists *      
!               in part of the upper Hessenberg matrix and the   *      
!               transforming matrices.                           *      
!               The number N(I,R+1), that is needed for the      *      
!               reduction, is stored in the (I,R) position.      *      
!     PERM:     INTEGER vector that stores the row- and column   *      
!               permutations performed during the reduction      *      
!                                                                *      
!     RETURN VALUE:                                              *      
!     =============                                              *      
!     0:        no error                                         *      
!                                                                *      
!     LOCAL VARIABLES:                                           *      
!     ================                                           *      
!     ZERO,ONE: floating-point constants                         *      
!     I,J,M:    counters                                         *      
!     X,Y:      auxiliary variables used for storing matrix      *      
!               elements and intermediate results                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: SWAP                                    *      
!                                                                *      
!                                                                *      
!  sources : Martin, R. S. and Wilkinson, J. H., see [MART70].   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Juergen Dietel                                     *      
!  date     : 04.10.1987                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER N, LOW, HIGH, PERM (1:N) 
      DOUBLEPRECISION MAT (1:LD, 1:N) 
      DOUBLEPRECISION ZERO, ONE 
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0) 
      INTEGER I, J, M 
      DOUBLEPRECISION X, Y 
!                                                                       
      DO 70 M = LOW + 1, HIGH - 1 
         I = M 
         X = ZERO 
         DO 10 J = M, HIGH 
            IF (ABS (MAT (J, M - 1) ) .GT.ABS (X) ) THEN 
               X = MAT (J, M - 1) 
               I = J 
            ENDIF 
   10    END DO 
         PERM (M) = I 
         IF (I.NE.M) THEN 
!                                                                       
!           Swap rows and columns of MAT                                
!                                                                       
            DO 20 J = M - 1, N 
   20       CALL SWAP (MAT (I, J), MAT (M, J) ) 
            DO 30 J = 1, HIGH 
   30       CALL SWAP (MAT (J, I), MAT (J, M) ) 
         ENDIF 
         IF (X.NE.ZERO) THEN 
            DO 60 I = M + 1, HIGH 
               Y = MAT (I, M - 1) 
               IF (Y.NE.ZERO) THEN 
                  Y = Y / X 
                  MAT (I, M - 1) = Y 
                  DO 40 J = M, N 
   40             MAT (I, J) = MAT (I, J) - Y * MAT (M, J) 
                  DO 50 J = 1, HIGH 
   50             MAT (J, M) = MAT (J, M) + Y * MAT (J, I) 
               ENDIF 
   60       END DO 
         ENDIF 
   70 END DO 
      ELMHES = 0 
      END FUNCTION ELMHES                           
!                                                                       
!                                                                       
      INTEGER FUNCTION ELMTRA (LD, N, LOW, HIGH, MAT, PERM, H) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     PROGRAM OBJECTIVE:                                         *      
!     ==================                                         *      
!     Form the matrix of accumulated transformations from the    *      
!     information left by procedure ELMHES in the lower triangle *      
!     of the Hessenberg matrix H - in MAT(1:N,1:N) and in the    *      
!     INTEGER vector PERM(1:N). Store in the array H(1:N,1:N).   *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     LD:       leading dimension of the matrix as defined in    *      
!               the calling program                              *      
!     N:        order of the Hessenberg matrix H                 *      
!     LOW,HIGH: INTEGER numbers, that were produced by procedure *      
!               BALAN (if it was used; otherwise set LOW:=1,     *      
!               HIGH:=N.)                                        *      
!     PERM:     INTEGER N-vector produced by ELMHES              *      
!     MAT:      (N,N)-matrix, that was produced by ELMHES and    *      
!               contains the Hessenberg matrix H and the LR      *      
!               multipliers used in order to create H from the   *      
!               given matrix A                                   *      
!                                                                *      
!     OUTPUT PARAMETER:                                          *      
!     =================                                          *      
!     H:        (N,N)-matrix that describes the similarity       *      
!               transformation of A to Hessenberg form H         *      
!                                                                *      
!     RETURN VALUE:                                              *      
!     =============                                              *      
!     0:             no error                                    *      
!                                                                *      
!     LOCAL FACTORS:                                             *      
!     ==============                                             *      
!     ZERO,ONE: floating-point constants                         *      
!     I,J,K:    index variables                                  *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!                                                                *      
!  sources : Peters, G. and Wilkinson, J. H., see [PETE70].      *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Juergen Dietel                                     *      
!  date     : 04.10.1987                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER LD, N, LOW, HIGH, PERM (1:N) 
      DOUBLEPRECISION MAT (1:LD, 1:N), H (1:LD, 1:N) 
      DOUBLEPRECISION ZERO, ONE 
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0) 
      INTEGER I, J, K 
!                                                                       
      DO 20 I = 1, N 
         DO 10 J = 1, N 
   10    H (I, J) = ZERO 
         H (I, I) = ONE 
   20 END DO 
!                                                                       
      DO 50 I = HIGH - 1, LOW + 1, - 1 
         J = PERM (I) 
         DO 30 K = I + 1, HIGH 
   30    H (K, I) = MAT (K, I - 1) 
         IF (I.NE.J) THEN 
            DO 40 K = I, HIGH 
               H (I, K) = H (J, K) 
               H (J, K) = ZERO 
   40       END DO 
            H (J, I) = ONE 
         ENDIF 
   50 END DO 
      ELMTRA = 0 
      END FUNCTION ELMTRA                           
!                                                                       
!                                                                       
      INTEGER FUNCTION ORTHES (LD, N, LOW, HIGH, MAT, D, EPSM) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This function transform the matrix MAT to upper Hessenberg    *      
!  using Householder transforms.                                 *      
!  The essential transform information is stored in the otherwise*      
!  unused lower triangle of MAT and in the vector D.             *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! LD       leading dimension of the matrix MAT, as defined in the*      
!          calling program                                       *      
! N        order of the matrix MAT                               *      
! LOW  \   the rows 1 to LOW-1 and HIGH+1 to N contain the       *      
! HIGH  >  isolated eigenvalues, i. e., those eigenvalues that   *      
!      /   have unit vectors e_i as eigenvectors.                *      
! MAT      the original matrix MAT(1:N,1:N)                      *      
! EPSM     machine constant                                      *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! MAT      the desired Hessenberg matrix with part of the        *      
!          transform information below the subdiagonal           *      
! D        N-vector with the remainder of the transform info     *      
!                                                                *      
! RETURN VALUE :                                                 *      
! ==============                                                 *      
! Error code. Only  0 (no error) is possible here.               *      
!                                                                *      
! LOCAL VARIABLES:                                               *      
! ================                                               *      
! I,J,M  Loop variables                                          *      
! S      Euclidean norm sigma of a column below the subdiagonal  *      
!        v of MAT, which must be transformed into a multiple of  *      
!        e1 = (1,0,...,0); (v = (v1,...,v(HIGH-M+1))             *      
! X      initially leading element of v, then summation value    *      
!        inside the Householder transformation                   *      
! Y      initially  sigma^2, then ||u||^2, where                 *      
!        u := v +- sigma * e1                                    *      
! EPS    accuracy bound to check transformation                  *      
! ZERO   double precision zero 0.0D0                             *      
!                                                                *      
!*****************************************************************      
! Reference: R.S. Martin and J.H. Wilkinson, Num. Math. 12 (1968)*      
!            pp. 359, 360, see [WILK71], contrib. II/13          *      
! Author:    Juergen Dietel, Computer Center, RWTH Aachen        *      
! Date:      7. 15. 1993                                         *      
!*****************************************************************      
!                                                                       
      INTEGER LD, N, LOW, HIGH 
      DOUBLEPRECISION MAT (1:LD, 1:N), D (1:N), EPSM 
!                                                                       
      INTEGER I, J, M 
      DOUBLEPRECISION S 
      DOUBLEPRECISION X 
      DOUBLEPRECISION Y 
      DOUBLEPRECISION EPS 
!                                                                       
      DOUBLEPRECISION ZERO 
      PARAMETER (ZERO = 0.0D0) 
!                                                                       
      EPS = 128 * EPSM 
!                                                                       
      DO 80 M = LOW + 1, HIGH - 1 
!                                                                       
         Y = ZERO 
         DO 10 I = HIGH, M, - 1 
            X = MAT (I, M - 1) 
            D (I) = X 
            Y = Y + X * X 
   10    END DO 
         IF (Y.LE.EPS) THEN 
            S = ZERO 
         ELSE 
!                                                                       
            IF (X.GE.ZERO) THEN 
               S = - SQRT (Y) 
            ELSE 
               S = SQRT (Y) 
            ENDIF 
            Y = Y - X * S 
            D (M) = X - S 
!                                                                       
!           Multiply MAT on the left by  (E-(u * uT)/y)                 
!                                                                       
            DO 40 J = M, N 
               X = ZERO 
               DO 20 I = HIGH, M, - 1 
                  X = X + D (I) * MAT (I, J) 
   20          END DO 
               X = X / Y 
               DO 30 I = M, HIGH 
                  MAT (I, J) = MAT (I, J) - X * D (I) 
   30          END DO 
   40       END DO 
!                                                                       
!           Multiply MAT on the right by (E-(u * uT)/y)                 
!                                                                       
            DO 70 I = 1, HIGH 
               X = ZERO 
               DO 50 J = HIGH, M, - 1 
                  X = X + D (J) * MAT (I, J) 
   50          END DO 
               X = X / Y 
               DO 60 J = M, HIGH 
                  MAT (I, J) = MAT (I, J) - X * D (J) 
   60          END DO 
   70       END DO 
!                                                                       
         ENDIF 
!                                                                       
         MAT (M, M - 1) = S 
   80 END DO 
!                                                                       
      ORTHES = 0 
!                                                                       
      END FUNCTION ORTHES                           
!                                                                       
!                                                                       
      INTEGER FUNCTION ORTTRA (LD, N, LOW, HIGH, MAT, D, V) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  Reconstruct the transformation matrix V of the Householder    *      
!  reductions of MAT to Hessenberg form from the information     *      
!  stored in the lower triangle of MAT and in D.                 *      
!  The contents of D is lost.                                    *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! LD       leading dimension of the matrices  MAT and V, as      *      
!          defined in the calling routine                        *      
! N        size of the matrix MAT                                *      
! LOW  \   the rows  1 to LOW-1 and the rows  HIGH+1  to N       *      
! HIGH  >  contain the isolated eigenvalues, i.e. those with     *      
!      /   unit vectors e_i as eigenvectors.                     *      
! MAT      (1:N,1:N) matrix, that ORTHES has reduced to upper    *      
!          Hessenberg form; partially filled with transformation *      
!          information                                           *      
! D        N-vector with the remaining transform information     *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! D        modified input vector                                 *      
! V        (1:N,1:N)-matrix, giving the similarity transform     *      
!          from A to the upper Hessenberg matrix in MAT          *      
!                                                                *      
! RETURN VALUE :                                                 *      
! ==============                                                 *      
! Error code. Only 0 possible here (no error).                   *      
!                                                                *      
! LOCAL VARIABLES:                                               *      
! ================                                               *      
! I,J,M  loop variables                                          *      
! X      summation variable during Householder transformation    *      
! Y      sigma  or  sigma * (v1 +- sigma)                        *      
! ZERO   double precision  0.0D0                                 *      
! ONE    double precision  1.0D0                                 *      
!                                                                *      
!*****************************************************************      
! Reference: G. Peters and J.H. Wilkinson, Num. Math. 16 (1970)  *      
!            p. 191, see [WILK71], contrib. II/15                *      
! Author:    Juergen Dietel, Computer Center, RWTH Aachen        *      
! Date:      7. 15. 1993                                         *      
!*****************************************************************      
!                                                                       
      INTEGER LD, N, LOW, HIGH 
      DOUBLEPRECISION MAT (1:LD, 1:N), D (1:N), V (1:LD, 1:N) 
!                                                                       
      INTEGER I, J, M 
      DOUBLEPRECISION X 
      DOUBLEPRECISION Y 
!                                                                       
      DOUBLEPRECISION ZERO, ONE 
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0) 
!                                                                       
!     start with identity matrix in V                                   
!                                                                       
      DO 20 I = 1, N 
         DO 10 J = 1, N 
            V (I, J) = ZERO 
   10    END DO 
         V (I, I) = ONE 
   20 END DO 
!                                                                       
!     Perform the transformations that have reduced MAT to              
!     Hessenberg form on the unit matrix stored in V in order           
!     to generate the transformation matrix.                            
!                                                                       
      DO 70 M = HIGH - 1, LOW + 1, - 1 
         Y = MAT (M, M - 1) 
!                                                                       
         IF (Y.NE.ZERO) THEN 
            Y = Y * D (M) 
            DO 30 I = M + 1, HIGH 
               D (I) = MAT (I, M - 1) 
   30       END DO 
            DO 60 J = M, HIGH 
               X = ZERO 
               DO 40 I = M, HIGH 
                  X = X + D (I) * V (I, J) 
   40          END DO 
               X = X / Y 
               DO 50 I = M, HIGH 
                  V (I, J) = V (I, J) + X * D (I) 
   50          END DO 
   60       END DO 
         ENDIF 
   70 END DO 
!                                                                       
      ORTTRA = 0 
!                                                                       
      END FUNCTION ORTTRA                           
!                                                                       
!                                                                       
      INTEGER FUNCTION HQR2 (VEC, LD, N, LOW, HIGH, H, VALR, VALI,      &
      EIVEC, CNT, EPS)                                                  
!                                                                       
!*****************************************************************      
!                                                                *      
!     PROGRAM OBJECTIVE:                                         *      
!     ==================                                         *      
!     Finds the eigenvalues and eigenvectors (if VEC = .TRUE.) of*      
!     a real matrix, which has been reduced to upper Hessenberg  *      
!     form and is stored in the array H(1:N,1:N) with the accu-  *      
!     mulated transformations stored in the array EIVEC(1:N,1:N).*      
!     The real and imaginary parts of the eigenvalues are placed *      
!     in the two vectors VALR(1:N), VALI(1:N) while the eigen-   *      
!     vectors are stored in the array EIVEC(1:N,1:N), where only *      
!     one complex eigenvector corresponding to the eigenvalue    *      
!     with positive imaginary part is stored for a complex con-  *      
!     jugate eigenvalue pair. LOW and HIGH are two INTEGER       *      
!     numbers produced during balancing MAT, so that eigenvalues *      
!     are isolated in positions 1 to LOW-1 and HIGH+1 to N. If   *      
!     no initial balancing was performed, set LOW:=1, HIGH:=N.   *      
!     The subroutine is aborted with an error message if any one *      
!     of the eigenvalues requires more than MAXSTP iterations.   *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     VEC:      parameter indicating eigenvector computation:    *      
!                 =  .TRUE.  : compute eigenvectors              *      
!                 =  .FALSE. : compute eigenvalues only          *      
!     LD:       leading dimension of H and EIVEC, as defined in  *      
!               the calling program                              *      
!     N:        order of the Hessenberg matrix H                 *      
!     LOW,HIGH: INTEGER numbers produced by BALAN, if it was     *      
!               used. Otherwise set LOW:=1, HIGH:=N.             *      
!     EPS:      machine constant                                 *      
!     H:        (N,N) matrix containing H with its relevant      *      
!               components                                       *      
!     EIVEC:    (N,N) matrix containing the array that describes *      
!               the similarity transformation of A to H..        *      
!               (This was procuded by ELMTRA or ORTTRA.)         *      
!               If H is the original matrix, define EIVEC := I,  *      
!               the identity matrix.                             *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     H:           the upper triangle of this (N,N) matrix       *      
!                  contains the eigenvectors of the quasi-trian- *      
!                  gular matrix, that is produced by the QR steps*      
!     VALR,VALI:   two N-vectors, with the real and imaginary    *      
!                  parts of the eigenvalues                      *      
!     CNT:         INTEGER N-vector, with the count of iterations*      
!                  for each eigenvalue. If two eigenvalues are   *      
!                  found as a complex conjugate pair, the number *      
!                  of iterations is entered with a positive sign *      
!                  for the first and with a negative sign for the*      
!                  second eigenvalue.                            *      
!     EIVEC:       (N,N) matrix, where (for VEC = .TRUE.) the    *      
!                  non-normalized eigenvectors of the original   *      
!                  matrix are stored in case H is not the        *      
!                  original matrix.                              *      
!                  If the I-th eigenvalue is real, the I-th      *      
!                  column of EIVEC contains the corresponding    *      
!                  real eigenvector. If the eigenvalues numbered *      
!                  I and I+1 form a complex conjugate eigenvalue *      
!                  pair for MAT, then the I-th and (I+1)-st      *      
!                  columns of EIVEC contain the real and         *      
!                  imaginary parts of the eigenvector for the    *      
!                  corresponding eigenvalue with positive real   *      
!                  part.                                         *      
!                                                                *      
!                                                                *      
!     RETURN VALUES:                                             *      
!     ==============                                             *      
!     0:           no error                                      *      
!     1:           parameters N, LOW or HIGH are unacceptable    *      
!     2:           all eigenvectors are the zero vector.         *      
!     3:           the maximum number of QR-steps has been       *      
!                  exceeded.                                     *      
!                                                                *      
!     LOCAL VARIABLES:                                           *      
!     ================                                           *      
!     ZERO,ONE,TWO,PT75,PT4375: floating-point constants         *      
!     MAXSTP:                   maximum number of allowed steps  *      
!     I,J,K,L,M,N,NA,EN:        index variables                  *      
!     ITER:                     step counter                     *      
!     P,Q,R,S,T,W,X,Y,Z,NORM,                                    *      
!       RA,SA,VR,VI:            auxiliary variables              *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: COMDIV                                  *      
!                                                                *      
!                                                                *      
!  sources : Peters, G. and Wilkinson, J. H., see [PETE70].      *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Juergen Dietel                                     *      
!  date     : 07.14.1993                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      LOGICAL VEC 
      INTEGER LD, N, LOW, HIGH, CNT (1:N) 
      DOUBLEPRECISION H (1:LD, 1:N), EIVEC (1:LD, 1:N), VALR (1:N) 
      DOUBLEPRECISION VALI (1:N), EPS 
      DOUBLEPRECISION ZERO, ONE, TWO, PT75, PT4375 
      INTEGER MAXSTP 
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, PT75 = 0.75D0, &
      PT4375 = 0.4375D0, MAXSTP = 100)                                  
      INTEGER I, J, K, L, M, NA, EN, ITER 
      DOUBLEPRECISION P, Q, R, S, T, W, X, Y, Z, NORM, RA, SA, VR, VI 
!                                                                       
!     error 1: one of the parameters N, LOW or HIGH has an              
!              unacceptable value                                       
!                                                                       
      IF (N.LT.1.OR.LOW.LT.1.OR.HIGH.GT.N) THEN 
         HQR2 = 1 
         RETURN 
      ENDIF 
!                                                                       
!     initialize the isolated eigenvalues found during balancing:       
!                                                                       
      DO 10 I = 1, N 
         IF (I.LT.LOW.OR.I.GT.HIGH) THEN 
            VALR (I) = H (I, I) 
            VALI (I) = ZERO 
            CNT (I) = 0 
         ELSE 
            CNT (I) = - 1 
         ENDIF 
   10 END DO 
!                                                                       
      EN = HIGH 
      T = ZERO 
   15 IF (EN.LT.LOW) GOTO 333 
      ITER = 0 
      NA = EN - 1 
!                                                                       
!           search for a single small subdiagonal element:              
!                                                                       
   20 DO 30 L = EN, LOW + 1, - 1 
   30 IF (ABS (H (L, L - 1) ) .LE.EPS * (ABS (H (L - 1, L - 1) )        &
      + ABS (H (L, L) ) ) ) GOTO 40                                     
   40 X = H (EN, EN) 
      IF (L.EQ.EN) THEN 
!                                                                       
!              found one root:                                          
!                                                                       
         VALR (EN) = X + T 
         H (EN, EN) = VALR (EN) 
         VALI (EN) = ZERO 
         CNT (EN) = ITER 
         EN = NA 
         GOTO 15 
      ENDIF 
!                                                                       
      Y = H (NA, NA) 
      W = H (EN, NA) * H (NA, EN) 
      IF (L.EQ.NA) THEN 
!                                                                       
!              found two roots:                                         
!                                                                       
         P = (Y - X) / TWO 
         Q = P * P + W 
         Z = SQRT (ABS (Q) ) 
         H (EN, EN) = X + T 
         X = H (EN, EN) 
         H (NA, NA) = Y + T 
         CNT (EN) = - ITER 
         CNT (NA) = ITER 
         IF (Q.GE.ZERO) THEN 
!                                                                       
!                 found a real pair:                                    
!                                                                       
            IF (P.LT.ZERO) Z = - Z 
            Z = P + Z 
            VALR (NA) = X + Z 
            VALR (EN) = X - W / Z 
            VALI (NA) = ZERO 
            VALI (EN) = ZERO 
            X = H (EN, NA) 
            R = SQRT (X * X + Z * Z) 
                                                                        
            IF (VEC) THEN 
               P = X / R 
               Q = Z / R 
!                                                                       
!                    row modification:                                  
!                                                                       
               DO 50 J = NA, N 
                  Z = H (NA, J) 
                  H (NA, J) = Q * Z + P * H (EN, J) 
                  H (EN, J) = Q * H (EN, J) - P * Z 
   50          END DO 
!                                                                       
!                    column modification:                               
!                                                                       
               DO 60 I = 1, EN 
                  Z = H (I, NA) 
                  H (I, NA) = Q * Z + P * H (I, EN) 
                  H (I, EN) = Q * H (I, EN) - P * Z 
   60          END DO 
!                                                                       
!                    accumulate transformations:                        
!                                                                       
               DO 70 I = LOW, HIGH 
                  Z = EIVEC (I, NA) 
                  EIVEC (I, NA) = Q * Z + P * EIVEC (I, EN) 
                  EIVEC (I, EN) = Q * EIVEC (I, EN) - P * Z 
   70          END DO 
            ENDIF 
         ELSE 
!                                                                       
!                 complex conjugate pair:                               
!                                                                       
            VALR (NA) = X + P 
            VALR (EN) = VALR (NA) 
            VALI (NA) = Z 
            VALI (EN) = - Z 
         ENDIF 
         EN = EN - 2 
         GOTO 15 
      ENDIF 
!                                                                       
      IF (ITER.EQ.MAXSTP) THEN 
!                                                                       
!           error 3: maximal number of iterations exceeded:             
!                                                                       
         CNT (EN) = MAXSTP + 1 
         HQR2 = 3 
         RETURN 
      ENDIF 
      IF (MOD (ITER, 10) .EQ.0.AND.ITER.NE.0) THEN 
!                                                                       
!              use exceptional shift:                                   
!                                                                       
         T = T + X 
         DO 80 I = LOW, EN 
   80    H (I, I) = H (I, I) - X 
         S = ABS (H (EN, NA) ) + ABS (H (NA, EN - 2) ) 
         X = PT75 * S 
         Y = X 
         W = - PT4375 * S * S 
      ENDIF 
      ITER = ITER + 1 
!                                                                       
!        search for two consecutive small subdiagonal elements:         
!                                                                       
      DO 90 M = EN - 2, L, - 1 
         Z = H (M, M) 
         R = X - Z 
         S = Y - Z 
         P = (R * S - W) / H (M + 1, M) + H (M, M + 1) 
         Q = H (M + 1, M + 1) - Z - R - S 
         R = H (M + 2, M + 1) 
         S = ABS (P) + ABS (Q) + ABS (R) 
         P = P / S 
         Q = Q / S 
         R = R / S 
         IF (M.EQ.L) GOTO 100 
         IF (ABS (H (M, M - 1) ) * (ABS (Q) + ABS (R) ) .LE.EPS * ABS ( &
         P) * (ABS (H (M - 1, M - 1) ) + ABS (Z) + ABS (H (M + 1, M + 1)&
         ) ) ) GOTO 100                                                 
   90 END DO 
  100 DO 110 I = M + 2, EN 
  110 H (I, I - 2) = ZERO 
      DO 120 I = M + 3, EN 
  120 H (I, I - 3) = ZERO 
!                                                                       
!           double QR-step involving rows L to EN and                   
!           columns M to EN of the complete array:                      
!                                                                       
      DO 200 K = M, NA 
         IF (K.NE.M) THEN 
            P = H (K, K - 1) 
            Q = H (K + 1, K - 1) 
            IF (K.NE.NA) THEN 
               R = H (K + 2, K - 1) 
            ELSE 
               R = ZERO 
            ENDIF 
            X = ABS (P) + ABS (Q) + ABS (R) 
            IF (X.EQ.ZERO) GOTO 200 
            P = P / X 
            Q = Q / X 
            R = R / X 
         ENDIF 
         S = SQRT (P * P + Q * Q + R * R) 
         IF (P.LT.ZERO) S = - S 
         IF (K.NE.M) THEN 
            H (K, K - 1) = - S * X 
         ELSEIF (L.NE.M) THEN 
            H (K, K - 1) = - H (K, K - 1) 
         ENDIF 
         P = P + S 
         X = P / S 
         Y = Q / S 
         Z = R / S 
         Q = Q / P 
         R = R / P 
!                                                                       
!              row modification:                                        
!                                                                       
         DO 130 J = K, N 
            P = H (K, J) + Q * H (K + 1, J) 
            IF (K.NE.NA) THEN 
               P = P + R * H (K + 2, J) 
               H (K + 2, J) = H (K + 2, J) - P * Z 
            ENDIF 
            H (K + 1, J) = H (K + 1, J) - P * Y 
            H (K, J) = H (K, J) - P * X 
  130    END DO 
         J = MIN (K + 3, EN) 
!                                                                       
!              column modification:                                     
!                                                                       
         DO 140 I = 1, J 
            P = X * H (I, K) + Y * H (I, K + 1) 
            IF (K.NE.NA) THEN 
               P = P + Z * H (I, K + 2) 
               H (I, K + 2) = H (I, K + 2) - P * R 
            ENDIF 
            H (I, K + 1) = H (I, K + 1) - P * Q 
            H (I, K) = H (I, K) - P 
  140    END DO 
!                                                                       
         IF (VEC) THEN 
!                                                                       
!                 accumulate transformations:                           
!                                                                       
            DO 150 I = LOW, HIGH 
               P = X * EIVEC (I, K) + Y * EIVEC (I, K + 1) 
               IF (K.NE.NA) THEN 
                  P = P + Z * EIVEC (I, K + 2) 
                  EIVEC (I, K + 2) = EIVEC (I, K + 2) - P * R 
               ENDIF 
               EIVEC (I, K + 1) = EIVEC (I, K + 1) - P * Q 
               EIVEC (I, K) = EIVEC (I, K) - P 
  150       END DO 
         ENDIF 
  200 END DO 
      GOTO 20 
!                                                                       
!                                                                       
  333 IF (.NOT.VEC) THEN 
         HQR2 = 0 
         RETURN 
      ENDIF 
!                                                                       
!                                                                       
!     all eigenvalues have been found; now transform back:              
!                                                                       
!     find the 1-norm of H :                                            
!                                                                       
      NORM = ZERO 
      K = 1 
      DO 201 I = 1, N 
         DO 101 J = K, N 
  101    NORM = NORM + ABS (H (I, J) ) 
  201 K = I 
      IF (NORM.EQ.ZERO) THEN 
!        Fehler 2: 1-Norm von H ist gleich 0:                           
         HQR2 = 2 
         RETURN 
      ENDIF 
!                                                                       
!     back transformation:                                              
!                                                                       
      DO 207 EN = N, 1, - 1 
         P = VALR (EN) 
         Q = VALI (EN) 
         NA = EN - 1 
         IF (Q.EQ.ZERO) THEN 
!                                                                       
!           real vector:                                                
!                                                                       
            M = EN 
            H (EN, EN) = ONE 
            DO 63 I = NA, 1, - 1 
               W = H (I, I) - P 
               R = H (I, EN) 
               DO 38 J = M, NA 
   38          R = R + H (I, J) * H (J, EN) 
               IF (VALI (I) .LT.ZERO) THEN 
                  Z = W 
                  S = R 
               ELSE 
                  M = I 
                  IF (VALI (I) .EQ.ZERO) THEN 
                     IF (W.NE.ZERO) THEN 
                        H (I, EN) = - R / W 
                     ELSE 
                        H (I, EN) = - R / (EPS * NORM) 
                     ENDIF 
                  ELSE 
!                                                                       
!                    solve the linear system:                           
!                    [ W   X ] [ H(I,EN)   ]   [ -R ]                   
!                    [       ] [           ] = [    ]                   
!                    [ Y   Z ] [ H(I+1,EN) ]   [ -S ]                   
!                                                                       
                     X = H (I, I + 1) 
                     Y = H (I + 1, I) 
                     Q = (VALR (I) - P) * (VALR (I) - P) + VALI (I)     &
                     * VALI (I)                                         
                     T = (X * S - Z * R) / Q 
                     H (I, EN) = T 
                     IF (ABS (X) .GT.ABS (Z) ) THEN 
                        H (I + 1, EN) = ( - R - W * T) / X 
                     ELSE 
                        H (I + 1, EN) = ( - S - Y * T) / Z 
                     ENDIF 
                  ENDIF 
               ENDIF 
   63       END DO 
         ELSEIF (Q.LT.ZERO) THEN 
!                                                                       
!           complexer eigenvector for LAMBDA = P - I * Q :              
!                                                                       
            M = NA 
            IF (ABS (H (EN, NA) ) .GT.ABS (H (NA, EN) ) ) THEN 
               H (NA, NA) = - (H (EN, EN) - P) / H (EN, NA) 
               H (NA, EN) = - Q / H (EN, NA) 
            ELSE 
               CALL COMDIV ( - H (NA, EN), ZERO, H (NA, NA) - P, Q, H ( &
               NA, NA), H (NA, EN) )                                    
            ENDIF 
            H (EN, NA) = ONE 
            H (EN, EN) = ZERO 
            DO 190 I = NA - 1, 1, - 1 
               W = H (I, I) - P 
               RA = H (I, EN) 
               SA = ZERO 
               DO 75 J = M, NA 
                  RA = RA + H (I, J) * H (J, NA) 
                  SA = SA + H (I, J) * H (J, EN) 
   75          END DO 
               IF (VALI (I) .LT.ZERO) THEN 
                  Z = W 
                  R = RA 
                  S = SA 
               ELSE 
                  M = I 
                  IF (VALI (I) .EQ.ZERO) THEN 
                     CALL COMDIV ( - RA, - SA, W, Q, H (I, NA), H (I,   &
                     EN) )                                              
                  ELSE 
!                                                                       
!                    solve the complex linear system:                   
!           [ W+Q*I   X   ] [H(I,NA)+H(I,EN)*I    ]   [-RA-SA*I]        
!           [             ] [                     ] = [        ]        
!           [   Y   Z+Q*I ] [H(I+1,NA)+H(I+1,EN)*I]   [-R-S*I  ]        
!                                                                       
                     X = H (I, I + 1) 
                     Y = H (I + 1, I) 
                     VR = (VALR (I) - P) * (VALR (I) - P) + VALI (I)    &
                     * VALI (I) - Q * Q                                 
                     VI = TWO * Q * (VALR (I) - P) 
                     IF (VR.EQ.ZERO.AND.VI.EQ.ZERO) VR = EPS * NORM *   &
                     (ABS (W) + ABS (Q) + ABS (X) + ABS (Y) + ABS (Z) ) 
                     CALL COMDIV (X * R - Z * RA + Q * SA, X * S - Z *  &
                     SA - Q * RA, VR, VI, H (I, NA), H (I, EN) )        
                     IF (ABS (X) .GT.ABS (Z) + ABS (Q) ) THEN 
                        H (I + 1, NA) = ( - RA - W * H (I, NA) + Q * H (&
                        I, EN) ) / X                                    
                        H (I + 1, EN) = ( - SA - W * H (I, EN) - Q * H (&
                        I, NA) ) / X                                    
                     ELSE 
                        CALL COMDIV ( - R - Y * H (I, NA), - S - Y * H (&
                        I, EN), Z, Q, H (I + 1, NA), H (I + 1, EN) )    
                     ENDIF 
                  ENDIF 
               ENDIF 
  190       END DO 
         ENDIF 
  207 END DO 
!                                                                       
!     find eigenvectors for isolated eigenvalues:                       
!                                                                       
      DO 230 I = 1, N 
         IF (I.LT.LOW.OR.I.GT.HIGH) THEN 
            DO 220 J = I + 1, N 
  220       EIVEC (I, J) = H (I, J) 
         ENDIF 
  230 END DO 
!                                                                       
!     multiply by transformation matrix in order to                     
!     obtain eigenvectors of the original matrix MAT:                   
!                                                                       
      DO 300 J = N, LOW, - 1 
         IF (J.LE.HIGH) THEN 
            M = J 
         ELSE 
            M = HIGH 
         ENDIF 
         L = J - 1 
         IF (VALI (J) .LT.ZERO) THEN 
            DO 330 I = LOW, HIGH 
               Y = ZERO 
               Z = ZERO 
               DO 320 K = LOW, M 
                  Y = Y + EIVEC (I, K) * H (K, L) 
                  Z = Z + EIVEC (I, K) * H (K, J) 
  320          END DO 
               EIVEC (I, L) = Y 
               EIVEC (I, J) = Z 
  330       END DO 
         ELSE 
            IF (VALI (J) .EQ.ZERO) THEN 
               DO 350 I = LOW, HIGH 
                  Z = ZERO 
                  DO 340 K = LOW, M 
  340             Z = Z + EIVEC (I, K) * H (K, J) 
  350          EIVEC (I, J) = Z 
            ENDIF 
         ENDIF 
  300 END DO 
!                                                                       
!     Return  0: no error                                               
!                                                                       
      HQR2 = 0 
      END FUNCTION HQR2                             
!                                                                       
!                                                                       
      SUBROUTINE COMDIV (AR, AI, BR, BI, RESR, RESI) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     PROGRAM OBJECTIVE:                                         *      
!     ==================                                         *      
!     complex division: RESR+I*RESI := (AR+I*AI)/(BR+I*BI).      *      
!     (this procedure should not be called if BR=BI=0.)          *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     AR,AI:     real and imaginary parts of the numerator       *      
!     BR,BI:     real and imaginary parts of the denominator     *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     RESR,RESI: real and imaginary parts of the quotient        *      
!                                                                *      
!     LOCAL VARIABLES:                                           *      
!     ================                                           *      
!     ZERO:              floating-point constant 0               *      
!     TEMP1,TEMP2,TEMP3: auxiliary variables for intermediate    *      
!                        values                                  *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!                                                                *      
!  sources : Martin, R. S. and Wilkinson, J. H., see [MART68].   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Juergen Dietel                                     *      
!  date     : 04.10.1987                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION AR, AI, BR, BI, RESR, RESI 
      DOUBLEPRECISION ZERO 
      PARAMETER (ZERO = 0.0D0) 
      DOUBLEPRECISION TEMP1, TEMP2, TEMP3 
!                                                                       
      IF (BR.EQ.ZERO.AND.BI.EQ.ZERO) THEN 
         RESR = ZERO 
         RESI = ZERO 
         RETURN 
      ENDIF 
      IF (ABS (BR) .GT.ABS (BI) ) THEN 
         TEMP1 = BI / BR 
         TEMP2 = TEMP1 * BI + BR 
         TEMP3 = (AR + TEMP1 * AI) / TEMP2 
         RESI = (AI - TEMP1 * AR) / TEMP2 
         RESR = TEMP3 
      ELSE 
         TEMP1 = BR / BI 
         TEMP2 = TEMP1 * BR + BI 
         TEMP3 = (TEMP1 * AR + AI) / TEMP2 
         RESI = (TEMP1 * AI - AR) / TEMP2 
         RESR = TEMP3 
      ENDIF 
      END SUBROUTINE COMDIV                         
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION COMABS (AR, AI) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     PROGRAM OBJECTIVE:                                         *      
!     ==================                                         *      
!     Determine the absolute value of the complex number AR+I*AI:*      
!     COMABS:=SQRT(AR*AR+AI*AI)                                  *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     AR,AI: the real and imaginary parts of the complex number  *      
!            whose absolute value is to be determined            *      
!                                                                *      
!     OUTPUT PARAMETER:                                          *      
!     =================                                          *      
!     none                                                       *      
!                                                                *      
!     RETURN VALUE:                                              *      
!     =============                                              *      
!     value of the complex parameter                             *      
!                                                                *      
!     LOCAL VARIABLES:                                           *      
!     ================                                           *      
!     ZERO,ONE:    constants                                     *      
!     TEMP1,TEMP2: auxiliary variables                           *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: SWAP                                    *      
!                                                                *      
!                                                                *      
!  sources : Martin, R. S. and Wilkinson, J. H., see [MART68].   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Juergen Dietel                                     *      
!  date     : 04.10.1987                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION AR, AI 
      DOUBLEPRECISION ZERO, ONE 
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0) 
      DOUBLEPRECISION TEMP1, TEMP2 
!                                                                       
      TEMP1 = ABS (AR) 
      TEMP2 = ABS (AI) 
      IF (AR.EQ.ZERO.AND.AI.EQ.ZERO) THEN 
         COMABS = ZERO 
         RETURN 
      ENDIF 
      IF (TEMP2.GT.TEMP1) CALL SWAP (TEMP1, TEMP2) 
      IF (TEMP2.EQ.ZERO) THEN 
         COMABS = TEMP1 
      ELSE 
         COMABS = TEMP1 * SQRT (ONE+ (TEMP2 / TEMP1) **2) 
      ENDIF 
      END FUNCTION COMABS                           
!                                                                       
!                                                                       
      INTEGER FUNCTION NORMAL (LD, N, V, WI) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     PROGRAM OBJECTIVE:                                         *      
!     ==================                                         *      
!     NORMAL normalizes a set of eigenvectors in the maximum norm*      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     LD:     leading dimension of array V as defined in the     *      
!             calling program                                    *      
!     N:      the order of the square array V                    *      
!     V:      (N,N) matrix of type DOUBLE PRECISION that contains*      
!             the eigenvectors column-wise (see its description  *      
!             in EIGEN under EIVEC)                              *      
!     WI:     N-vector WI(1:N) of type DOUBLE PRECISION. Its     *      
!             components are the imaginary parts of the          *      
!             corresponding eigenvalues                          *      
!                                                                *      
!     OUTPUT PARAMETER:                                          *      
!     =================                                          *      
!     V:      array containing the normalized eigenvectors       *      
!                                                                *      
!     LOCAL VARIABLES:                                           *      
!     ================                                           *      
!     ZERO,ONE: floating-point constants 0 and 1                 *      
!     I,J:      indexing variables                               *      
!     MAXI:     auxiliary variable for determining the norm of   *      
!               a real vector                                    *      
!     TR,TI:    auxiliary variables for determining the norm of  *      
!               a complex vector                                 *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: COMABS, COMDIV                          *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Juergen Dietel                                     *      
!  date     : 04.10.1987                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER LD, N 
      DOUBLEPRECISION V (1:LD, 1:N), WI (1:N) 
      DOUBLEPRECISION ZERO, ONE 
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0) 
      INTEGER I, J 
      DOUBLEPRECISION MAXI, TR, TI, COMABS 
!                                                                       
      J = 1 
   10 IF (J.GT.N) GOTO 80 
      IF (WI (J) .EQ.ZERO) THEN 
         MAXI = V (1, J) 
         DO 15 I = 2, N 
   15    IF (ABS (V (I, J) ) .GT.ABS (MAXI) ) MAXI = V (I, J) 
         IF (MAXI.NE.ZERO) THEN 
            MAXI = ONE / MAXI 
            DO 20 I = 1, N 
   20       V (I, J) = V (I, J) * MAXI 
         ENDIF 
         J = J + 1 
      ELSE 
         TR = V (1, J) 
         TI = V (1, J + 1) 
         DO 30 I = 2, N 
            IF (COMABS (V (I, J), V (I, J + 1) ) .GT.COMABS (TR, TI) )  &
            THEN                                                        
               TR = V (I, J) 
               TI = V (I, J + 1) 
            ENDIF 
   30    END DO 
         IF (TR.NE.ZERO.OR.TI.NE.ZERO) THEN 
            DO 40 I = 1, N 
   40       CALL COMDIV (V (I, J), V (I, J + 1), TR, TI, V (I, J),      &
            V (I, J + 1) )                                              
         ENDIF 
         J = J + 2 
      ENDIF 
      GOTO 10 
   80 NORMAL = 0 
      END FUNCTION NORMAL                           
!                                                                       
!                                                                       
      SUBROUTINE SWAP (X, Y) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     PROGRAM OBJECTIVE:                                         *      
!     ==================                                         *      
!     SWAP interchanges the values of the two DOUBLE PRECISION   *      
!     variables X and Y.                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION X, Y 
      DOUBLEPRECISION TEMP 
      TEMP = X 
      X = Y 
      Y = TEMP 
      END SUBROUTINE SWAP                           
