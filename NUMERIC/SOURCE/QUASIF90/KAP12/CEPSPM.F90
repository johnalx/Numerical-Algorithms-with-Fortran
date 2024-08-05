      SUBROUTINE CEPSPM (AP, N, B, KPVT, RCOND, Z, WK) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Condition estimate for a symmetric matrix AP which is given*      
!     in factored and condensed form from SUBROUTINE ZSPMMK.     *      
!     If a condition estimate is not required, the SUBROUTINE    *      
!     ZSPMOK is more time efficient.                             *      
!     In order to solve A*X = B, a subsequent call of SESSPM     *      
!     is necessary with WK as input matrix.                      *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     AP    DOUBLE PRECISION vector AP(1:N*(N+1)/2), that        *      
!           contains the symmetric matrix A in condensed form.   *      
!           The columns of the upper triangle of A are stored in *      
!           sequence in the vector AP as subvectors of length    *      
!           N*(N+1)/2                                            *      
!     N     dimension of the matrix A                            *      
!     B     DOUBLE PRECISION vector B(1:N) containing the        *      
!           right-hand side of the system of equations A*X = B.  *      
!           The right-hand side is required for the SUBROUTINE   *      
!           ZSPMMK, since the decomposition and possible alter-  *      
!           ation of the right-hand side is performed in this    *      
!           subroutine in order to save arithmetic operations.   *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     AP    DOUBLE PRECISION vector AP(1:N*(N+1)/2), representing*      
!           a block diagonal matrix that contains the factors of *      
!           a decomposition in condensed form. The decomposition *      
!           is given as  A = U*D*TRANS(U). Here U is the product *      
!           of the permutation matrix and a unit upper triangular*      
!           matrix, TRANS(U) denotes the transpose of U and D is *      
!           a block diagonal matrix composed of 1 x 1 and 2 x 2  *      
!           blocks                                               *      
!     KPVT  INTEGER vector KPVT(1:N) containing the PIVOT indices*      
!     RCOND DOUBLE PRECISION estimate of the reciprocal condi-   *      
!           tion number of A.  For the system A*X = B, with      *      
!           relative input errors of size EPSILON in A and B,    *      
!           the relative error in the solution X will have the   *      
!           size EPSILON/RCOND.                                  *      
!           If RCOND is smaller than the machine constant, then  *      
!           the matrix A is numerically singular.                *      
!     Z     DOUBLE PRECISION auxiliary vector Z(1:N). Usually    *      
!           the contents of Z has no significance. In case A is  *      
!           close to being singular then Z is an approximated    *      
!           null-vector for A, i.e.,                             *      
!                   NORM(A*Z) = RCOND*NORM(A)*NORM(Z).           *      
!     B     DOUBLE PRECISION vector B(1:N) containing the right- *      
!           side of the system of equations  A*X = B  for use in *      
!           the SUBROUTINE SESSPM in the necessary form.         *      
!     WK    DOUBLE PRECISION vector WK(1:N*(N+1)/2)  (input      *      
!           parameter for SESSPM), auxiliary vector containing   *      
!           entries that are required for solving the system of  *      
!           equations  A*X = B. This way the solution of  A*X = B*      
!           can be determined in one, instead of two elimination *      
!           steps by the SUBROUTINE SESSPM.                      *      
!                                                                *      
!                                                                *      
!     condensed form                                             *      
!                                                                *      
!           The following code condenses the upper triangular    *      
!           part of a symmetric matrix A to a vector AP          *      
!                                                                *      
!                       K = 0                                    *      
!                       DO 20 J = 1, N                           *      
!                          DO 10 I= 1, N                         *      
!                             K = K + 1                          *      
!                             AP(K) = A(I,J)                     *      
!                    10    CONTINUE                              *      
!                    20 CONTINUE                                 *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!     subroutines required: ZSPMMK, PCOSOL, PCOLTG, VECADD,      *      
!                           VECMWC, SCAPRO, ABSSUM, INDMAX,      *      
!                           VECXCH                               *      
!                                                                *      
!                                                                *      
!     source : Linpack User's Guide, SIAM, Philadelphia, 1979    *      
!                                                                *      
!              The source was converted to FORTRAN 77. In some   *      
!              instances it had to be modified and adjusted for  *      
!              the requirements of our specific calling programs.*      
!              This program and the corresponding subprograms    *      
!              are not compatible with the originals from the    *      
!              Linpack User's Guide.                             *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!     authors :  Michael Groenheim, Ina Hinze                    *      
!     date    :  10.25.1989                                      *      
!     source  :  FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER N, KPVT (N) 
      DOUBLEPRECISION AP (N * (N + 1) / 2), Z (N), B (N), WK (N *       &
      (N + 1) / 2)                                                      
      DOUBLEPRECISION RCOND 
      DOUBLEPRECISION EK, ANORM, S, ABSSUM, YNORM 
      INTEGER I, IJ, IERR, J, JM1, J1 
!                                                                       
!     determine the norm of A making use of symmetry                    
!                                                                       
      J1 = 1 
      DO 20 J = 1, N 
         Z (J) = ABSSUM (J, AP (J1) ) 
         IJ = J1 
         J1 = J1 + J 
         JM1 = J - 1 
         IF (JM1.GE.1) THEN 
            DO 10 I = 1, JM1 
               Z (I) = Z (I) + DABS (AP (IJ) ) 
               IJ = IJ + 1 
   10       END DO 
         ENDIF 
   20 END DO 
      ANORM = 0.0D0 
      DO 30 J = 1, N 
         ANORM = DMAX1 (ANORM, Z (J) ) 
   30 END DO 
      IF (ANORM.EQ.0.0D0) THEN 
         RCOND = 0.0D0 
         RETURN 
      ENDIF 
!                                                                       
!     factor the matrix by using SUBROUTINE ZSPMMK                      
!                                                                       
      CALL ZSPMMK (AP, N, B, KPVT, WK, IERR) 
!                                                                       
!     RCOND = 1/(NORM(A)*(ESTIMATE OF THE NORM(INVERSE(A)))).           
!     ESTIMATE = NORM(Z)/NORM(Y). Here A*Z = Y and A*Y = E.             
!     The elements of E are chosen in such a way that the elements      
!     of W become as large as possible, where U*D*W = E.                
!                                                                       
      EK = 1.0D0 
      DO 40 J = 1, N 
         Z (J) = 0.0D0 
   40 END DO 
!                                                                       
!     solve U*D*W = E                                                   
!                                                                       
      CALL PCOSOL (1, AP, N, Z, KPVT, S, EK, YNORM) 
!                                                                       
!     solve TRANS(U) * Y = W                                            
!                                                                       
      CALL PCOLTG (AP, N, Z, KPVT, S) 
      S = 1.0D0 / ABSSUM (N, Z) 
      CALL VECMWC (N, S, Z) 
!                                                                       
      YNORM = 1.0D0 
!                                                                       
!     solve U*D*V = Y                                                   
!                                                                       
      CALL PCOSOL (2, AP, N, Z, KPVT, S, EK, YNORM) 
!                                                                       
!     solve TRANS(U) * Z = V                                            
!                                                                       
      CALL PCOLTG (AP, N, Z, KPVT, S) 
      YNORM = S * YNORM 
!                                                                       
!     set ZNORM = 1.0                                                   
!                                                                       
      S = 1.0D0 / ABSSUM (N, Z) 
      CALL VECMWC (N, S, Z) 
      YNORM = S * YNORM 
!                                                                       
      RCOND = YNORM / ANORM 
      RETURN 
      END SUBROUTINE CEPSPM                         
!                                                                       
!                                                                       
      SUBROUTINE ZSPMMK (AP, N, B, KPVT, WK, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Factorization of a symmetric matrix, given in condensed    *      
!     format, by elimination that relies on symmetric pivoting.  *      
!     If a condition estimate is not required, the SUBROUTINE    *      
!     ZSPMOK will perform the task faster.                       *      
!     In order to solve A*X = B the subsequent call of SESSPM is *      
!     necessary.                                                 *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     AP    DOUBLE PRECISION vector AP(1:(N*(N+1)/2)), contain-  *      
!           ing the symmetric matrix A in condensed form.        *      
!           The columns of its upper triangle are stored sequen- *      
!           tially in the vector AP of length N*(N+1)/2          *      
!     N     dimension of the matrix A                            *      
!     B     DOUBLE PRECISION vector B(1:N), the right-hand side. *      
!           The right-hand side is required in SUBROUTINE ZSPMMK *      
!           in order to save arithmetic operations.              *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     AP    DOUBLE PRECISION vector AP(1:N*(N+1)/2), representing*      
!           a block diagonal matrix that contains the factors of *      
!           a decomposition in condensed form. The decomposition *      
!           is given as  A = U*D*TRANS(U). Here U is the product *      
!           of the permutation matrix and a unit upper triangular*      
!           matrix, TRANS(U) denotes the transpose of U and D is *      
!           a block diagonal matrix composed of 1 x 1 and 2 x 2  *      
!           blocks                                               *      
!     KPVT  INTEGER vector KPVT(1:N) containing the PIVOT indices*      
!     B     DOUBLE PRECISION vector B(1:N) containing the right- *      
!           side of the system of equations  A*X = B  for use in *      
!           the SUBROUTINE SESSPM in the necessary form.         *      
!     WK    DOUBLE PRECISION vector WK(1:N*(N+1)/2)  (input      *      
!           parameter for SESSPM), auxiliary vector containing   *      
!           entries that are required for solving the system of  *      
!           equations  A*X = B. This way the solution of  A*X = B*      
!           can be determined in one, instead of two elimination *      
!           steps by the SUBROUTINE SESSPM.                      *      
!     IERR  error parameter                                      *      
!           = 0, everything is o.k                               *      
!           = K, the K-th PIVOT block is singular.               *      
!                For the current subroutine, this does not denote*      
!                an error; however, it indicates that the        *      
!                SUBROUTINE SESSPM may encounter division by zero*      
!                                                                *      
!                                                                *      
!     condensed form                                             *      
!                                                                *      
!           The following code condenses the upper triangular    *      
!           part of a symmetric matrix A to a vector AP          *      
!                                                                *      
!                       K = 0                                    *      
!                       DO 20 J = 1, N                           *      
!                          DO 10 I= 1, N                         *      
!                             K = K + 1                          *      
!                             AP(K) = A(I,J)                     *      
!                    10    CONTINUE                              *      
!                    20 CONTINUE                                 *      
!                    20 CONTINUE                                 *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!     subroutines required: VECADD, VECXCH, INDMAX               *      
!                                                                *      
!                                                                *      
!     source : Linpack User's Guide, SIAM, Philadelphia, 1979    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!     authors :  Michael Groenheim, Ina Hinze                    *      
!     date    :  10.25.1989                                      *      
!     source  :  FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER N, KPVT (N), IERR 
      DOUBLEPRECISION AP (N * (N + 1) / 2), B (N), WK (N * (N + 1)      &
      / 2)                                                              
      DOUBLEPRECISION D, D1, D2, T 
      DOUBLEPRECISION ABSAKK, ALPHA, COLMAX, ROWMAX 
      INTEGER INDMAX, IJ, IJJ, IK, IKM1, IM, IMAX, IMAXP1, IMIM, IMJ,   &
      IMK                                                               
      INTEGER J, JJ, JK, JKM1, JMAX, JMIM, K, KK, KM1, KM1K, KM1KM1,    &
      KM2, KSTEP                                                        
      LOGICAL SWAP 
!                                                                       
!     store condensed vector AP in auxiliary vector WK                  
!                                                                       
      DO 5 I = 1, N * (N + 1) / 2 
         WK (I) = AP (I) 
    5 END DO 
!                                                                       
!     initialization                                                    
!                                                                       
!     ALPHA is used to determine the PIVOT block size                   
!                                                                       
      ALPHA = (1.0D0 + SQRT (17.0D0) ) / 8.0D0 
!                                                                       
      IERR = 0 
!                                                                       
!     main loop over K; K runs from N back to 1                         
!                                                                       
      K = N 
      IK = (N * (N - 1) ) / 2 
   10 CONTINUE 
!                                                                       
!     leave the loop if K=0 or K=1                                      
!                                                                       
      IF (K.EQ.0) RETURN 
      IF (K.LE.1) THEN 
         KPVT (1) = 1 
         IF (AP (1) .EQ.0.0D0) IERR = 1 
         RETURN 
      ENDIF 
!                                                                       
!     this part of the program determines the elimination method to be  
!     used. After this part has been executed, KSTEP is set to be the si
!     of the PIVOT block and SWAP is set to .TRUE., if swapping is      
!     necessary.                                                        
!                                                                       
      KM1 = K - 1 
      KK = IK + K 
      ABSAKK = DABS (AP (KK) ) 
!                                                                       
!     determine the largest off-diagonal element in                     
!     magnitude in column K                                             
!                                                                       
      IMAX = INDMAX (K - 1, AP (IK + 1) ) 
      IMK = IK + IMAX 
      COLMAX = DABS (AP (IMK) ) 
      IF (ABSAKK.LT.ALPHA * COLMAX) THEN 
!                                                                       
!        determine the largest off-diagonal element                     
!        in magnitude in row IMAX                                       
!                                                                       
         ROWMAX = 0.0D0 
         IMAXP1 = IMAX + 1 
         IM = IMAX * (IMAX - 1) / 2 
         IMJ = IM + 2 * IMAX 
         DO 20 J = IMAXP1, K 
            ROWMAX = DMAX1 (ROWMAX, DABS (AP (IMJ) ) ) 
            IMJ = IMJ + J 
   20    END DO 
         IF (IMAX.NE.1) THEN 
            JMAX = INDMAX (IMAX - 1, AP (IM + 1) ) 
            JMIM = JMAX + IM 
            ROWMAX = DMAX1 (ROWMAX, DABS (AP (JMIM) ) ) 
         ENDIF 
         IMIM = IMAX + IM 
         IF (DABS (AP (IMIM) ) .LT.ALPHA * ROWMAX) THEN 
            IF (ABSAKK.LT.ALPHA * COLMAX * (COLMAX / ROWMAX) ) THEN 
               KSTEP = 2 
               SWAP = IMAX.NE.KM1 
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
      IF (DMAX1 (ABSAKK, COLMAX) .EQ.0.0D0) THEN 
!                                                                       
!        column K is the zero column. Modify IERR and reiterate loop    
!                                                                       
         KPVT (K) = K 
         IERR = K 
         IK = IK - (K - 1) 
         IF (KSTEP.EQ.2) IK = IK - (K - 2) 
         K = K - KSTEP 
         GOTO 10 
      ENDIF 
      IF (KSTEP.EQ.2) THEN 
!                                                                       
!        2 x 2 PIVOT block                                              
!                                                                       
         KM1K = IK + K - 1 
         IKM1 = IK - (K - 1) 
         IF (SWAP) THEN 
!                                                                       
!           perform swap                                                
!                                                                       
            CALL VECXCH (IMAX, AP (IM + 1), AP (IKM1 + 1) ) 
            CALL VECXCH (IMAX, WK (IM + 1), WK (IKM1 + 1) ) 
            IMJ = IKM1 + IMAX 
            DO 30 JJ = IMAX, KM1 
               J = KM1 + IMAX - JJ 
               JKM1 = IKM1 + J 
               T = AP (JKM1) 
               T1 = WK (JKM1) 
               AP (JKM1) = AP (IMJ) 
               WK (JKM1) = WK (IMJ) 
               AP (IMJ) = T 
               WK (IMJ) = T1 
               IMJ = IMJ - (J - 1) 
   30       END DO 
            T = AP (KM1K) 
            AP (KM1K) = AP (IMK) 
            AP (IMK) = T 
            T = WK (KM1K) 
            WK (KM1K) = WK (IMK) 
            WK (IMK) = T 
            T = B (K - 1) 
            B (K - 1) = B (IMAX) 
            B (IMAX) = T 
         ENDIF 
                                                                        
!                                                                       
!        perform elimination                                            
!                                                                       
         KM2 = K - 2 
         IF (KM2.NE.0) THEN 
            KM1KM1 = IKM1 + K - 1 
            D = AP (KM1K) * AP (KM1K) - AP (KK) * AP (KM1KM1) 
            IJ = IK - (K - 1) - (K - 2) 
            DO 40 JJ = 1, KM2 
               J = KM1 - JJ 
               JK = IK + J 
               JKM1 = IKM1 + J 
               D1 = (AP (KM1KM1) * AP (JK) - AP (JKM1) * AP (KM1K) )    &
               / D                                                      
               D2 = (AP (KK) * AP (JKM1) - AP (JK) * AP (KM1K) )        &
               / D                                                      
               CALL VECADD (J, D1, AP (IK + 1), AP (IJ + 1) ) 
               CALL VECADD (1, D1, B (K), B (J) ) 
               CALL VECADD (J, D2, AP (IKM1 + 1), AP (IJ + 1) ) 
               CALL VECADD (1, D2, B (K - 1), B (J) ) 
               IF (IKM1.EQ.1) WK (IKM1) = AP (IKM1) 
               DO 45 ID = 1, J 
                  WK (IJ + ID) = AP (IJ + ID) 
   45          END DO 
               AP (JK) = D1 
               AP (JKM1) = D2 
               IJJ = IJ + J 
               IJ = IJ - (J - 1) 
   40       END DO 
         ENDIF 
!                                                                       
!        set up PIVOT vector                                            
!                                                                       
         KPVT (K) = 1 - K 
         IF (SWAP) KPVT (K) = - IMAX 
         KPVT (K - 1) = KPVT (K) 
      ELSE 
!                                                                       
!        1 x 1 PIVOT block                                              
!                                                                       
         IF (SWAP) THEN 
!                                                                       
!           perform swap                                                
!                                                                       
            CALL VECXCH (IMAX, AP (IM + 1), AP (IK + 1) ) 
            CALL VECXCH (IMAX, WK (IM + 1), WK (IK + 1) ) 
            IMJ = IK + IMAX 
            DO 50 JJ = IMAX, K 
               J = K + IMAX - JJ 
               JK = IK + J 
               T = AP (JK) 
               AP (JK) = AP (IMJ) 
               AP (IMJ) = T 
               T = WK (JK) 
               WK (JK) = WK (IMJ) 
               WK (IMJ) = T 
               IMJ = IMJ - (J - 1) 
   50       END DO 
            T = B (K) 
            B (K) = B (IMAX) 
            B (IMAX) = T 
         ENDIF 
!                                                                       
!        perform elimination                                            
!                                                                       
         IJ = IK - (K - 1) 
         DO 60 JJ = 1, KM1 
            J = K - JJ 
            JK = IK + J 
            D1 = - AP (JK) / AP (KK) 
            CALL VECADD (J, D1, AP (IK + 1), AP (IJ + 1) ) 
            CALL VECADD (1, D1, B (K), B (J) ) 
            DO 70 ID = 1, J 
               WK (IJ + ID) = AP (IJ + ID) 
   70       END DO 
            IJJ = IJ + J 
            AP (JK) = D1 
            IJ = IJ - (J - 1) 
   60    END DO 
!                                                                       
!        set up PIVOT vector                                            
!                                                                       
         KPVT (K) = K 
         IF (SWAP) KPVT (K) = IMAX 
      ENDIF 
      IK = IK - (K - 1) 
      IF (KSTEP.EQ.2) IK = IK - (K - 2) 
      K = K - KSTEP 
      GOTO 10 
      END SUBROUTINE ZSPMMK                         
!                                                                       
!                                                                       
      SUBROUTINE VECADD (N, SA, SX, SY) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Multiplies the vector SX by the constant SA and then adds  *      
!     the vector SX to the vector SY.                            *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!     subroutines required: none                                 *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!     authors :  Michael Groenheim, Ina Hinze                    *      
!     date    :  10.25.1989                                      *      
!     source  :  FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION SX (N), SY (N), SA 
      INTEGER I, N 
!                                                                       
      IF (N.LE.0) RETURN 
!                                                                       
      IF (SA.EQ.0.0D0) RETURN 
!                                                                       
      DO 10 I = 1, N 
         SY (I) = SY (I) + SA * SX (I) 
   10 END DO 
      RETURN 
      END SUBROUTINE VECADD                         
!                                                                       
!                                                                       
      SUBROUTINE VECXCH (N, SX, SY) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Swaps the vectors SX and SY.                               *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!     subroutines required: none                                 *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!     authors :  Michael Groenheim, Ina Hinze                    *      
!     date    :  10.25.1989                                      *      
!     source  :  FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION SX (N), SY (N), DUMMY 
      INTEGER I, N 
!                                                                       
      DO 10 I = 1, N 
         DUMMY = SX (I) 
         SX (I) = SY (I) 
         SY (I) = DUMMY 
   10 END DO 
      RETURN 
      END SUBROUTINE VECXCH                         
!                                                                       
!                                                                       
      INTEGER FUNCTION INDMAX (N, SX) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Determines the index of the element of the vector SX that  *      
!     has the largest entry in magnitude.                        *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!     subroutines required: none                                 *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!     authors :  Michael Groenheim, Ina Hinze                    *      
!     date    :  10.25.1989                                      *      
!     source  :  FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION SX (N), SMAX 
      INTEGER I, N 
!                                                                       
      INDMAX = 0 
!                                                                       
      IF (N.LT.1) RETURN 
!                                                                       
      INDMAX = 1 
!                                                                       
      IF (N.EQ.1) RETURN 
!                                                                       
      SMAX = DABS (SX (1) ) 
!                                                                       
      DO 10 I = 2, N 
         IF (DABS (SX (I) ) .GT.SMAX) THEN 
            INDMAX = I 
            SMAX = DABS (SX (I) ) 
         ENDIF 
   10 END DO 
      RETURN 
      END FUNCTION INDMAX                           
!                                                                       
!                                                                       
      SUBROUTINE PCOSOL (IPOINT, AP, N, Z, KPVT, S, EK, YNORM) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Subroutine of CEPSPM.                                      *      
!     This subroutine helps to solve a system of equations       *      
!                   U*D*W  = E  or  U*D*V = Y.                   *      
!     It is called twice by ZSPM.. . The method of solving the   *      
!     equations above is identical except for small variations.  *      
!     In order to combine both tasks in one subroutine, we       *      
!     introduce the parameter IPOINT, which, depending on the    *      
!     kind of call, initiates branches at the relevant places in *      
!     this subroutine.                                           *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     IPOINT  flag                                               *      
!             = 1  , in case PCOSOL is called for solving        *      
!                    U*D*W = E                                   *      
!             = 2  , in case PCOSOL is called for solving        *      
!                    U*D*V = Y                                   *      
!     N       dimension of the matrix A                          *      
!     KPVT    INTEGER vector KVPT(1:N) containing the PIVOT      *      
!             indices                                            *      
!     Z       DOUBLE PRECISION auxiliary vector Z(1:N), needed   *      
!             for solving  U*D*W = E  or  U*D*V = Y.             *      
!     AP      vector containing the factors of the decomposition *      
!             of the symmetric matrix A in condensed form.       *      
!     S       DOUBLE PRECISION auxiliary variable                *      
!     EK      DOUBLE PRECISION auxiliary variable                *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     S       DOUBLE PRECISION auxiliary variable                *      
!     YNORM   DOUBLE PRECISION norm of Y                         *      
!     EK      DOUBLE PRECISION auxiliary variable                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!     subroutines required: VECADD, VECMWC                       *      
!                                                                *      
!                                                                *      
!     source : Linpack User's Guide, SIAM, Philadelphia, 1979    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!     authors :  Michael Groenheim, Ina Hinze                    *      
!     date    :  10.25.1989                                      *      
!     source  :  FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER N, IPOINT, KPVT (N) 
      DOUBLEPRECISION AP (N * (N + 1) / 2), Z (N) 
      DOUBLEPRECISION D, D1, D2, EK, T 
      DOUBLEPRECISION S, YNORM 
      INTEGER IK, IKM1, K, KK, KM1K, KM1KM1, KP, KPS, KSTEP 
      K = N 
      IK = N * (N - 1) / 2 
   10 IF (K.EQ.0) RETURN 
      KK = IK + K 
      IKM1 = IK - (K - 1) 
      KSTEP = 1 
      IF (KPVT (K) .LT.0) KSTEP = 2 
      IF (IPOINT.EQ.1) THEN 
!                                                                       
!        call was for solving U*D*W = E                                 
!                                                                       
         KP = IABS (KPVT (K) ) 
         KPS = K + 1 - KSTEP 
         IF (KP.NE.KPS) THEN 
            T = Z (KPS) 
            Z (KPS) = Z (KP) 
            Z (KP) = T 
         ENDIF 
         IF (Z (K) .NE.0.0D0) EK = DSIGN (EK, Z (K) ) 
         Z (K) = Z (K) + EK 
         CALL VECADD (K - KSTEP, Z (K), AP (IK + 1), Z (1) ) 
      ENDIF 
      IF (IPOINT.EQ.2) THEN 
!                                                                       
!        call was for solving U*D*V = Y                                 
!                                                                       
         IF (K.NE.KSTEP) THEN 
            KP = IABS (KPVT (K) ) 
            KPS = K + 1 - KSTEP 
            IF (KP.NE.KPS) THEN 
               T = Z (KPS) 
               Z (KPS) = Z (KP) 
               Z (KP) = T 
            ENDIF 
            CALL VECADD (K - KSTEP, Z (K), AP (IK + 1), Z (1) ) 
            IF (KSTEP.EQ.2) THEN 
               CALL VECADD (K - KSTEP, Z (K - 1), AP (IKM1 + 1),        &
               Z (1) )                                                  
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
!     1 x 1 PIVOT block                                                 
!                                                                       
      IF (KSTEP.EQ.1) THEN 
         IF (DABS (Z (K) ) .GT.DABS (AP (KK) ) ) THEN 
            S = DABS (AP (KK) ) / DABS (Z (K) ) 
            CALL VECMWC (N, S, Z) 
            IF (IPOINT.EQ.1) EK = S * EK 
            IF (IPOINT.EQ.2) YNORM = S * YNORM 
         ENDIF 
         IF (AP (KK) .NE.0.0D0) Z (K) = Z (K) / AP (KK) 
         IF (AP (KK) .EQ.0.0D0) Z (K) = 1.0D0 
         K = K - KSTEP 
         IK = IK - K 
         GOTO 10 
      ENDIF 
!                                                                       
!     2 x 2 PIVOT block                                                 
!                                                                       
      IF (KSTEP.EQ.2) THEN 
         IF (IPOINT.EQ.1) THEN 
            IF (Z (K - 1) .NE.0.0D0) EK = DSIGN (EK, Z (K - 1) ) 
            Z (K - 1) = Z (K - 1) + EK 
            CALL VECADD (K - KSTEP, Z (K - 1), AP (IKM1 + 1), Z (1) ) 
         ENDIF 
         KM1K = IK + K - 1 
         KM1KM1 = IKM1 + K - 1 
         D1 = AP (KM1KM1) * Z (K) - AP (KM1K) * Z (K - 1) 
         D2 = AP (KK) * Z (K - 1) - AP (KM1K) * Z (K) 
         D = AP (KK) * AP (KM1KM1) - AP (KM1K) * AP (KM1K) 
         Z (K) = D1 / D 
         Z (K - 1) = D2 / D 
         K = K - KSTEP 
         IK = IK - K - (K + 1) 
         GOTO 10 
      ENDIF 
      END SUBROUTINE PCOSOL                         
!                                                                       
!                                                                       
      SUBROUTINE VECMWC (N, SA, SX) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Multiplies the vector SX by the constant SA.               *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!     subroutines required: none                                 *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!     authors :  Michael Groenheim, Ina Hinze                    *      
!     date    :  10.25.1989                                      *      
!     source  :  FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION SA, SX (N) 
      INTEGER I, N 
!                                                                       
      IF (N.LE.0) RETURN 
!                                                                       
      DO 10 I = 1, N 
         SX (I) = SA * SX (I) 
   10 END DO 
      RETURN 
      END SUBROUTINE VECMWC                         
!                                                                       
!                                                                       
      SUBROUTINE PCOLTG (AP, N, Z, KPVT, S) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Subroutine of CEPSPM.                                      *      
!     This subroutine helps to solve the system the equations    *      
!     TRANS(U) * Y = W   or  TRANS(U) * Z = V.                   *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     N     dimension of the matrix A                            *      
!     Z     DOUBLE PRECISION auxiliary vector Z(1:N), needed     *      
!           to solve TRANS(U) * Y = W  or  TRANS(U) * Z = V      *      
!     KPVT  INTEGER vector KPVT(1:N) containing the PIVOT        *      
!           indices                                              *      
!     AP    vector containing the factors of the decomposition   *      
!           of the symmetric matrix A in condensed form          *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETER:                                          *      
!     =================                                          *      
!     S     DOUBLE PRECISION auxiliary variable                  *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!     subroutines required: SCAPRO, VECMWC, ABSSUM               *      
!                                                                *      
!                                                                *      
!     source : Linpack User's Guide, SIAM, Philadelphia, 1979    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!     authors :  Michael Groenheim, Ina Hinze                    *      
!     date    :  10.25.1989                                      *      
!     source  :  FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER N, KPVT (N) 
      DOUBLEPRECISION AP (N * (N + 1) / 2), Z (N) 
      DOUBLEPRECISION SCAPRO, T, S, ABSSUM 
      INTEGER IK, IKP1, K, KP, KSTEP 
      S = 1.0D0 / ABSSUM (N, Z) 
      CALL VECMWC (N, S, Z) 
      K = 1 
      IK = 0 
   10 IF (K.LE.N) THEN 
         KSTEP = 1 
         IF (KPVT (K) .LT.0) KSTEP = 2 
         IF (K.NE.1) THEN 
            Z (K) = Z (K) + SCAPRO (K - 1, AP (IK + 1), Z (1) ) 
            IKP1 = IK + K 
            IF (KSTEP.EQ.2) Z (K + 1) = Z (K + 1) + SCAPRO (K - 1, AP ( &
            IKP1 + 1), Z (1) )                                          
            KP = IABS (KPVT (K) ) 
            IF (KP.NE.K) THEN 
               T = Z (K) 
               Z (K) = Z (KP) 
               Z (KP) = T 
            ENDIF 
         ENDIF 
         IK = IK + K 
         IF (KSTEP.EQ.2) IK = IK + (K + 1) 
         K = K + KSTEP 
         GOTO 10 
      ENDIF 
      RETURN 
      END SUBROUTINE PCOLTG                         
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION SCAPRO (N, SX, SY) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Determines the scalar product of two vectors SX and SY.    *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!     subroutines required: none                                 *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!     authors :  Michael Groenheim, Ina Hinze                    *      
!     date    :  10.25.1989                                      *      
!     source  :  FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION SX (N), SY (N) 
      INTEGER I, N 
!                                                                       
      SCAPRO = 0.0D0 
!                                                                       
      DO 10 I = 1, N 
         SCAPRO = SCAPRO + SX (I) * SY (I) 
   10 END DO 
      RETURN 
      END FUNCTION SCAPRO                           
!                                                                       
!                                                                       
      DOUBLEPRECISION FUNCTION ABSSUM (N, SX) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Forms the sum of the absolute values of the entries in SX. *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!     subroutine required: none                                  *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!     authors :  Michael Groenheim, Ina Hinze                    *      
!     date    :  10.25.1989                                      *      
!     source  :  FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION SX (N) 
      INTEGER I, N 
!                                                                       
      ABSSUM = 0.0D0 
!                                                                       
      IF (N.LE.0) RETURN 
!                                                                       
      DO 10 I = 1, N 
         ABSSUM = ABSSUM + DABS (SX (I) ) 
   10 END DO 
      RETURN 
      END FUNCTION ABSSUM                           
