      SUBROUTINE ZSPMOK (AP, N, B, KPVT, INFO) 
!                                                                       
!*****************************************************************      
!                                                                *      
!     Decomposition of a symmetric matrix A given in condensed   *      
!     form using symmetric pivoting for elimination.             *      
!                                                                *      
!     To solve A*X = B one call of LGSSPM is required.           *      
!                                                                *      
!                                                                *      
!     INPUT PARAMETERS:                                          *      
!     =================                                          *      
!     AP    DOUBLE PRECISION vector AP(1:(N*(N+1)/2)), containing*      
!           the symmetric matrix A in condensed form.            *      
!           The columns of its upper triangle are stored cosecu- *      
!           tively in this vector of length N*(N+1)/2.           *      
!     N     Dimension of the matrix A                            *      
!     B     DOUBLE PRECISION vector B(1:N), the right-hand side B*      
!           of the linear system A*X = B.                        *      
!           The right-hand side is needed in SUBROUTINE ZSPMMK   *      
!           since the elimination and updating steps are per-    *      
!           formed there for both A and B in order to save time. *      
!                                                                *      
!                                                                *      
!     OUTPUT PARAMETERS:                                         *      
!     ==================                                         *      
!     AP    DOUBLE PRECISION vector AP(1:N*(N+1)/2) which holds  *      
!           the essential information needed to solve A*X = B.   *      
!           Hence A*X = B can be solved via SUBROUTINE LGSSPM in *      
!           one pass instead of two.                             *      
!     KPVT  INTEGER vector KPVT(1:N) of the PIVOT indices.       *      
!     B     DOUBLE PRECISION vector B(1:N), the right-hand side  *      
!           of the linear system A*X = B in a form suitable for  *      
!           applying SUBROUTINE LGSSPM.                          *      
!     INFO  error parameter                                      *      
!           = 0, all is o.k.                                     *      
!           = k, the K-th PIVOT block is numerically singular.   *      
!                This error message is irrelevant for this sub-  *      
!                routine. However, the SUBROUTINE LGSSPM might   *      
!                encounter division by zero.                     *      
!                                                                *      
!                                                                *      
!     condensed form                                             *      
!                                                                *      
!           The following code condenses the upper triangle of a *      
!           symmetric matrix A to a vector AP.                   *      
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
!     required subroutines: VECADD, VECXCH, INDMAX               *      
!                                                                *      
!                                                                *      
!     Source : Linpack User's Guide , SIAM, Philadelphia, 1979   *      
!                                                                *      
!              The sourcecode in FORTRAN 4 was translated into   *      
!              FORTRAN 77. Several details were modified and     *      
!              adapted to the needs of the calling program.      *      
!              Hence this and the related subroutines will not   *      
!              be compatible with the original codes in LINPACK. *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!     Authors   :  Michael Gr”nheim, Ina Hinze                   *      
!     Date      :  10.25.1989                                    *      
!     Source    :  FORTRAN 77                                    *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER N, KPVT ( * ), INFO 
      DOUBLEPRECISION AP ( * ), B ( * ) 
      DOUBLEPRECISION D, D1, D2, T 
      DOUBLEPRECISION ABSAKK, ALPHA, SPAMAX, ZEIMAX 
      INTEGER INDMAX, IJ, IJJ, IK, IKM1, IM, IMAX, IMAXP1, IMIM, IMJ,   &
      IMK                                                               
      INTEGER J, JJ, JK, JKM1, JMAX, JMIM, K, KK, KM1, KM1K, KM1KM1,    &
      KM2, KSTEP                                                        
      LOGICAL SWAP 
!                                                                       
!     Initializing                                                      
!                                                                       
!     ALPHA is used to determine the PIVOT block size                   
!                                                                       
      ALPHA = (1.0D0 + SQRT (17.0D0) ) / 8.0D0 
!                                                                       
      INFO = 0 
!                                                                       
!     Main loop over K, K from N to 1                                   
!                                                                       
      K = N 
      IK = (N * (N - 1) ) / 2 
   10 CONTINUE 
!                                                                       
!     leave this loop if K=0 or K=1                                     
!                                                                       
      IF (K.EQ.0) RETURN 
      IF (K.LE.1) THEN 
         KPVT (1) = 1 
         IF (AP (1) .EQ.0.0D0) INFO = 1 
         RETURN 
      ENDIF 
!                                                                       
!     This section determines the kind of elimination employed.         
!     At the end of this section, KSTEP is assigned the size of the     
!     PIVOT block and SWAP is set to .TRUE. in case a swap has taken pla
!                                                                       
      KM1 = K - 1 
      KK = IK + K 
      ABSAKK = DABS (AP (KK) ) 
!                                                                       
!     Compute the largest off-diagonal element in magnitude             
!     in column K                                                       
!                                                                       
      IMAX = INDMAX (K - 1, AP (IK + 1) ) 
      IMK = IK + IMAX 
      SPAMAX = DABS (AP (IMK) ) 
      IF (ABSAKK.LT.ALPHA * SPAMAX) THEN 
!                                                                       
!          Compute the largest off-diagonal element in magnitude        
!          in row IMAX                                                  
!                                                                       
         ZEIMAX = 0.0D0 
         IMAXP1 = IMAX + 1 
         IM = IMAX * (IMAX - 1) / 2 
         IMJ = IM + 2 * IMAX 
         DO 20 J = IMAXP1, K 
            ZEIMAX = DMAX1 (ZEIMAX, DABS (AP (IMJ) ) ) 
            IMJ = IMJ + J 
   20    END DO 
         IF (IMAX.NE.1) THEN 
            JMAX = INDMAX (IMAX - 1, AP (IM + 1) ) 
            JMIM = JMAX + IM 
            ZEIMAX = DMAX1 (ZEIMAX, DABS (AP (JMIM) ) ) 
         ENDIF 
         IMIM = IMAX + IM 
         IF (DABS (AP (IMIM) ) .LT.ALPHA * ZEIMAX) THEN 
            IF (ABSAKK.LT.ALPHA * SPAMAX * (SPAMAX / ZEIMAX) ) THEN 
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
      IF (DMAX1 (ABSAKK, SPAMAX) .EQ.0.0D0) THEN 
!                                                                       
!         Column  K is the zero vector. Record error in INFO            
!         and repeat loop                                               
!                                                                       
         KPVT (K) = K 
         INFO = K 
         IK = IK - (K - 1) 
         IF (KSTEP.EQ.2) IK = IK - (K - 2) 
         K = K - KSTEP 
         GOTO 10 
      ENDIF 
      IF (KSTEP.EQ.2) THEN 
!                                                                       
!         2 x 2 PIVOT block                                             
!                                                                       
         KM1K = IK + K - 1 
         IKM1 = IK - (K - 1) 
         IF (SWAP) THEN 
!                                                                       
!               swap                                                    
!                                                                       
            CALL VECXCH (IMAX, AP (IM + 1), AP (IKM1 + 1) ) 
            IMJ = IKM1 + IMAX 
            DO 30 JJ = IMAX, KM1 
               J = KM1 + IMAX - JJ 
               JKM1 = IKM1 + J 
               T = AP (JKM1) 
               AP (JKM1) = AP (IMJ) 
               AP (IMJ) = T 
               IMJ = IMJ - (J - 1) 
   30       END DO 
            T = AP (KM1K) 
            AP (KM1K) = AP (IMK) 
            AP (IMK) = T 
            T = B (K - 1) 
            B (K - 1) = B (IMAX) 
            B (IMAX) = T 
         ENDIF 
                                                                        
!                                                                       
!           Perform elimination                                         
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
               IJJ = IJ + J 
               IJ = IJ - (J - 1) 
   40       END DO 
         ENDIF 
!                                                                       
!           Set up PIVOT vector K                                       
!                                                                       
         KPVT (K) = 1 - K 
         IF (SWAP) KPVT (K) = - IMAX 
         KPVT (K - 1) = KPVT (K) 
      ELSE 
!                                                                       
!           1 x 1 PIVOT block                                           
!                                                                       
         IF (SWAP) THEN 
!                                                                       
!                swap                                                   
!                                                                       
            CALL VECXCH (IMAX, AP (IM + 1), AP (IK + 1) ) 
            IMJ = IK + IMAX 
            DO 50 JJ = IMAX, K 
               J = K + IMAX - JJ 
               JK = IK + J 
               T = AP (JK) 
               AP (JK) = AP (IMJ) 
               AP (IMJ) = T 
               IMJ = IMJ - (J - 1) 
   50       END DO 
            T = B (K) 
            B (K) = B (IMAX) 
            B (IMAX) = T 
         ENDIF 
!                                                                       
!           Perform elimination                                         
!                                                                       
         IJ = IK - (K - 1) 
         DO 60 JJ = 1, KM1 
            J = K - JJ 
            JK = IK + J 
            D1 = - AP (JK) / AP (KK) 
            CALL VECADD (J, D1, AP (IK + 1), AP (IJ + 1) ) 
            CALL VECADD (1, D1, B (K), B (J) ) 
            IJJ = IJ + J 
            IJ = IJ - (J - 1) 
   60    END DO 
!                                                                       
!           Adjust PIVOT vector                                         
!                                                                       
         KPVT (K) = K 
         IF (SWAP) KPVT (K) = IMAX 
      ENDIF 
      IK = IK - (K - 1) 
      IF (KSTEP.EQ.2) IK = IK - (K - 2) 
      K = K - KSTEP 
      GOTO 10 
      END SUBROUTINE ZSPMOK                         
