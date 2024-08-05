      SUBROUTINE DEBOOR (N, M, DP, K, KV, T, IR, D, E, X, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This program computes the point corresponding to the parameter*      
!  T of a uniform b spline whose nodes are known by using the    *      
!  algorithm of DE BOOR.                                         *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  N      : N+1 is the number of DE BOOR points, N >= 2          *      
!  M      : Dimension of the space for the DE BOOR points; M >= 2*      
!  DP     : DOUBLE PRECISION array DP(0:N,1:M) with the coordi-  *      
!           nates of the DE BOOR points                          *      
!  K      : Order of the B spline, 3 <= K <= N+1                 *      
!  KV     : INTEGER vectorKV(1:N+1+K), with the node vector of   *      
!           order K                                              *      
!  T      : Parameter value for which we want to determine a     *      
!           point on the B spline                                *      
!  IR     : Index of the element of the node vector with         *      
!              KV(IR) <= T <= KV(IR+1); IR >= K-1                *      
!                                                                *      
!                                                                *      
!  AUX VECTORS:                                                  *      
!  ============                                                  *      
!  D,E    : DOUBLE PRECISION arrays ..(1:K,1:M)                  *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  X      : DOUBLE PRECISION vector X(1:M), with the coordinates *      
!           of the computed point on the B spline                *      
!  IERR   : Error parameter                                      *      
!           IERR=0, all ok                                       *      
!           IERR=1, error: N < 2 or K < 3 or K > N+1             *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  Required subroutines: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author      : Gisela Engeln-Muellges                          *      
!  Date        : 11.30.91                                        *      
!  Source code : FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      DOUBLEPRECISION DP (0:N, 1:M), D (1:K, 1:M), E (1:K, 1:M),        &
      X (1:M), T, ZA                                                    
      INTEGER KV (1:N + K - 1) 
!                                                                       
      IERR = 0 
!                                                                       
!     Check input:g: N >= 2, 3 <= K <= N+1                              
!                                                                       
      IF (N.LT.2.OR.K.LT.3.OR.K.GT.N + 1) THEN 
         IERR = 1 
         RETURN 
      ENDIF 
!                                                                       
!     Choose a set of K consecutive points from DP                      
!     Use the DE BOOR algorithm; store in D                             
!                                                                       
      DO 10 I = IR - K + 1, IR 
         DO 20 J = 1, M 
            D (I - IR + K, J) = DP (I, J) 
   20    END DO 
   10 END DO 
!                                                                       
!     Compute entries of E                                              
!                                                                       
      DO 30 L = 1, K - 1 
         DO 40 J = L + 1, K 
            ZA = T - DBLE (KV (J + IR - K) ) 
            NE = KV (J + IR - L) - KV (J + IR - K) 
            IF (NE.EQ.0) THEN 
               ALPHA = 0.0D0 
            ELSE 
               ALPHA = ZA / DBLE (NE) 
            ENDIF 
            DO 50 I = 1, M 
               E (J, I) = D (J - 1, I) + ALPHA * (D (J, I) - D (J - 1,  &
               I) )                                                     
   50       END DO 
   40    END DO 
         DO 60 I = 1, K 
            DO 70 J = 1, M 
               D (I, J) = E (I, J) 
   70       END DO 
   60    END DO 
   30 END DO 
!                                                                       
!     Store coordinates of the computed point (in aux vector D) in X    
!                                                                       
      DO 80 I = 1, M 
         X (I) = D (K, I) 
   80 END DO 
      RETURN 
      END SUBROUTINE DEBOOR                         
