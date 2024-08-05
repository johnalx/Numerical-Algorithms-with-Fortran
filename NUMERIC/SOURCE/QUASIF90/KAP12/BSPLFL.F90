      SUBROUTINE BSPLFL (KK, NU, NV, M, N, KU, KV, KVU, KVV, DP, XP,    &
      UPOL, D1, D0, E1, E0, HILF, X, IERR)                              
!********************************************************************   
!                                                                   *   
!   This subroutine generates a mesh for a b spline surface from    *   
!   given u and v curves.                                           *   
!                                                                   *   
!   INPUT PARAMETERS                                                *   
!   ================                                                *   
!   KK      : Order                                                 *   
!   NU      : Number of nodes on a u curve, NU >= 2; NU - 1 denotes *   
!             the number of u intervals                             *   
!   NV      : Number of nodes on a v curve, NV >= 2; NV - 1 denotes *   
!             the number of v intervals                             *   
!   M       : M+1 is the number of v polygons; M >= 2               *   
!   N       : N+1 is the number of u polygons; N >= 2               *   
!   KU      : Order of the u curves, 2 <= KU <= M+1                 *   
!   KV      : Order of the V curves, 2 <= KV <= N+1                 *   
!   KVU     : INTEGER vector KVU(1:KU+M-1); the node vector KVU for *   
!             the open u curves of order KU                         *   
!   KVV     : INTEGER vector KVV(1:KV+N-1); the node vector KVV for *   
!             the open v curves of order KV                         *   
!   DP      : 3 dimensional DOUBLE PRECISION array DP(0:M,0:N,1:KK) *   
!             containing the DE BOOR polytope with the u and v poly-*   
!             gons                                                  *   
!                                                                   *   
!                                                                   *   
!   AUX ARRAYS                                                      *   
!   ===========                                                     *   
!   UPOL   : 2 dim. DOUBLE PRECISION array UPOL(1:M,1:KK)           *   
!   D0,E0  : 2 dim. DOUBLE PRECISION arays ..(1:KU,1:KK)            *   
!   D1,E1  : 2 dim. DOUBLE PRECISION arrays ..(1:KV,1:KK)           *   
!   DUMMY  : 2 dim. DOUBLE PRECISION array DUMMY(0:N,1:KK)          *   
!   X      : DOUBLE PRECISION vector X(1:KK)                        *   
!                                                                   *   
!                                                                   *   
!   OUTPUT PARAMETERS                                               *   
!   =================                                               *   
!   XP      : 3 dim. DOUBLE PRECISION array XP(1:N,1:NN,1:KK) with  *   
!             the generated mesh points                             *   
!   IERR    : Error parameter                                       *   
!             IERR=0 : all is ok                                    *   
!             IERR=1 : invalid input                                *   
!                                                                   *   
!-------------------------------------------------------------------*   
!                                                                   *   
!   Required subroutines:   DEBOOR, KNOTVO                          *   
!                                                                   *   
!********************************************************************   
!                                                                   *   
!   Authors    : Reinhold Wodicka, Bjoern Terwege                   *   
!   Date       : 6.12.1995                                          *   
!   Sourcecode : FORTRAN 77                                         *   
!                                                                   *   
!********************************************************************   
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION DP (0:M, 0:N, 1:KK), KVU (1:KU + M - 1), KVV (1:KV + N -&
      1), XP (1:NU, 1:NV, 1:M), X (1:KK), UPOL (0:M, 1:KK), D1 (1:KV, 1:&
      KK), E1 (1:KV, 1:KK), D0 (1:KU, 1:KK), E0 (1:KU, 1:KK), HILF (0:N,&
      1:KK)                                                             
!                                                                       
!  Stopping criteria:                                                   
!                                                                       
      IERR = 0 
      IF (NU.LT.2.OR.NV.LT.2.OR.M.LT.2.OR.N.LT.2) THEN 
         IERR = 1 
         RETURN 
      ENDIF 
!                                                                       
!  Compute node vectors                                                 
!                                                                       
      CALL KNOTVO (N, KV, KVV) 
      CALL KNOTVO (M, KU, KVU) 
!                                                                       
!  Compute stepsizes                                                    
!                                                                       
      DU = DBLE (M + 2 - KU) / DBLE (NU - 1) 
      DV = DBLE (N + 2 - KV) / DBLE (NV - 1) 
!                                                                       
!  Starting parameter                                                   
!                                                                       
      VJ = DBLE (KV - 1) 
!                                                                       
!  Starting index                                                       
!                                                                       
      IS = KV - 1 
      DO 10 J = 1, NV 
!                                                                       
!  Generate one u curve for each vj, j=1, ...., Nv                      
!                                                                       
         DO 20 L = 0, M 
!                                                                       
!  Start by constructing a u polygon for vj                             
!                                                                       
!  Prepare v polygon by storing separately                              
!                                                                       
            DO 1 L1 = 0, N 
               DO 2 L2 = 1, KK 
                  HILF (L1, L2) = DP (L, L1, L2) 
    2          END DO 
    1       END DO 
            CALL DEBOOR (N, KK, HILF, KV, KVV, VJ, IS, D1, E1, X, IERR) 
            DO 22 II = 1, 3 
               UPOL (L, II) = X (II) 
   22       END DO 
   20    END DO 
!                                                                       
!  v polygon uPol has been computed; it shall be used for the u curve   
!                                                                       
!                                                                       
!  Starting parameter                                                   
!                                                                       
         UI = DBLE (KU - 1) 
!                                                                       
!  Starting index                                                       
!                                                                       
         IR = KU - 1 
         DO 30 I = 1, NU 
            CALL DEBOOR (M, KK, UPOL, KU, KVU, UI, IR, D0, E0, X, IERR) 
            DO 33 II = 1, KK 
               XP (I, J, II) = X (II) 
   33       END DO 
!                                                                       
!  nex parameter, m+1 is the maximal u parameter                        
!                                                                       
            UI = DMIN1 (UI + DU, DBLE (M + 1) ) 
!                                                                       
!  Find next index                                                      
!                                                                       
   35       IF (UI.GT.IR + 1) THEN 
               IR = IR + 1 
               GOTO 35 
            ENDIF 
   30    END DO 
!                                                                       
!  Nu points of the u curve for v = vj have been found                  
!                                                                       
!                                                                       
!  next v parameter;  n+1 is  maximal v parameter                       
!                                                                       
         VJ = DMIN1 (VJ + DV, DBLE (N + 1) ) 
!                                                                       
!  Find next index                                                      
!                                                                       
   40    IF (VJ.GT.IS + 1) THEN 
            IS = IS + 1 
            GOTO 40 
         ENDIF 
!                                                                       
!  Nv u curves have been computed                                       
!                                                                       
   10 END DO 
      RETURN 
      END SUBROUTINE BSPLFL                         
