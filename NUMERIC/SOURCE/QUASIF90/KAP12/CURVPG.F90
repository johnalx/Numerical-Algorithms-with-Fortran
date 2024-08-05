      SUBROUTINE CURVPG (NP, N, M, K, KV, DP, X, XP, D, E, IERR) 
                                                                        
!********************************************************************   
!                                                                   *   
!   The subroutine CURVPG computes NP points of a closed uniform    *   
!   b spline curve of order K.                                      *   
!                                                                   *   
!   INPUT PARAMETERS                                                *   
!   ================                                                *   
!   NP      : desired number of points on the curve                 *   
!   N       : number of DE BOOR pointsis N - 1, N >= 2              *   
!   M       : Dimension of DE BOOR points (M >= 2)                  *   
!   K       : Order of the b spline, 3 <= K <= N+1                  *   
!   KV      : INTEGER vector KV(1:N+2*(K-1)) with the nodes         *   
!   DP      : DOUBLE PRECISION array DP(0:N+K-1,1:M); only N+1      *   
!             points used on input                                  *   
!                                                                   *   
!   AUX ARRAYS                                                      *   
!   ===========                                                     *   
!   X       : DOUBLE PRECISION vector X(1:M)                        *   
!   D,E     : DOUBLE PRECISION arrays ..(1:K,1:M)                   *   
!                                                                   *   
!   OUTPUT PARAMETERS                                               *   
!   =================                                               *   
!   XP      : DOUBLE PRECISION array XP(1:NP,1:M) with the computed *   
!             points; each row of XP contains the M coordinates of  *   
!             one point                                             *   
!   IERR    : error parameter                                       *   
!             IERR=0, all o k                                       *   
!             IERR=1, Input incorrect                               *   
!                                                                   *   
!-------------------------------------------------------------------*   
!                                                                   *   
!   Required subroutines: KNOTVG, DEBOOR                            *   
!                                                                   *   
!********************************************************************   
!                                                                   *   
!   Author      : Reinhold Wodicka, Bj”rn Terwege                    *  
!   Date        : 06.06.1995                                        *   
!   Source code : FORTRAN 77                                        *   
!                                                                   *   
!********************************************************************   
!                                                                       
!                                                                       
!  Declarations                                                         
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION KV (1:N + 2 * (K - 1) ), DP (0:N + K - 1, 1:M), X (1:M),&
      XP (1:NP, 1:M), D (1:K, 1:M), E (1:K, 1:M)                        
      IERR = 0 
!                                                                       
!  Check input                                                          
!                                                                       
      IF (N.LT.2.OR.K.LT.3.OR.K.GT.N + 1) THEN 
         IERR = 1 
         RETURN 
      ENDIF 
!                                                                       
!  Append the extra DE BOOR points to DP                                
!                                                                       
      J = 0 
      DO 10 I = N + 1, N + K - 1 
         DO 20 J = 1, M 
            DP (I, J) = DP (I - N - 1, J) 
   20    END DO 
   10 END DO 
!                                                                       
!  Call SUBROUTINE KNOTVG to compute the nodes                          
!                                                                       
      CALL KNOTVG (N, K, KV) 
!                                                                       
!  Determine step size                                                  
!                                                                       
      DT = DBLE (N + 1) / DBLE (NP - 1) 
      T = DBLE (K - 1) 
      IR = K - 1 
!                                                                       
!  compute points of the curve                                          
!                                                                       
      DO 30 I = 1, NP 
         CALL DEBOOR (N + K - 1, M, DP, K, KV, T, IR, D, E, X, IERR) 
!                                                                       
!  Save the coordinates of the newly computed point X in row I of XP    
!                                                                       
         DO 40 L = 1, M 
            XP (I, L) = X (L) 
   40    END DO 
         T = DMIN1 (T + DT, DBLE (N + K) ) 
   25    IF (T.GT.DBLE (IR + 1) ) THEN 
            IR = IR + 1 
            GOTO 25 
         ENDIF 
   30 END DO 
      RETURN 
      END SUBROUTINE CURVPG                         
