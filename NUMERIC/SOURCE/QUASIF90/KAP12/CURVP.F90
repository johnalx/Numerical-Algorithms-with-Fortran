      SUBROUTINE CURVP (NP, N, M, K, KV, DP, X, XP, D, E, IERR, IX) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  The subroutine CURVP computes at most NP+1 points on an open  *      
!  uniform B-Spline curve of order K.                            *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!  NP     : number of desired points on the curve:               *      
!           maximally NP+1 points and at least 2*(N-K+2)+1 points*      
!           are computed, i.e., NP=MAX(NP,2*(N-K+2))             *      
!  N      : N+1 is the number of DE BOOR points                  *      
!  M      : dimension of the space for the DE BOOR points        *      
!           ( M >= 2)                                            *      
!  K      : order of the B-Spline curve (3 <= K <= N+1)          *      
!  KV     : INTEGER vector KV(1:N+K-1) for the node vectors      *      
!  DP     : DOUBLE PRECISION array DP(0:N,1:M) containing the    *      
!           N+1 DE BOOR points, N >= 2                           *      
!                                                                *      
!                                                                *      
!  AUXILIARY VARIABLES:                                          *      
!  ====================                                          *      
!  X      : DOUBLE PRECISION vector X(1:M)                       *      
!  D,E    : DOUBLE PRECISION 2 dimensional arrays D(1:K,1:M),    *      
!           E(1:K,1:M)                                           *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  XP     : DOUBLE PRECISION array XP(1:NP+1,1:M) containing the *      
!           computed points on the curve (each row contains the M*      
!           coordinates of a point; their maximal number is NP+1)*      
!  IERR   : error parameter                                      *      
!           IERR=0, everything is o.k.                           *      
!           IERR=1, input conditions violated                    *      
!  IX     : Number of points stored in XP                        *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  required subroutines: DEBOOR, KNOTVO                          *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : Gisela Engeln-MÅllges                              *      
!  Date     : 11.30.1991                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
!     Declarations                                                      
!                                                                       
      DOUBLEPRECISION DP (0:N, 12:M), X (1:M), XP (1:NP + 1, 1:M),      &
      D (1:K, 1:M), E (1:K, 1:M), H, T                                  
      INTEGER KV (1:N + K - 1), IX 
!                                                                       
!     Checking the input data                                           
!                                                                       
      IF (N.LT.2.OR.K.LT.3.OR.K.GT.N + 1) THEN 
         IERR = 1 
         RETURN 
      ENDIF 
!                                                                       
!     Call SUBROUTINE KNOTVO in order to compute the nodes              
!                                                                       
      CALL KNOTVO (N, K, KV) 
!                                                                       
!     Compute the number of points between two adjacent nodes           
!     of the node vector and the step size H                            
!                                                                       
      NI = INT (NP / (N - K + 2) ) 
      NI = MAX (NI, 2) 
      H = 1.0D0 / DBLE (NI) 
!                                                                       
!     Call SUBROUTINE DEBOOR in order to compute the coordinates        
!     of the first point of the curve; store in the first row of XP     
!                                                                       
      CALL DEBOOR (N, M, DP, K, KV, DBLE (K - 1), K - 1, D, E, X, IERR) 
!                                                                       
      DO 10 I = 1, M 
         XP (1, I) = X (I) 
   10 END DO 
!                                                                       
!     Compute the subsequent points of the curve                        
!                                                                       
      IX = 1 
      DO 20 IR = K - 1, N 
         T = DBLE (IR) 
         DO 30 L = 1, NI 
            T = T + H 
            CALL DEBOOR (N, M, DP, K, KV, T, IR, D, E, X, IERR) 
            IX = IX + 1 
!                                                                       
!           Store the coordinates of the computed point X on the curve  
!           in row IX of XP, 2 <= IX <= NP+1                            
!                                                                       
            DO 40 I = 1, M 
               XP (IX, I) = X (I) 
   40       END DO 
   30    END DO 
   20 END DO 
      RETURN 
      END SUBROUTINE CURVP                          
