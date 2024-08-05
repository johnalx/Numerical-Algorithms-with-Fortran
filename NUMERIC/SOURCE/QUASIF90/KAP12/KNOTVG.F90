      SUBROUTINE KNOTVG (N, K, KV) 
!********************************************************************   
!                                                                   *   
!   This program computes the node vector KV of a closed uniform B  *   
!   spline.                                                         *   
!                                                                   *   
!   INPUT PARAMETERS                                                *   
!   ================                                                *   
!   N       : N+1 is the number of DE BOOR points, N >= 2           *   
!   K       : degree of the B spline curve, 2 <= K <= N+1           *   
!                                                                   *   
!   OUTPUT PARAMETERS                                               *   
!   =================                                               *   
!   KV      : INTEGER vector KV(1:N+K-1), with the node vector of   *   
!             order K                                               *   
!                                                                   *   
!-------------------------------------------------------------------*   
!                                                                   *   
!   Required subroutines: none                                      *   
!                                                                   *   
!********************************************************************   
!                                                                   *   
!   Authors     : Reinhold Wodicka, Bjoern Terwege                  *   
!   Date        : 05.07.1995                                        *   
!   Source code : FORTRAN 77                                        *   
!                                                                   *   
!********************************************************************   
!                                                                       
      INTEGER KV (1:N + 2 * (K - 1) ) 
!                                                                       
!     Create the node vector KV of a closed uniform B spline            
!                                                                       
      DO 10 J = 1, N + 2 * (K - 1) 
         KV (J) = j 
   10 END DO 
      RETURN 
      END SUBROUTINE KNOTVG                         
