      SUBROUTINE KNOTVO (N, K, KV) 
!                                                                       
!*****************************************************************      
!                                                                *      
!                                                                *      
!   This program computes the node vector KV of an open uniform  *      
!   B spline.                                                    *      
!                                                                *      
!   INPUT PARAMETERS                                             *      
!   ================                                             *      
!   N       : N+1 is the number of DE BOOR points, N >= 2        *      
!   K       : degree of the B spline curve, 3 <= K <= N+1        *      
!                                                                *      
!   OUTPUT PARAMETERS                                            *      
!   =================                                            *      
!   KV      : INTEGER vector KV(1:N+K+1), with the node vector   *      
!             of order K                                         *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!   Required subroutines: none                                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author      : Gisela Engeln-Muellge                           *     
!  Date        : 11.30.91                                        *      
!  Source code : FORTRAN 77                                      *      
!                                                                *      
!*****************************************************************      
!                                                                       
      INTEGER KV (1:N + 1 + K) 
!                                                                       
!      Create the node vector KV of an open uniform B spline            
!                                                                       
      DO 10 J = 1, K - 1 
         KV (J) = K - 1 
   10 END DO 
      DO 20 J = K, N 
         KV (J) = J 
   20 END DO 
      DO 30 J = N + 1, N + K - 1 
         KV (J) = N + 1 
   30 END DO 
      RETURN 
      END SUBROUTINE KNOTVO                         
