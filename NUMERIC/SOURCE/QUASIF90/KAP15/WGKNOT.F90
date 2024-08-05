      SUBROUTINE WGKNOT (N, WORK, IERR) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  This subroutine determines the weights and the nodes of the   *      
!  CLENSHAW-CURTIS quadrature formula of local error order N+3   *      
!  for the reference interval [-1,1].                            *      
!                                                                *      
!  INPUT PARAMETER:                                              *      
!  ================                                              *      
!  N    : N+1 denotes the number of nodes and weights,           *      
!         N >= 2, and N must be even.                            *      
!                                                                *      
!  OUTPUT PARAMETERS:                                            *      
!  ==================                                            *      
!  WORK : 2-dimensional array WORK(0:N,2) containing the weights *      
!         and nodes of the quadrature formula for the interval   *      
!         [-1,1]; the weights appear in the first column, with   *      
!         nodes in the second column of WORK.                    *      
!  IERR : error parameter.                                       *      
!           IERR = 0: everything o.k.                            *      
!           IERR = 1: condition for N is not met                 *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required : none                                   *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author    : Gisela Engeln-Muellges                            *      
!  date      : 05.10.1989                                        *      
!  source    : FORTRAN 77                                        *      
!                                                                *      
!*****************************************************************      
!                                                                       
!  declarations                                                         
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION WORK (0:N, 2) 
!                                                                       
!  initializing                                                         
!                                                                       
      PI = 4.0D0 * DATAN (1.0D0) 
!                                                                       
!  checking the validity of N                                           
!                                                                       
      IERR = 1 
      IF (N.LT.2.OR.MOD (N, 2) .NE.0) RETURN 
!                                                                       
!  determine the weights and nodes                                      
!                                                                       
      IERR = 0 
      DUMMY = N * N - 1 
      WORK (0, 1) = 1.0D0 / DUMMY 
      WORK (N, 1) = WORK (0, 1) 
      DUMMY1 = PI / N 
      WORK (0, 2) = 1.0D0 
      WORK (N, 2) = - 1.0D0 
      DO 10 K = 1, N - 1 
         DUMMY2 = 2.0D0 * (DUMMY - ( - 1) **K) / (N * DUMMY) 
         DUMMY3 = 0.0D0 
         DO 20 L = 1, N / 2 - 1 
            DUMMY3 = DUMMY3 + DCOS (2.0D0 * L * K * DUMMY1) / (4 * L *  &
            L - 1)                                                      
   20    END DO 
         WORK (K, 1) = DUMMY2 - 4.0D0 / N * DUMMY3 
         WORK (K, 2) = DCOS (K * DUMMY1) 
   10 END DO 
      RETURN 
      END SUBROUTINE WGKNOT                         
