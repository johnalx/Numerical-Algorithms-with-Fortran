      SUBROUTINE SLPRE (X, W, IWFL, FCT, LDA, M, N, F, A) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  The SUBROUTINE SLPRE forms the matrix A made up of function   *      
!  values as needed for a discrete linear least square problem.  *      
!                                                                *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  -----------------                                             *      
!                                                                *      
!  The parameters X, W, IWFL, FCT, LDA, M, and N are of the same *      
!  as those in SUBROUTINE SLFIT                                  *      
!                                                                *      
!                                                                *      
!  AUXILIARY PARAMETERS:                                         *      
!  ---------------------                                         *      
!                                                                *      
!  F  (N+1)-vector F(0:N) used for calling SUBROUTINE  FCT       *      
!                                                                *      
!                                                                *      
!  OUTPUT PARAMETER:                                             *      
!  -----------------                                             *      
!                                                                *      
!  A  2-dim. array A(0:LDA,0:N) containing the values of the     *      
!     model functions:                                           *      
!     A(I,K) = value of the K-th model function at node X(I)     *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Ilona Westermann                                   *      
!  date     : 09.01.1987                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
      DIMENSION A (0:LDA, 0:N), F (0:N), X (0:M), W (0:M) 
      IF (IWFL.EQ.0) THEN 
         DO 10 I = 0, M 
            WI = DSQRT (W (I) ) 
            CALL FCT (X (I), N, F) 
            DO 10 K = 0, N 
               A (I, K) = F (K) * WI 
   10    CONTINUE 
      ELSE 
         DO 20 I = 0, M 
            CALL FCT (X (I), N, F) 
            DO 20 K = 0, N 
               A (I, K) = F (K) 
   20    CONTINUE 
      ENDIF 
      RETURN 
      END SUBROUTINE SLPRE                          
