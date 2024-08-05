      DOUBLEPRECISION FUNCTION SENORM (X, N) 
!                                                                       
!*****************************************************************      
!                                                                *      
!  The FUNCTION-subroutine SENORM computes the square of the     *      
!  euclidean norm of a real vector X of length N+1               *      
!                                                                *      
!  INPUT PARAMETERS:                                             *      
!  =================                                             *      
!                                                                *      
!  X    (N+1)-vector X(0:N) containing the vector whose norm     *      
!       shall be computed                                        *      
!  N    number of components of the vector X                     *      
!                                                                *      
!  OUTPUT PARAMETER:                                             *      
!  =================                                             *      
!                                                                *      
!  none                                                          *      
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
      INTEGER N 
      DOUBLEPRECISION X (0:N) 
      SENORM = 0.0D0 
      DO 10 I = 0, N 
         SENORM = SENORM + X (I) * X (I) 
   10 END DO 
      RETURN 
      END FUNCTION SENORM                           
