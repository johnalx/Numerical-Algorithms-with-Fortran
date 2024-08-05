      DOUBLEPRECISION FUNCTION FENORM (N, X) 
!                                                                       
!*****************************************************************      
!                                                                *      
! Determining the euclidean norm of a vector X.                  *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! X       : N-vector X(1:N); vector for which the euclidean norm *      
!           is to be determined                                  *      
! N       : size of the vector X                                 *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETER:                                              *      
! =================                                              *      
! FENORM  : euclidean norm of the vector X                       *      
!                                                                *      
!                                                                *      
! LOCAL VARIABLES:                                               *      
! ================                                               *      
! SCPROD  : dot product of the vector X with itself              *      
! I       : iteration variable                                   *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  author   : Thomas Eul                                         *      
!  date     : 05.03.1985                                         *      
!  source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
      IMPLICIT DOUBLEPRECISION (A - H, O - Z) 
!                                                                       
!   Input declarations                                                  
!                                                                       
      INTEGER N 
      DOUBLEPRECISION X (N) 
!                                                                       
!   local variables                                                     
!                                                                       
      DOUBLEPRECISION SCPROD 
      INTEGER I 
!                                                                       
      SCPROD = 0.0D0 
      FENORM = 0.0D0 
      DO 10 I = 1, N 
         SCPROD = SCPROD+X (I) * X (I) 
   10 END DO 
      FENORM = DSQRT (SCPROD) 
!                                                                       
      RETURN 
      END FUNCTION FENORM                           
