      DOUBLEPRECISION FUNCTION VMNORM (VECTOR, N) 
!                                                                       
!*****************************************************************      
!                                                                *      
! This FUNCTION routine determines the maximum norm of a vector. *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! VECTOR  : DOUBLE PRECISION vector VECTOR(1:N)                  *      
! N       : size of the vector                                   *      
!                                                                *      
!                                                                *      
! LOCAL VARIABLE:                                                *      
! ===============                                                *      
! I       : loop variable                                        *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  subroutines required: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author  : Volker KrÅger                                     *        
!  Date    : 07.08.1990                                        *        
!  Source  : FORTRAN 77                                        *        
!                                                                *      
!*****************************************************************      
!                                                                       
! Declaration                                                           
!                                                                       
      DOUBLEPRECISION VECTOR (N) 
!                                                                       
! Compute maximum norm                                                  
!                                                                       
      VMNORM = DABS (VECTOR (1) ) 
      DO 10 I = 2, N 
         VMNORM = DMAX1 (VMNORM, DABS (VECTOR (I) ) ) 
   10 END DO 
      RETURN 
      END FUNCTION VMNORM                           
