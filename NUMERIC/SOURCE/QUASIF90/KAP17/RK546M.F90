      SUBROUTINE RK546M (M, COEFF, QG) 
!                                                                       
!*****************************************************************      
!                                                                *      
! This subroutine determines the coefficients for the embedding  *      
! formula RK5(4)6M (method of Dormand und Prince).               *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! M       : Dimension of the matrix COEFF depending on the chosen*      
!           embedding formula                                    *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! COEFF   : 2-dimensional DOUBLE PRECISION array COEFF(1:16,1:M) *      
!           with the coefficients for the embedding formula      *      
! QG      : DOUBLE PRECISION value for the global error order    *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  required subroutines: none                                    *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : Volker KrÅger                                      *      
!  Date     : 26.04.1993                                         *      
!  Source   : FORTRAN 77                                         *      
!  Sources  : J. R. Dormand P. J. Prince :                       *      
!             A family of embedded Runge-Kutta formulae          *      
!             page 19-26                                         *      
!             Journal of Computational and Applied Mathematics,  *      
!             volume 6, no 1, 1980                               *      
!                                                                *      
!*****************************************************************      
!                                                                       
! Declarations                                                          
!                                                                       
      DOUBLEPRECISION COEFF (16, M), QG 
!                                                                       
! Initialize QG                                                         
!                                                                       
      QG = 4.0D0 
!                                                                       
! Determine the matrix elements in COEFF                                
!                                                                       
!         a - values (see scheme )                                      
!                                                                       
      COEFF (2, 1) = 1.0D0 / 5.0D0 
      COEFF (3, 1) = 3.0D0 / 10.0D0 
      COEFF (4, 1) = 3.0D0 / 5.0D0 
      COEFF (5, 1) = 2.0D0 / 3.0D0 
      COEFF (6, 1) = 1.0D0 
!                                                                       
!         b - values (see scheme)                                       
!                                                                       
      COEFF (2, 2) = 1.0D0 / 5.0D0 
      COEFF (3, 2) = 3.0D0 / 40.0D0 
      COEFF (4, 2) = 3.0D0 / 10.0D0 
      COEFF (5, 2) = 226.0D0 / 729.0D0 
      COEFF (6, 2) = - 181.0D0 / 270.0D0 
      COEFF (3, 3) = 9.0D0 / 40.0D0 
      COEFF (4, 3) = - 9.0D0 / 10.0D0 
      COEFF (5, 3) = - 25.0D0 / 27.0D0 
      COEFF (6, 3) = 5.0D0 / 2.0D0 
      COEFF (4, 4) = 6.0D0 / 5.0D0 
      COEFF (5, 4) = 880.0D0 / 729.0D0 
      COEFF (6, 4) = - 266.0D0 / 297.0D0 
      COEFF (5, 5) = 55.0D0 / 729.0D0 
      COEFF (6, 5) = - 91.0D0 / 27.0D0 
      COEFF (6, 6) = 189.0D0 / 55.0D0 
!                                                                       
!         A tilde values (see scheme)                                   
!                                                                       
      COEFF (2, 3) = 19.0D0 / 216.0D0 
      COEFF (2, 5) = 1000.0D0 / 2079.0D0 
      COEFF (2, 6) = - 125.0D0 / 216.0D0 
      COEFF (3, 4) = 81.0D0 / 88.0D0 
      COEFF (3, 5) = 5.0D0 / 56.0D0 
!                                                                       
!          A - values (see scheme)                                      
!                                                                       
      COEFF (1, 1) = 31.0D0 / 540.0D0 
      COEFF (1, 3) = 190.0D0 / 297.0D0 
      COEFF (1, 4) = - 145.0D0 / 108.0D0 
      COEFF (1, 5) = 351.0D0 / 220.0D0 
      COEFF (1, 6) = 1.0D0 / 20.0D0 
      RETURN 
      END SUBROUTINE RK546M                         
