      SUBROUTINE RK547S (M, COEFF, QG) 
!                                                                       
!*****************************************************************      
!                                                                *      
! This subroutine determines the coefficients for the embedding  *      
! formula RK5(4)7S (method of Dormand und Prince).               *      
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
      COEFF (2, 1) = 2.0D0 / 9.0D0 
      COEFF (3, 1) = 1.0D0 / 3.0D0 
      COEFF (4, 1) = 5.0D0 / 9.0D0 
      COEFF (5, 1) = 2.0D0 / 3.0D0 
      COEFF (6, 1) = 1.0D0 
      COEFF (7, 1) = 1.0D0 
!                                                                       
!         b - values (see scheme)                                       
!                                                                       
      COEFF (2, 2) = 2.0D0 / 9.0D0 
      COEFF (3, 2) = 1.0D0 / 12.0D0 
      COEFF (4, 2) = 55.0D0 / 324.0D0 
      COEFF (5, 2) = 83.0D0 / 330.0D0 
      COEFF (6, 2) = - 19.0D0 / 28.0D0 
      COEFF (7, 2) = 19.0D0 / 200.0D0 
      COEFF (3, 3) = 1.0D0 / 4.0D0 
      COEFF (4, 3) = - 25.0D0 / 108.0D0 
      COEFF (5, 3) = - 13.0D0 / 22.0D0 
      COEFF (6, 3) = 9.0D0 / 4.0D0 
      COEFF (4, 4) = 50.0D0 / 81.0D0 
      COEFF (5, 4) = 61.0D0 / 66.0D0 
      COEFF (6, 4) = 1.0D0 / 7.0D0 
      COEFF (7, 4) = 3.0D0 / 5.0D0 
      COEFF (5, 5) = 9.0D0 / 110.0D0 
      COEFF (6, 5) = - 27.0D0 / 7.0D0 
      COEFF (7, 5) = - 243.0D0 / 400.0D0 
      COEFF (6, 6) = 22.0D0 / 7.0D0 
      COEFF (7, 6) = 33.0D0 / 40.0D0 
      COEFF (7, 7) = 7.0D0 / 80.0D0 
!                                                                       
!         A tilde values (see scheme)                                   
!                                                                       
      COEFF (2, 3) = 19.0D0 / 200.0D0 
      COEFF (2, 5) = 3.0D0 / 5.0D0 
      COEFF (2, 6) = - 243.0D0 / 400.0D0 
      COEFF (2, 7) = 33.0D0 / 40.0D0 
      COEFF (3, 4) = 7.0D0 / 80.0D0 
!                                                                       
!          A - values (see scheme)                                      
!                                                                       
                                                                        
      COEFF (1, 1) = 431.0D0 / 5000.0D0 
      COEFF (1, 3) = 333.0D0 / 500.0D0 
      COEFF (1, 4) = - 7857.0D0 / 10000.0D0 
      COEFF (1, 5) = 957.0D0 / 1000.0D0 
      COEFF (1, 6) = 193.0D0 / 2000.0D0 
      COEFF (1, 7) = - 1.0D0 / 50.0D0 
      RETURN 
      END SUBROUTINE RK547S                         
