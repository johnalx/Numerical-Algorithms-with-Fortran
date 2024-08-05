      SUBROUTINE RK658S (M, COEFF, QG) 
!                                                                       
!*****************************************************************      
!                                                                *      
! This subroutine determines the coefficients for the embedding  *      
! formula RK6(5)8S (method of Dormand und Prince).               *      
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
!             A reconsideration of some embedded Runge Kutta     *      
!             formulae                                           *      
!             page 203-211                                       *      
!             Journal of Computational and Applied Mathematics,  *      
!             15, 1986, North-Holland                            *      
!                                                                *      
!*****************************************************************      
!                                                                       
! Declarations                                                          
!                                                                       
      DOUBLEPRECISION COEFF (16, M), QG 
!                                                                       
! Initialize QG                                                         
!                                                                       
      QG = 5.0D0 
!                                                                       
! Determine the matrix elements in COEFF                                
!                                                                       
!         a - values (see scheme )                                      
!                                                                       
      COEFF (2, 1) = 1.0D0 / 4.0D0 
      COEFF (3, 1) = 3.0D0 / 10.0D0 
      COEFF (4, 1) = 6.0D0 / 7.0D0 
      COEFF (5, 1) = 3.0D0 / 5.0D0 
      COEFF (6, 1) = 4.0D0 / 5.0D0 
      COEFF (7, 1) = 1.0D0 
      COEFF (8, 1) = 1.0D0 
!                                                                       
!         b - values (see scheme)                                       
!                                                                       
      COEFF (2, 2) = 1.0D0 / 4.0D0 
      COEFF (3, 2) = 3.0D0 / 25.0D0 
      COEFF (4, 2) = 102.0D0 / 343.0D0 
      COEFF (5, 2) = - 3.0D0 / 100.0D0 
      COEFF (6, 2) = 37.0D0 / 225.0D0 
      COEFF (7, 2) = 11.0D0 / 648.0D0 
      COEFF (8, 2) = 796.0D0 / 1701.0D0 
      COEFF (3, 3) = 9.0D0 / 50.0D0 
      COEFF (4, 3) = - 1368.0D0 / 343.0D0 
      COEFF (5, 3) = 36.0D0 / 25.0D0 
      COEFF (6, 3) = - 48.0D0 / 25.0D0 
      COEFF (7, 3) = 14.0D0 / 3.0D0 
      COEFF (8, 3) = - 352.0D0 / 63.0D0 
      COEFF (4, 4) = 1560.0D0 / 343.0D0 
      COEFF (5, 4) = - 12.0D0 / 13.0D0 
      COEFF (6, 4) = 872.0D0 / 351.0D0 
      COEFF (7, 4) = - 10193.0D0 / 2106.0D0 
      COEFF (8, 4) = 134096.0D0 / 22113.0D0 
      COEFF (5, 5) = 147.0D0 / 1300.0D0 
      COEFF (6, 5) = 49.0D0 / 1053.0D0 
      COEFF (7, 5) = - 30331.0D0 / 50544.0D0 
      COEFF (8, 5) = - 78281.0D0 / 75816.0D0 
      COEFF (6, 6) = 2.0D0 / 81.0D0 
      COEFF (7, 6) = 1025.0D0 / 1944.0D0 
      COEFF (8, 6) = - 9425.0D0 / 20412.0D0 
      COEFF (7, 7) = 59.0D0 / 48.0D0 
      COEFF (8, 7) = 781.0D0 / 504.0D0 
!                                                                       
!         A tilde values (see scheme)                                   
!                                                                       
      COEFF (2, 3) = 29.0D0 / 324.0D0 
      COEFF (2, 5) = 3400.0D0 / 7371.0D0 
      COEFF (2, 6) = - 16807.0D0 / 25272.0D0 
      COEFF (2, 7) = - 125.0D0 / 1944.0D0 
      COEFF (2, 8) = 25.0D0 / 24.0D0 
      COEFF (3, 4) = 1.0D0 / 84.0D0 
      COEFF (3, 5) = 1.0D0 / 8.0D0 
!                                                                       
!          A - values (see scheme)                                      
!                                                                       
      COEFF (1, 1) = 2041.0D0 / 21600.0D0 
      COEFF (1, 3) = 748.0D0 / 1755.0D0 
      COEFF (1, 4) = - 2401.0D0 / 46800.0D0 
      COEFF (1, 5) = 11.0D0 / 108.0D0 
      COEFF (1, 6) = 59.0D0 / 160.0D0 
      COEFF (1, 7) = 3.0D0 / 50.0D0 
      RETURN 
      END SUBROUTINE RK658S                         
