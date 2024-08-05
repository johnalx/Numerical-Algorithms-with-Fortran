      SUBROUTINE RK658M (M, COEFF, QG) 
!                                                                       
!*****************************************************************      
!                                                                *      
! This subroutine determines the coefficients for the embedding  *      
! formula RK6(5)8M (method of Dormand und Prince).               *      
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
!             High order embedded Runge-Kutta formulae           *      
!             page 67-75                                         *      
!             Journal of Computational and Applied Mathematics,  *      
!             volume 7, no 1, 1981                               *      
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
      COEFF (2, 1) = 1.0D0 / 10.0D0 
      COEFF (3, 1) = 2.0D0 / 9.0D0 
      COEFF (4, 1) = 3.0D0 / 7.0D0 
      COEFF (5, 1) = 3.0D0 / 5.0D0 
      COEFF (6, 1) = 4.0D0 / 5.0D0 
      COEFF (7, 1) = 1.0D0 
      COEFF (8, 1) = 1.0D0 
!                                                                       
!         b - values (see scheme)                                       
!                                                                       
      COEFF (2, 2) = 1.0D0 / 10.0D0 
      COEFF (3, 2) = - 2.0D0 / 81.0D0 
      COEFF (4, 2) = 615.0D0 / 1372.0D0 
      COEFF (5, 2) = 3243.0D0 / 5500.0D0 
      COEFF (6, 2) = - 26492.0D0 / 37125.0D0 
      COEFF (7, 2) = 5561.0D0 / 2376.0D0 
      COEFF (8, 2) = 465467.0D0 / 266112.0D0 
      COEFF (3, 3) = 20.0D0 / 81.0D0 
      COEFF (4, 3) = - 270.0D0 / 343.0D0 
      COEFF (5, 3) = - 54.0D0 / 55.0D0 
      COEFF (6, 3) = 72.0D0 / 55.0D0 
      COEFF (7, 3) = - 35.0D0 / 11.0D0 
      COEFF (8, 3) = - 2945.0D0 / 1232.0D0 
      COEFF (4, 4) = 1053.0D0 / 1372.0D0 
      COEFF (5, 4) = 50949.0D0 / 71500.0D0 
      COEFF (6, 4) = 2808.0D0 / 23375.0D0 
      COEFF (7, 4) = - 24117.0D0 / 31603.0D0 
      COEFF (8, 4) = - 561.0201D4 / 1415.8144D4 
      COEFF (5, 5) = 4998.0D0 / 17875.0D0 
      COEFF (6, 5) = - 24206.0D0 / 37125.0D0 
      COEFF (7, 5) = 899983.0D0 / 200772.0D0 
      COEFF (8, 5) = 1051.3573D4 / 321.2352D4 
      COEFF (6, 6) = 338.0D0 / 459.0D0 
      COEFF (7, 6) = - 5225.0D0 / 1836.0D0 
      COEFF (8, 6) = - 424325.0D0 / 205632.0D0 
      COEFF (7, 7) = 3925.0D0 / 4056.0D0 
      COEFF (8, 7) = 376225.0D0 / 454272.0D0 
!                                                                       
!         A tilde values (see scheme)                                   
!                                                                       
      COEFF (2, 3) = 61.0D0 / 864.0D0 
      COEFF (2, 5) = 98415.0D0 / 321776.0D0 
      COEFF (2, 6) = 16807.0D0 / 146016.0D0 
      COEFF (2, 7) = 1375.0D0 / 7344.0D0 
      COEFF (2, 8) = 1375.0D0 / 5408.0D0 
      COEFF (3, 4) = - 37.0D0 / 1120.0D0 
      COEFF (3, 5) = 1.0D0 / 10.0D0 
!                                                                       
!          A - values (see scheme)                                      
!                                                                       
      COEFF (1, 1) = 821.0D0 / 10800.0D0 
      COEFF (1, 3) = 19683.0D0 / 71825.0D0 
      COEFF (1, 4) = 175273.0D0 / 912600.0D0 
      COEFF (1, 5) = 395.0D0 / 3672.0D0 
      COEFF (1, 6) = 785.0D0 / 2704.0D0 
      COEFF (1, 7) = 3.0D0 / 50.0D0 
      RETURN 
      END SUBROUTINE RK658M                         
