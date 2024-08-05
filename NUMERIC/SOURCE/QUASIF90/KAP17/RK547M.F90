      SUBROUTINE RK547M (M, COEFF, QG) 
!                                                                       
!*****************************************************************      
!                                                                *      
! This subroutine determines the coefficients for the embedding  *      
! formula RK5(4)7M (method of Dormand und Prince).               *      
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
      COEFF (4, 1) = 4.0D0 / 5.0D0 
      COEFF (5, 1) = 8.0D0 / 9.0D0 
      COEFF (6, 1) = 1.0D0 
      COEFF (7, 1) = 1.0D0 
!                                                                       
!         b - values (see scheme)                                       
!                                                                       
      COEFF (2, 2) = 1.0D0 / 5.0D0 
      COEFF (3, 2) = 3.0D0 / 40.0D0 
      COEFF (4, 2) = 44.0D0 / 45.0D0 
      COEFF (5, 2) = 19372.0D0 / 6561.0D0 
      COEFF (6, 2) = 9017.0D0 / 3168.0D0 
      COEFF (7, 2) = 35.0D0 / 384.0D0 
      COEFF (3, 3) = 9.0D0 / 40.0D0 
      COEFF (4, 3) = - 56.0D0 / 15.0D0 
      COEFF (5, 3) = - 25360.0D0 / 2187.0D0 
      COEFF (6, 3) = - 355.0D0 / 33.0D0 
      COEFF (4, 4) = 32.0D0 / 9.0D0 
      COEFF (5, 4) = 64448.0D0 / 6561.0D0 
      COEFF (6, 4) = 46732.0D0 / 5247.0D0 
      COEFF (7, 4) = 500.0D0 / 1113.0D0 
      COEFF (5, 5) = - 212.0D0 / 729.0D0 
      COEFF (6, 5) = 49.0D0 / 176.0D0 
      COEFF (7, 5) = 125.0D0 / 192.0D0 
      COEFF (6, 6) = - 5103.0D0 / 18656.0D0 
      COEFF (7, 6) = - 2187.0D0 / 6784.0D0 
      COEFF (7, 7) = 11.0D0 / 84.0D0 
!                                                                       
!         A tilde values (see scheme)                                   
!                                                                       
      COEFF (2, 3) = 35.0D0 / 384.0D0 
      COEFF (2, 5) = 500.0D0 / 1113.0D0 
      COEFF (2, 6) = 125.0D0 / 192.0D0 
      COEFF (2, 7) = - 2187.0D0 / 6784.0D0 
      COEFF (3, 4) = 11.0D0 / 84.0D0 
!                                                                       
!          A - values (see scheme)                                      
!                                                                       
      COEFF (1, 1) = 5179.0D0 / 57600.0D0 
      COEFF (1, 3) = 7571.0D0 / 16695.0D0 
      COEFF (1, 4) = 393.0D0 / 640.0D0 
      COEFF (1, 5) = - 92097.0D0 / 339200.0D0 
      COEFF (1, 6) = 187.0D0 / 2100.0D0 
      COEFF (1, 7) = 1.0D0 / 40.0D0 
      RETURN 
      END SUBROUTINE RK547M                         
