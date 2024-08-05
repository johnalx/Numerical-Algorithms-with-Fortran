      SUBROUTINE RK547C (M, COEFF, QG) 
!                                                                       
!*****************************************************************      
!                                                                *      
! This subroutine determines the coefficients for the embedding  *      
! formula RK5(4)7C (method of Dormand und Prince).               *      
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
!              A reconsideration of some embedded Runge Kutta    *      
!              formulae                                          *      
!              page 203-211                                      *      
!              Journal of Computational and Applied Mathematics, *      
!              15, 1986, North-Holland                           *      
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
      COEFF (4, 1) = 6.0D0 / 13.0D0 
      COEFF (5, 1) = 2.0D0 / 3.0D0 
      COEFF (6, 1) = 1.0D0 
      COEFF (7, 1) = 1.0D0 
!                                                                       
!         b - values (see scheme)                                       
!                                                                       
      COEFF (2, 2) = 1.0D0 / 5.0D0 
      COEFF (3, 2) = 3.0D0 / 40.0D0 
      COEFF (4, 2) = 264.0D0 / 2197.0D0 
      COEFF (5, 2) = 932.0D0 / 3645.0D0 
      COEFF (6, 2) = - 367.0D0 / 513.0D0 
      COEFF (7, 2) = 35.0D0 / 432.0D0 
      COEFF (3, 3) = 9.0D0 / 40.0D0 
      COEFF (4, 3) = - 90.0D0 / 2197.0D0 
      COEFF (5, 3) = - 14.0D0 / 27.0D0 
      COEFF (6, 3) = 30.0D0 / 19.0D0 
      COEFF (4, 4) = 840.0D0 / 2197.0D0 
      COEFF (5, 4) = 3256.0D0 / 5103.0D0 
      COEFF (6, 4) = 9940.0D0 / 5643.0D0 
      COEFF (7, 4) = 8500.0D0 / 14553.0D0 
      COEFF (5, 5) = 7436.0D0 / 25515.0D0 
      COEFF (6, 5) = - 29575.0D0 / 8208.0D0 
      COEFF (7, 5) = - 28561.0D0 / 84672.0D0 
      COEFF (6, 6) = 6615.0D0 / 3344.0D0 
      COEFF (7, 6) = 405.0D0 / 704.0D0 
      COEFF (7, 7) = 19.0D0 / 196.0D0 
!                                                                       
!         A tilde values (see scheme)                                   
!                                                                       
      COEFF (2, 3) = 35.0D0 / 432.0D0 
      COEFF (2, 5) = 8500.0D0 / 14553.0D0 
      COEFF (2, 6) = - 28561.0D0 / 84672.0D0 
      COEFF (2, 7) = 405.0D0 / 704.0D0 
      COEFF (3, 4) = 19.0D0 / 196.0D0 
!                                                                       
!          A - values (see scheme)                                      
!                                                                       
      COEFF (1, 1) = 11.0D0 / 108.0D0 
      COEFF (1, 3) = 6250.0D0 / 14553.0D0 
      COEFF (1, 4) = - 2197.0D0 / 21168.0D0 
      COEFF (1, 5) = 81.0D0 / 176.0D0 
      COEFF (1, 6) = 171.0D0 / 1960.0D0 
      COEFF (1, 7) = 1.0D0 / 40.0D0 
      RETURN 
      END SUBROUTINE RK547C                         
