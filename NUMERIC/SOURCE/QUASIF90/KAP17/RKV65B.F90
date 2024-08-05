      SUBROUTINE RKV65B (M, COEFF, QG) 
!                                                                       
!*****************************************************************      
!                                                                *      
! This subroutine determines the coefficients for the embedding  *      
! formula RKV6(5)9B (method of Verner).                          *      
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
!  Sources  : J. H. Verner :                                     *      
!             Some Runge-Kutta formular pairs                    *      
!             page 496-511                                       *      
!             SIAM J. NUMER. ANAL., Vol. 28, No. 2, 1991         *      
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
      COEFF (2, 1) = 1.0D0 / 8.0D0 
      COEFF (3, 1) = 1.0D0 / 6.0D0 
      COEFF (4, 1) = 1.0D0 / 4.0D0 
      COEFF (5, 1) = 1.0D0 / 2.0D0 
      COEFF (6, 1) = 3.0D0 / 5.0D0 
      COEFF (7, 1) = 4.0D0 / 5.0D0 
      COEFF (8, 1) = 1.0D0 
      COEFF (9, 1) = 1.0D0 
!                                                                       
!         b - values (see scheme)                                       
!                                                                       
      COEFF (2, 2) = 1.0D0 / 8.0D0 
      COEFF (3, 2) = 1.0D0 / 18.0D0 
      COEFF (4, 2) = 1.0D0 / 16.0D0 
      COEFF (5, 2) = 1.0D0 / 4.0D0 
      COEFF (6, 2) = 134.0D0 / 625.0D0 
      COEFF (7, 2) = - 98.0D0 / 1875.0D0 
      COEFF (8, 2) = 9.0D0 / 50.0D0 
      COEFF (9, 2) = 11.0D0 / 144.0D0 
      COEFF (3, 3) = 1.0D0 / 9.0D0 
      COEFF (4, 4) = 3.0D0 / 16.0D0 
      COEFF (5, 4) = - 3.0D0 / 4.0D0 
      COEFF (6, 4) = - 333.0D0 / 625.0D0 
      COEFF (7, 4) = 12.0D0 / 625.0D0 
      COEFF (8, 4) = 21.0D0 / 25.0D0 
      COEFF (5, 5) = 1.0D0 
      COEFF (6, 5) = 476.0D0 / 625.0D0 
      COEFF (7, 5) = 107.36D2 / 131.25D2 
      COEFF (8, 5) = - 2924.0D0 / 1925.0D0 
      COEFF (9, 5) = 256.0D0 / 693.0D0 
      COEFF (6, 6) = 98.0D0 / 625.0D0 
      COEFF (7, 6) = - 1936.0D0 / 1875.0D0 
      COEFF (8, 6) = 74.0D0 / 25.0D0 
      COEFF (7, 7) = 22.0D0 / 21.0D0 
      COEFF (8, 7) = - 15.0D0 / 7.0D0 
      COEFF (9, 7) = 125.0D0 / 504.0D0 
      COEFF (8, 8) = 15.0D0 / 22.0D0 
      COEFF (9, 8) = 125.0D0 / 528.0D0 
      COEFF (9, 9) = 5.0D0 / 72.0D0 
!                                                                       
!         A tilde values (see scheme)                                   
!                                                                       
      COEFF (2, 3) = 11.0D0 / 144.0D0 
      COEFF (2, 6) = 256.0D0 / 693.0D0 
      COEFF (2, 8) = 125.0D0 / 504.0D0 
      COEFF (2, 9) = 125.0D0 / 528.0D0 
      COEFF (3, 4) = 5.0D0 / 72.0D0 
!                                                                       
!          A - values (see scheme)                                      
!                                                                       
      COEFF (1, 1) = 1.0D0 / 18.0D0 
      COEFF (1, 4) = 32.0D0 / 63.0D0 
      COEFF (1, 5) = - 2.0D0 / 3.0D0 
      COEFF (1, 6) = 125.0D0 / 126.0D0 
      COEFF (1, 8) = - 5.0D0 / 63.0D0 
      COEFF (1, 9) = 4.0D0 / 21.0D0 
      RETURN 
      END SUBROUTINE RKV65B                         
