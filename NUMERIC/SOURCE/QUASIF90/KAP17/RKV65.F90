      SUBROUTINE RKV65 (M, COEFF, QG) 
!                                                                       
!*****************************************************************      
!                                                                *      
! This subroutine determines the coefficients for the embedding  *      
! formula RKV6(5) (method of Verner).                            *      
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
!             Explicit Runge-Kutta methods with estimates of     *      
!             the local truncation error                         *      
!             page 772-790                                       *      
!             SIAM J. NUMER. ANAL., Vol. 15, No. 4, 1978         *      
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
      COEFF (2, 1) = 1.0D0 / 18.0D0 
      COEFF (3, 1) = 1.0D0 / 6.0D0 
      COEFF (4, 1) = 2.0D0 / 9.0D0 
      COEFF (5, 1) = 2.0D0 / 3.0D0 
      COEFF (6, 1) = 1.0D0 
      COEFF (7, 1) = 8.0D0 / 9.0D0 
      COEFF (8, 1) = 1.0D0 
!                                                                       
!         b - values (see scheme)                                       
!                                                                       
      COEFF (2, 2) = 1.0D0 / 18.0D0 
      COEFF (3, 2) = - 1.0D0 / 12.0D0 
      COEFF (4, 2) = - 2.0D0 / 81.0D0 
      COEFF (5, 2) = 40.0D0 / 33.0D0 
      COEFF (6, 2) = - 369.0D0 / 73.0D0 
      COEFF (7, 2) = - 8716.0D0 / 891.0D0 
      COEFF (8, 2) = 3015.0D0 / 256.0D0 
      COEFF (3, 3) = 1.0D0 / 4.0D0 
      COEFF (4, 3) = 4.0D0 / 27.0D0 
      COEFF (5, 3) = - 4.0D0 / 11.0D0 
      COEFF (6, 3) = 72.0D0 / 73.0D0 
      COEFF (7, 3) = 656.0D0 / 297.0D0 
      COEFF (8, 3) = - 9.0D0 / 4.0D0 
      COEFF (4, 4) = 8.0D0 / 81.0D0 
      COEFF (5, 4) = - 56.0D0 / 11.0D0 
      COEFF (6, 4) = 5380.0D0 / 219.0D0 
      COEFF (7, 4) = 39520.0D0 / 891.0D0 
      COEFF (8, 4) = - 4219.0D0 / 78.0D0 
      COEFF (5, 5) = 54.0D0 / 11.0D0 
      COEFF (6, 5) = - 122.85D2 / 584.0D0 
      COEFF (7, 5) = - 416.0D0 / 11.0D0 
      COEFF (8, 5) = 5985.0D0 / 128.0D0 
      COEFF (6, 6) = 2695.0D0 / 1752.0D0 
      COEFF (7, 6) = 52.0D0 / 27.0D0 
      COEFF (8, 6) = - 539.0D0 / 384.0D0 
      COEFF (8, 8) = 693.0D0 / 3328.0D0 
!                                                                       
!         A tilde values (see scheme)                                   
!                                                                       
      COEFF (2, 3) = 57.0D0 / 640.0D0 
      COEFF (2, 5) = - 16.0D0 / 65.0D0 
      COEFF (2, 6) = 1377.0D0 / 2240.0D0 
      COEFF (2, 7) = 121.0D0 / 320.0D0 
      COEFF (3, 4) = 891.0D0 / 8320.0D0 
      COEFF (3, 5) = 2.0D0 / 35.0D0 
!                                                                       
!          A - values (see scheme)                                      
!                                                                       
      COEFF (1, 1) = 3.0D0 / 80.0D0 
      COEFF (1, 3) = 4.0D0 / 25.0D0 
      COEFF (1, 4) = 243.0D0 / 1120.0D0 
      COEFF (1, 5) = 77.0D0 / 160.0D0 
      COEFF (1, 6) = 73.0D0 / 700.0D0 
      RETURN 
      END SUBROUTINE RKV65                          
