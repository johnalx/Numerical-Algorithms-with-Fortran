      SUBROUTINE RKF87 (M, COEFF, QG) 
!                                                                       
!*****************************************************************      
!                                                                *      
! This subroutine determines the coefficients for the embedding  *      
! formula RKF8(7) (method of Fehlberg).                          *      
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
!  Sources  : E. Fehlberg :                                      *      
!             Klassische Runge-Kutta-Formeln fÅnfter und         *      
!             siebenter Ordnung mit Schrittweiten-Kontrolle      *      
!             page 93-106                                        *      
!             Computing 4 1969                                   *      
!                                                                *      
!*****************************************************************      
!                                                                       
! Declarations                                                          
!                                                                       
      DOUBLEPRECISION COEFF (16, M), QG 
!                                                                       
! Initialize QG                                                         
!                                                                       
      QG = 7.0D0 
!                                                                       
! Determine the matrix elements in COEFF                                
!                                                                       
!         a - values (see scheme)                                       
!                                                                       
      COEFF (2, 1) = 2.0D0 / 27.0D0 
      COEFF (3, 1) = 1.0D0 / 9.0D0 
      COEFF (4, 1) = 1.0D0 / 6.0D0 
      COEFF (5, 1) = 5.0D0 / 12.0D0 
      COEFF (6, 1) = 1.0D0 / 2.0D0 
      COEFF (7, 1) = 5.0D0 / 6.0D0 
      COEFF (8, 1) = 1.0D0 / 6.0D0 
      COEFF (9, 1) = 2.0D0 / 3.0D0 
      COEFF (10, 1) = 1.0D0 / 3.0D0 
      COEFF (11, 1) = 1.0D0 
      COEFF (13, 1) = 1.0D0 
!                                                                       
!         b - values (see scheme)                                       
!                                                                       
      COEFF (2, 2) = 2.0D0 / 27.0D0 
      COEFF (3, 2) = 1.0D0 / 36.0D0 
      COEFF (4, 2) = 1.0D0 / 24.0D0 
      COEFF (5, 2) = 5.0D0 / 12.0D0 
      COEFF (6, 2) = 1.0D0 / 20.0D0 
      COEFF (7, 2) = - 25.0D0 / 108.0D0 
      COEFF (8, 2) = 31.0D0 / 300.0D0 
      COEFF (9, 2) = 2.0D0 
      COEFF (10, 2) = - 91.0D0 / 108.0D0 
      COEFF (11, 2) = 2383.0D0 / 4100.0D0 
      COEFF (12, 2) = 3.0D0 / 205.0D0 
      COEFF (13, 2) = - 1777.0D0 / 4100.0D0 
      COEFF (3, 3) = 1.0D0 / 12.0D0 
      COEFF (4, 4) = 1.0D0 / 8.0D0 
      COEFF (5, 4) = - 25.0D0 / 16.0D0 
      COEFF (5, 5) = 25.0D0 / 16.0D0 
      COEFF (6, 5) = 1.0D0 / 4.0D0 
      COEFF (7, 5) = 125.0D0 / 108.0D0 
      COEFF (9, 5) = - 53.0D0 / 6.0D0 
      COEFF (10, 5) = 23.0D0 / 108.0D0 
      COEFF (11, 5) = - 341.0D0 / 164.0D0 
      COEFF (13, 5) = - 341.0D0 / 164.0D0 
      COEFF (6, 6) = 1.0D0 / 5.0D0 
      COEFF (7, 6) = - 65.0D0 / 27.0D0 
      COEFF (8, 6) = 61.0D0 / 225.0D0 
      COEFF (9, 6) = 704.0D0 / 45.0D0 
      COEFF (10, 6) = - 976.0D0 / 135.0D0 
      COEFF (11, 6) = 4496.0D0 / 1025.0D0 
      COEFF (13, 6) = 4496.0D0 / 1025.0D0 
      COEFF (7, 7) = 125.0D0 / 54.0D0 
      COEFF (8, 7) = - 2.0D0 / 9.0D0 
      COEFF (9, 7) = - 107.0D0 / 9.0D0 
      COEFF (10, 7) = 311.0D0 / 54.0D0 
      COEFF (11, 7) = - 301.0D0 / 82.0D0 
      COEFF (12, 7) = - 6.0D0 / 41.0D0 
      COEFF (13, 7) = - 289.0D0 / 82.0D0 
      COEFF (8, 8) = 13.0D0 / 900.0D0 
      COEFF (9, 8) = 67.0D0 / 90.0D0 
      COEFF (10, 8) = - 19.0D0 / 60.0D0 
      COEFF (11, 8) = 2133.0D0 / 4100.0D0 
      COEFF (12, 8) = - 3.0D0 / 205.0D0 
      COEFF (13, 8) = 2193.0D0 / 4100.0D0 
      COEFF (9, 9) = 3.0D0 
      COEFF (10, 9) = 17.0D0 / 6.0D0 
      COEFF (11, 9) = 45.0D0 / 82.0D0 
      COEFF (12, 9) = - 3.0D0 / 41.0D0 
      COEFF (13, 9) = 51.0D0 / 82.0D0 
      COEFF (10, 10) = - 1.0D0 / 12.0D0 
      COEFF (11, 10) = 45.0D0 / 164.0D0 
      COEFF (12, 10) = 3.0D0 / 41.0D0 
      COEFF (13, 10) = 33.0D0 / 164.0D0 
      COEFF (11, 11) = 18.0D0 / 41.0D0 
      COEFF (12, 11) = 6.0D0 / 41.0D0 
      COEFF (13, 11) = 12.0D0 / 41.0D0 
      COEFF (13, 13) = 1.0D0 
!                                                                       
!         A tilde values (see scheme)                                   
!                                                                       
      COEFF (2, 8) = 34.0D0 / 105.0D0 
      COEFF (2, 9) = 9.0D0 / 35.0D0 
      COEFF (2, 10) = 9.0D0 / 35.0D0 
      COEFF (2, 11) = 9.0D0 / 280.0D0 
      COEFF (2, 12) = 9.0D0 / 280.0D0 
      COEFF (3, 4) = 41.0D0 / 840.0D0 
      COEFF (3, 5) = 41.0D0 / 840.0D0 
!                                                                       
!          A - values (see scheme)                                      
!                                                                       
      COEFF (1, 1) = 41.0D0 / 840.0D0 
      COEFF (1, 6) = 34.0D0 / 105.0D0 
      COEFF (1, 7) = 9.0D0 / 35.0D0 
      COEFF (1, 8) = 9.0D0 / 35.0D0 
      COEFF (1, 9) = 9.0D0 / 280.0D0 
      COEFF (1, 10) = 9.0D0 / 280.0D0 
      COEFF (1, 11) = 41.0D0 / 840.0D0 
      RETURN 
      END SUBROUTINE RKF87                          
