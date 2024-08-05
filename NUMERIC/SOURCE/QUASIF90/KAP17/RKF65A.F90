      SUBROUTINE RKF65A (M, COEFF, QG) 
!                                                                       
!*****************************************************************      
!                                                                *      
! This subroutine determines the coefficients for the embedding  *      
! formula RKF6(5)A (method of Fehlberg).                          *     
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
      QG = 5.0D0 
!                                                                       
! Determine the matrix elements in COEFF                                
!                                                                       
!         a - values (see scheme)                                       
!                                                                       
      COEFF (2, 1) = 1.0D0 / 6.0D0 
      COEFF (3, 1) = 4.0D0 / 15.0D0 
      COEFF (4, 1) = 2.0D0 / 3.0D0 
      COEFF (5, 1) = 4.0D0 / 5.0D0 
      COEFF (6, 1) = 1.0D0 
      COEFF (8, 1) = 1.0D0 
!                                                                       
!         b - values (see scheme)                                       
!                                                                       
      COEFF (2, 2) = 1.0D0 / 6.0D0 
      COEFF (3, 2) = 4.0D0 / 75.0D0 
      COEFF (4, 2) = 5.0D0 / 6.0D0 
      COEFF (5, 2) = - 8.0D0 / 5.0D0 
      COEFF (6, 2) = 361.0D0 / 320.0D0 
      COEFF (7, 2) = - 11.0D0 / 640.0D0 
      COEFF (8, 2) = 93.0D0 / 640.0D0 
      COEFF (3, 3) = 16.0D0 / 75.0D0 
      COEFF (4, 3) = - 8.0D0 / 3.0D0 
      COEFF (5, 3) = 144.0D0 / 25.0D0 
      COEFF (6, 3) = - 18.0D0 / 5.0D0 
      COEFF (8, 3) = - 18.0D0 / 5.0D0 
      COEFF (4, 4) = 5.0D0 / 2.0D0 
      COEFF (5, 4) = - 4.0D0 
      COEFF (6, 4) = 407.0D0 / 128.0D0 
      COEFF (7, 4) = 11.0D0 / 256.0D0 
      COEFF (8, 4) = 803.0D0 / 256.0D0 
      COEFF (5, 5) = 16.0D0 / 25.0D0 
      COEFF (6, 5) = - 11.0D0 / 80.0D0 
      COEFF (7, 5) = - 11.0D0 / 160.0D0 
      COEFF (8, 5) = - 11.0D0 / 160.0D0 
      COEFF (6, 6) = 55.0D0 / 128.0D0 
      COEFF (7, 6) = 11.0D0 / 256.0D0 
      COEFF (8, 6) = 99.0D0 / 256.0D0 
      COEFF (8, 8) = 1.0D0 
!                                                                       
!         A tilde values (see scheme)                                   
!                                                                       
      COEFF (2, 3) = 7.0D0 / 1408.0D0 
      COEFF (2, 5) = 1125.0D0 / 2816.0D0 
      COEFF (2, 6) = 9.0D0 / 32.0D0 
      COEFF (2, 7) = 125.0D0 / 768.0D0 
      COEFF (3, 4) = 5.0D0 / 66.0D0 
      COEFF (3, 5) = 5.0D0 / 66.0D0 
!                                                                       
!          A - values (see scheme)                                      
!                                                                       
      COEFF (1, 1) = 31.0D0 / 384.0D0 
      COEFF (1, 3) = 1125.0D0 / 2816.0D0 
      COEFF (1, 4) = 9.0D0 / 32.0D0 
      COEFF (1, 5) = 125.0D0 / 768.0D0 
      COEFF (1, 6) = 5.0D0 / 66.0D0 
      RETURN 
      END SUBROUTINE RKF65A                         
