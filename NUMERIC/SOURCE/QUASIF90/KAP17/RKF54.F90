      SUBROUTINE RKF54 (M, COEFF, QG) 
!                                                                       
!*****************************************************************      
!                                                                *      
! This subroutine determines the coefficients for the embedding  *      
! formula RKF5(4) (method of Fehlberg).                          *      
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
!             Klassische Runge-Kutta-Formeln vierter und         *      
!             niedrigerer Ordnung mit Schrittweiten-Kontrolle    *      
!             und ihrer Anwendung auf WÑrmeleitungsprobleme      *      
!             page 61-71                                         *      
!             Computing 6 1970                                   *      
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
!         a - values (see scheme)                                       
!                                                                       
      COEFF (2, 1) = 1.0D0 / 4.0D0 
      COEFF (3, 1) = 3.0D0 / 8.0D0 
      COEFF (4, 1) = 12.0D0 / 13.0D0 
      COEFF (5, 1) = 1.0D0 
      COEFF (6, 1) = 1.0D0 / 2.0D0 
!                                                                       
!         b - values (see scheme)                                       
!                                                                       
      COEFF (2, 2) = 1.0D0 / 4.0D0 
      COEFF (3, 2) = 3.0D0 / 32.0D0 
      COEFF (4, 2) = 1932.0D0 / 2197.0D0 
      COEFF (5, 2) = 439.0D0 / 216.0D0 
      COEFF (6, 2) = - 8.0D0 / 27.0D0 
      COEFF (3, 3) = 9.0D0 / 32.0D0 
      COEFF (4, 3) = - 7200.0D0 / 2197.0D0 
      COEFF (5, 3) = - 8.0D0 
      COEFF (6, 3) = 2.0D0 
      COEFF (4, 4) = 7296.0D0 / 2197.0D0 
      COEFF (5, 4) = 3680.0D0 / 513.0D0 
      COEFF (6, 4) = - 3544.0D0 / 2565.0D0 
      COEFF (5, 5) = - 845.0D0 / 4104.0D0 
      COEFF (6, 5) = 1859.0D0 / 4104.0D0 
      COEFF (6, 6) = - 11.0D0 / 40.0D0 
!                                                                       
!         A tilde values (see scheme)                                   
!                                                                       
      COEFF (2, 3) = 16.0D0 / 135.0D0 
      COEFF (2, 5) = 6656.0D0 / 12825.0D0 
      COEFF (2, 6) = 28561.0D0 / 56430.0D0 
      COEFF (3, 4) = - 9.0D0 / 50.0D0 
      COEFF (3, 5) = 2.0D0 / 55.0D0 
!                                                                       
!          A - values (see scheme)                                      
!                                                                       
      COEFF (1, 1) = 25.0D0 / 216.0D0 
      COEFF (1, 3) = 1408.0D0 / 2565.0D0 
      COEFF (1, 4) = 2197.0D0 / 4104.0D0 
      COEFF (1, 5) = - 1.0D0 / 5.0D0 
      RETURN 
      END SUBROUTINE RKF54                          
