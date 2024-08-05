      SUBROUTINE RKF43 (M, COEFF, QG) 
!                                                                       
!*****************************************************************      
!                                                                *      
! This subroutine determines the coefficients for the embedding  *      
! formula RKF4(3) (method of Fehlberg).                          *      
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
      QG = 3.0D0 
!                                                                       
! Determine the matrix elements in COEFF                                
!                                                                       
!         a - values (refer to the coefficient schemes in chapter 17.3.4
!                    ( similarly below )                                
!                                                                       
      COEFF (2, 1) = 2.0D0 / 7.0D0 
      COEFF (3, 1) = 7.0D0 / 15.0D0 
      COEFF (4, 1) = 35.0D0 / 38.0D0 
      COEFF (5, 1) = 1.0D0 
!                                                                       
!         b - values (see scheme)                                       
!                                                                       
      COEFF (2, 2) = 2.0D0 / 7.0D0 
      COEFF (3, 2) = 77.0D0 / 900.0D0 
      COEFF (4, 2) = 805.0D0 / 1444.0D0 
      COEFF (5, 2) = 79.0D0 / 490.0D0 
      COEFF (3, 3) = 343.0D0 / 900.0D0 
      COEFF (4, 3) = - 77175.0D0 / 54872.0D0 
      COEFF (4, 4) = 97125.0D0 / 54872.0D0 
      COEFF (5, 4) = 2175.0D0 / 3626.0D0 
      COEFF (5, 5) = 2166.0D0 / 9065.0D0 
!                                                                       
!         A tilde values (see scheme)                                   
!                                                                       
      COEFF (2, 3) = 79.0D0 / 490.0D0 
      COEFF (2, 5) = 2175.0D0 / 3626.0D0 
      COEFF (3, 4) = 2166.0D0 / 9065.0D0 
!                                                                       
!          A - values (see scheme)                                      
!                                                                       
      COEFF (1, 1) = 229.0D0 / 1470.0D0 
      COEFF (1, 3) = 1125.0D0 / 1813.0D0 
      COEFF (1, 4) = 13718.0D0 / 81585.0D0 
      COEFF (1, 5) = 1.0D0 / 18.0D0 
      RETURN 
      END SUBROUTINE RKF43                          
