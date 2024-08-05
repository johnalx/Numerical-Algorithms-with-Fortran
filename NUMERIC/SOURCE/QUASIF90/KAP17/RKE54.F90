      SUBROUTINE RKE54 (M, COEFF, QG) 
!                                                                       
!*****************************************************************      
!                                                                *      
! This subroutine determines the coefficients for the embedding  *      
! formula RKE5(4) (method of England).                           *      
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
!  Sources  : G. Engeln-MÅllges F. Reuter :                      *      
!             Formelsammlung zur Numerischen Mathematik mit      *      
!             Standard FORTRAN 77 - Programmen                   *      
!             page 334                                           *      
!             Wissenschaftsverlag  1988                          *      
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
      COEFF (2, 1) = 1.0D0 / 2.0D0 
      COEFF (3, 1) = 1.0D0 / 2.0D0 
      COEFF (4, 1) = 1.0D0 
      COEFF (5, 1) = 2.0D0 / 3.0D0 
      COEFF (6, 1) = 1.0D0 / 5.0D0 
!                                                                       
!         b - values (see scheme)                                       
!                                                                       
      COEFF (2, 2) = 1.0D0 / 2.0D0 
      COEFF (3, 2) = 1.0D0 / 4.0D0 
      COEFF (5, 2) = 7.0D0 / 27.0D0 
      COEFF (6, 2) = 28.0D0 / 625.0D0 
      COEFF (3, 3) = 1.0D0 / 4.0D0 
      COEFF (4, 3) = - 1.0D0 
      COEFF (5, 3) = 10.0D0 / 27.0D0 
      COEFF (6, 3) = - 125.0D0 / 625.0D0 
      COEFF (4, 4) = 2.0D0 
      COEFF (6, 4) = 546.0D0 / 625.0D0 
      COEFF (5, 5) = 1.0D0 / 27.0D0 
      COEFF (6, 5) = 54.0D0 / 625.0D0 
      COEFF (6, 6) = - 378.0D0 / 625.0D0 
!                                                                       
!         A tilde values (see scheme)                                   
!                                                                       
      COEFF (2, 3) = 14.0D0 / 336.0D0 
      COEFF (2, 6) = 35.0D0 / 336.0D0 
      COEFF (3, 4) = 162.0D0 / 336.0D0 
      COEFF (3, 5) = 125.0D0 / 336.0D0 
!                                                                       
!          A - values (see scheme)                                      
!                                                                       
      COEFF (1, 1) = 1.0D0 / 6.0D0 
      COEFF (1, 3) = 4.0D0 / 6.0D0 
      COEFF (1, 4) = 1.0D0 / 6.0D0 
      RETURN 
      END SUBROUTINE RKE54                          
