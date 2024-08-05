      SUBROUTINE HIHA5 (M, COEFF, QG) 
!                                                                       
!*****************************************************************      
!                                                                *      
! This subroutine determines the coefficients for the embedding  *      
! formula HIHA5 (method of Higham und Hall).                     *      
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
!  sources  : E. Hairer G. Wanner :                              *      
!             Solving Ordinary Differential Equations II         *      
!             page 28-31                                         *      
!             Springer-Verlag Berlin Heidelberg 1991             *      
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
      COEFF (2, 1) = 2.0D0 / 9.0D0 
      COEFF (3, 1) = 1.0D0 / 3.0D0 
      COEFF (4, 1) = 1.0D0 / 2.0D0 
      COEFF (5, 1) = 3.0D0 / 5.0D0 
      COEFF (6, 1) = 1.0D0 
      COEFF (7, 1) = 1.0D0 
!                                                                       
!         b - values (see scheme)                                       
!                                                                       
      COEFF (2, 2) = 2.0D0 / 9.0D0 
      COEFF (3, 2) = 1.0D0 / 12.0D0 
      COEFF (4, 2) = 1.0D0 / 8.0D0 
      COEFF (5, 2) = 91.0D0 / 500.0D0 
      COEFF (6, 2) = - 11.0D0 / 20.0D0 
      COEFF (7, 2) = 1.0D0 / 12.0D0 
      COEFF (3, 3) = 1.0D0 / 4.0D0 
      COEFF (5, 3) = - 27.0D0 / 100.0D0 
      COEFF (6, 3) = 27.0D0 / 20.0D0 
      COEFF (4, 4) = 3.0D0 / 8.0D0 
      COEFF (5, 4) = 78.0D0 / 125.0D0 
      COEFF (6, 4) = 12.0D0 / 5.0D0 
      COEFF (7, 4) = 27.0D0 / 32.0D0 
      COEFF (5, 5) = 8.0D0 / 125.0D0 
      COEFF (6, 5) = - 36.0D0 / 5.0D0 
      COEFF (7, 5) = - 4.0D0 / 3.0D0 
      COEFF (6, 6) = 5.0D0 
      COEFF (7, 6) = 125.0D0 / 96.0D0 
      COEFF (7, 7) = 5.0D0 / 48.0D0 
!                                                                       
!         A tilde values (see scheme)                                   
!                                                                       
      COEFF (2, 3) = 1.0D0 / 12.0D0 
      COEFF (2, 5) = 27.0D0 / 32.0D0 
      COEFF (2, 6) = - 4.0D0 / 3.0D0 
      COEFF (2, 7) = 125.0D0 / 96.0D0 
      COEFF (3, 4) = 5.0D0 / 48.0D0 
!                                                                       
!          A - values (see scheme)                                      
!                                                                       
                                                                        
      COEFF (1, 1) = 2.0D0 / 15.0D0 
      COEFF (1, 3) = 27.0D0 / 80.0D0 
      COEFF (1, 4) = - 2.0D0 / 15.0D0 
      COEFF (1, 5) = 25.0D0 / 48.0D0 
      COEFF (1, 6) = 1.0D0 / 24.0D0 
      COEFF (1, 7) = 1.0D0 / 10.0D0 
      RETURN 
      END SUBROUTINE HIHA5                          
