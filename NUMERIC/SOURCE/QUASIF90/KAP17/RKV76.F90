      SUBROUTINE RKV76 (M, COEFF, QG) 
!                                                                       
!*****************************************************************      
!                                                                *      
! This subroutine determines the coefficients for the embedding  *      
! formula RKV7(6) (method of Verner).                            *      
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
      QG = 6.0D0 
!                                                                       
! Determine the matrix elements in COEFF                                
!                                                                       
!         a - values (see scheme )                                      
!                                                                       
      COEFF (2, 1) = 1.0D0 / 12.0D0 
      COEFF (3, 1) = 1.0D0 / 6.0D0 
      COEFF (4, 1) = 1.0D0 / 4.0D0 
      COEFF (5, 1) = 3.0D0 / 4.0D0 
      COEFF (6, 1) = 16.0D0 / 17.0D0 
      COEFF (7, 1) = 1.0D0 / 2.0D0 
      COEFF (8, 1) = 1.0D0 
      COEFF (9, 1) = 2.0D0 / 3.0D0 
      COEFF (10, 1) = 1.0D0 
!                                                                       
!         b - values (see scheme)                                       
!                                                                       
      COEFF (2, 2) = 1.0D0 / 12.0D0 
      COEFF (4, 2) = 1.0D0 / 16.0D0 
      COEFF (5, 2) = 21.0D0 / 16.0D0 
      COEFF (6, 2) = 1344.688D3 / 2505.63D2 
      COEFF (7, 2) = - 559.0D0 / 384.0D0 
      COEFF (8, 2) = - 625.0D0 / 224.0D0 
      COEFF (9, 2) = - 1225.3D1 / 9914.4D1 
      COEFF (10, 2) = 3051.7D1 / 2512.0D0 
      COEFF (3, 3) = 1.0D0 / 6.0D0 
      COEFF (4, 4) = 3.0D0 / 16.0D0 
      COEFF (5, 4) = - 81.0D0 / 16.0D0 
      COEFF (6, 4) = - 1709.184D3 / 8352.1D1 
      COEFF (7, 4) = 6.0D0 
      COEFF (8, 4) = 12.0D0 
      COEFF (9, 4) = 16.0D0 / 27.0D0 
      COEFF (10, 4) = - 7296.0D0 / 157.0D0 
      COEFF (5, 5) = 9.0D0 / 2.0D0 
      COEFF (6, 5) = 1365.632D3 / 8352.1D1 
      COEFF (7, 5) = - 204.0D0 / 47.0D0 
      COEFF (8, 5) = - 456.0D0 / 47.0D0 
      COEFF (9, 5) = 16.0D0 / 459.0D0 
      COEFF (10, 5) = 2687.28D2 / 7379.0D0 
      COEFF (6, 6) = - 7820.8D1 / 2505.63D2 
      COEFF (7, 6) = 14.0D0 / 39.0D0 
      COEFF (8, 6) = 48.0D0 / 91.0D0 
      COEFF (9, 6) = 2907.2D1 / 1611.09D2 
      COEFF (10, 6) = 2472.0D0 / 2041.0D0 
      COEFF (7, 7) = - 4913.0D0 / 7820.8D1 
      COEFF (8, 7) = 1473.9D1 / 1368.64D2 
      COEFF (9, 7) = - 2023.0D0 / 7581.6D1 
      COEFF (10, 7) = - 3522.621D3 / 1074.3824D4 
      COEFF (8, 8) = 6.0D0 / 7.0D0 
      COEFF (9, 8) = 112.0D0 / 1239.3D1 
      COEFF (10, 8) = 132.0D0 / 157.0D0 
      COEFF (10, 10) = - 1239.3D1 / 4396.0D0 
!                                                                       
!         A tilde values (see scheme)                                   
!                                                                       
      COEFF (2, 3) = 2881.0D0 / 4032.0D1 
      COEFF (2, 6) = 1216.0D0 / 2961.0D0 
      COEFF (2, 7) = - 2624.0D0 / 4095.0D0 
      COEFF (2, 8) = 2413.7569D4 / 5748.2880D4 
      COEFF (2, 9) = - 4.0D0 / 21.0D0 
      COEFF (3, 4) = 4131.0D0 / 3920.0D0 
      COEFF (3, 5) = - 157.0D0 / 1260.0D0 
!                                                                       
!          A - values (see scheme)                                      
!                                                                       
      COEFF (1, 1) = 7.0D0 / 90.0D0 
      COEFF (1, 4) = 16.0D0 / 45.0D0 
      COEFF (1, 5) = 16.0D0 / 45.0D0 
      COEFF (1, 7) = 2.0D0 / 15.0D0 
      COEFF (1, 8) = 7.0D0 / 90.0D0 
      RETURN 
      END SUBROUTINE RKV76                          
