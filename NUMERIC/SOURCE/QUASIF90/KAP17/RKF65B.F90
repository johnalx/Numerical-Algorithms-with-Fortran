      SUBROUTINE RKF65B (M, COEFF, QG) 
!                                                                       
!*****************************************************************      
!                                                                *      
! This subroutine determines the coefficients for the embedding  *      
! formula RKF6(5)B (method of Fehlberg).                          *     
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
!             Some old and new Runge-Kutta formulas with         *      
!             stepsize control and their error coefficients      *      
!             page 265-270                                       *      
!             Computing 34 1985                                  *      
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
      COEFF (5, 1) = 65.0D0 / 81.0D0 
      COEFF (6, 1) = 1.0D0 
      COEFF (7, 1) = 1.0D0 / 15.0D0 
      COEFF (8, 1) = 1.0D0 
!                                                                       
!         b - values (see scheme)                                       
!                                                                       
      COEFF (2, 2) = 1.0D0 / 6.0D0 
      COEFF (3, 2) = 4.0D0 / 75.0D0 
      COEFF (4, 2) = 5.0D0 / 6.0D0 
      COEFF (5, 2) = - 43745.0D0 / 26244.0D0 
      COEFF (6, 2) = 16539.0D0 / 13780.0D0 
      COEFF (7, 2) = - 409063.0D0 / 195000.0D0 
      COEFF (8, 2) = 231213.0D0 / 193180.0D0 
      COEFF (3, 3) = 16.0D0 / 75.0D0 
      COEFF (4, 3) = - 8.0D0 / 3.0D0 
      COEFF (5, 3) = 353860.0D0 / 59049.0D0 
      COEFF (6, 3) = - 204.0D0 / 53.0D0 
      COEFF (7, 3) = 124.0D0 / 75.0D0 
      COEFF (8, 3) = - 2820.0D0 / 743.0D0 
      COEFF (4, 4) = 5.0D0 / 2.0D0 
      COEFF (5, 4) = - 493675.0D0 / 118098.0D0 
      COEFF (6, 4) = 232595.0D0 / 69006.0D0 
      COEFF (7, 4) = 25357.0D0 / 8680.0D0 
      COEFF (8, 4) = 9512.525D3 / 2902.158D3 
      COEFF (5, 5) = 155155.0D0 / 236196.0D0 
      COEFF (6, 5) = - 91.0D0 / 636.0D0 
      COEFF (7, 5) = - 8928.0D0 / 1375.0D0 
      COEFF (8, 5) = - 5113.0D0 / 80244.0D0 
      COEFF (6, 6) = 314928.0D0 / 747565.0D0 
      COEFF (7, 6) = 79184.709D3 / 19394.375D3 
      COEFF (8, 6) = 584348.904D3 / 156152.2235D4 
      COEFF (8, 8) = 30800.0D0 / 2989.089D3 
!                                                                       
!         A tilde values (see scheme)                                   
!                                                                       
      COEFF (2, 3) = 167.0D0 / 2080.0D0 
      COEFF (2, 5) = 91375.0D0 / 229152.0D0 
      COEFF (2, 6) = 41.0D0 / 144.0D0 
      COEFF (2, 7) = 43046.721D3 / 269010.560D3 
      COEFF (3, 4) = 125.0D0 / 150192.0D0 
      COEFF (3, 5) = 743.0D0 / 9856.0D0 
!                                                                       
!          A - values (see scheme)                                      
!                                                                       
      COEFF (1, 1) = 21.0D0 / 260.0D0 
      COEFF (1, 3) = 7625.0D0 / 19096.0D0 
      COEFF (1, 4) = 25.0D0 / 88.0D0 
      COEFF (1, 5) = 1594.323D3 / 9929.920D3 
      COEFF (1, 6) = 53.0D0 / 704.0D0 
      RETURN 
      END SUBROUTINE RKF65B                         
