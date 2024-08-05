      SUBROUTINE RK8713 (M, COEFF, QG) 
!                                                                       
!*****************************************************************      
!                                                                *      
! This subroutine determines the coefficients for the embedding  *      
! formula RK8(7)13M (method of Dormand und Prince).              *      
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
!  Sources  : J. R. Dormand P. J. Prince :                       *      
!             High order embedded Runge-Kutta formulae           *      
!             page 67-75                                         *      
!             Journal of Computational and Applied Mathematics,  *      
!             volume 7, no 1, 1981                               *      
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
!         a - values (see scheme )                                      
!                                                                       
      COEFF (2, 1) = 1.0D0 / 18.0D0 
      COEFF (3, 1) = 1.0D0 / 12.0D0 
      COEFF (4, 1) = 1.0D0 / 8.0D0 
      COEFF (5, 1) = 5.0D0 / 16.0D0 
      COEFF (6, 1) = 3.0D0 / 8.0D0 
      COEFF (7, 1) = 59.0D0 / 400.0D0 
      COEFF (8, 1) = 93.0D0 / 200.0D0 
      COEFF (9, 1) = 549002.3248D4 / 971916.9821D4 
      COEFF (10, 1) = 13.0D0 / 20.0D0 
      COEFF (11, 1) = 120114.6811D4 / 129901.9798D4 
      COEFF (12, 1) = 1.0D0 
      COEFF (13, 1) = 1.0D0 
!                                                                       
!         b - values (see scheme)                                       
!                                                                       
      COEFF (2, 2) = 1.0D0 / 18.0D0 
      COEFF (3, 2) = 1.0D0 / 48.0D0 
      COEFF (4, 2) = 1.0D0 / 32.0D0 
      COEFF (5, 2) = 5.0D0 / 16.0D0 
      COEFF (6, 2) = 3.0D0 / 80.0D0 
      COEFF (7, 2) = 2944.3841D4 / 61456.3906D4 
      COEFF (8, 2) = 1601.6141D4 / 94669.2911D4 
      COEFF (9, 2) = 3963.2708D4 / 57359.1083D4 
      COEFF (10, 2) = 24612.1993D4 / 134084.7787D4 
      COEFF (11, 2) = - 102846.8189D4 / 84618.0014D4 
      COEFF (12, 2) = 18589.2177D4 / 71811.6043D4 
      COEFF (13, 2) = 40386.3854D4 / 49106.3109D4 
      COEFF (3, 3) = 1.0D0 / 16.0D0 
      COEFF (4, 4) = 3.0D0 / 32.0D0 
      COEFF (5, 4) = - 75.0D0 / 64.0D0 
      COEFF (5, 5) = 75.0D0 / 64.0D0 
      COEFF (6, 5) = 3.0D0 / 16.0D0 
      COEFF (7, 5) = 7773.6538D4 / 69253.8347D4 
      COEFF (8, 5) = 6156.4180D4 / 15873.2637D4 
      COEFF (9, 5) = - 43363.6366D4 / 68370.1615D4 
      COEFF (10, 5) = - 376950.42795D5 / 152687.66246D5 
      COEFF (11, 5) = 847823.5783D4 / 50851.2852D4 
      COEFF (12, 5) = - 318509.4517D4 / 66710.7341D4 
      COEFF (13, 5) = - 506849.2393D4 / 43474.0067D4 
      COEFF (6, 6) = 3.0D0 / 20.0D0 
      COEFF (7, 6) = - 2869.3883D4 / 112500.0000D4 
      COEFF (8, 6) = 2278.9713D4 / 63344.5777D4 
      COEFF (9, 6) = - 42173.9975D4 / 261629.2301D4 
      COEFF (10, 6) = - 30912.1744D4 / 106122.7803D4 
      COEFF (11, 6) = 131172.9495D4 / 143242.2823D4 
      COEFF (12, 6) = - 47775.5414D4 / 109805.3517D4 
      COEFF (13, 6) = - 41142.1997D4 / 54304.3805D4 
      COEFF (7, 7) = 2312.4283D4 / 180000.0000D4 
      COEFF (8, 7) = 54581.5736D4 / 277105.7229D4 
      COEFF (9, 7) = 10030.2831D4 / 72342.3059D4 
      COEFF (10, 7) = - 1299.2083D4 / 49076.6935D4 
      COEFF (11, 7) = - 103041.29995D5 / 17013.04382D5 
      COEFF (12, 7) = - 70363.5378D4 / 23073.9211D4 
      COEFF (13, 7) = 65278.3627D4 / 91429.6604D4 
      COEFF (8, 8) = - 18019.3667D4 / 104330.7555D4 
      COEFF (9, 8) = 79020.4164D4 / 83981.3087D4 
      COEFF (10, 8) = 600594.3493D4 / 210894.7869D4 
      COEFF (11, 8) = - 487779.25059D5 / 30479.39560D5 
      COEFF (12, 8) = 573156.6787D4 / 102754.5527D4 
      COEFF (13, 8) = 111739.62825D5 / 9253.20556D5 
      COEFF (9, 9) = 80063.5310D4 / 378307.1287D4 
      COEFF (10, 9) = 39300.6217D4 / 139667.3457D4 
      COEFF (11, 9) = 153367.26248D5 / 10328.24649D5 
      COEFF (12, 9) = 523286.6602D4 / 85006.6563D4 
      COEFF (13, 9) = - 131589.90841D5 / 61847.27034D5 
      COEFF (10, 10) = 12387.2331D4 / 100102.9789D4 
      COEFF (11, 10) = - 454428.68181D5 / 33984.67696D5 
      COEFF (12, 10) = - 409366.4535D4 / 80868.8257D4 
      COEFF (13, 10) = 393664.7629D4 / 197804.9680D4 
      COEFF (11, 11) = 306599.3473D4 / 59717.2653D4 
      COEFF (12, 11) = 396213.7247D4 / 180595.7418D4 
      COEFF (13, 11) = - 16052.8059D4 / 68517.8525D4 
      COEFF (12, 12) = 6568.6358D4 / 48791.0083D4 
      COEFF (13, 12) = 24863.8103D4 / 141353.1060D4 
!                                                                       
!         A tilde values (see scheme)                                   
!                                                                       
      COEFF (2, 3) = 1400.5451D4 / 33548.0064D4 
      COEFF (2, 8) = - 5923.8493D4 / 106827.7825D4 
      COEFF (2, 9) = 18160.6767D4 / 75886.7731D4 
      COEFF (2, 10) = 56129.2985D4 / 79784.5732D4 
      COEFF (2, 11) = - 104189.1430D4 / 137134.3529D4 
      COEFF (2, 12) = 76041.7239D4 / 115116.5299D4 
      COEFF (2, 13) = 11882.0643D4 / 75113.8087D4 
      COEFF (3, 4) = - 52874.7749D4 / 222060.7170D4 
      COEFF (3, 5) = 1.0D0 / 4.0D0 
!                                                                       
!          A - values (see scheme)                                      
!                                                                       
      COEFF (1, 1) = 1345.1932D4 / 45517.6623D4 
      COEFF (1, 6) = - 80871.9846D4 / 97600.0145D4 
      COEFF (1, 7) = 175700.4468D4 / 564515.9321D4 
      COEFF (1, 8) = 65604.5339D4 / 26589.1186D4 
      COEFF (1, 9) = - 386757.4721D4 / 151851.7206D4 
      COEFF (1, 10) = 46588.5868D4 / 32273.6535D4 
      COEFF (1, 11) = 5301.1238D4 / 66751.6719D4 
      COEFF (1, 12) = 2.0D0 / 45.0D0 
      RETURN 
      END SUBROUTINE RK8713                         
