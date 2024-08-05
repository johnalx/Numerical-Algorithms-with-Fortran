      SUBROUTINE COEFFI (M, IFLAG, COEFF, QG) 
!                                                                       
!*****************************************************************      
!                                                                *      
! This subroutine determines the coefficients for the designated *      
! method.                                                        *      
!                                                                *      
!                                                                *      
! INPUT PARAMETERS:                                              *      
! =================                                              *      
! M       : Dimension of the matrix COEFF depending on the chosen*      
!           embedding formula                                    *      
! IFLAG   : classifies the various embedding formulas,           *      
!           1 <= IFLAG <= 22                                     *      
!                                                                *      
!                                                                *      
! OUTPUT PARAMETERS:                                             *      
! ==================                                             *      
! COEFF   : 2-dimensional DOUBLE PRECISION array COEFF(1:16,1:M) *      
!           with the coefficients for the embedding formula      *      
! QG      : DOUBLE PRECISION value for the global error order    *      
!                                                                *      
!                                                                *      
! LOCAL VARIABLES:                                               *      
! ================                                               *      
! I       : ] loop                                               *      
! J       : ]    counters                                        *      
!                                                                *      
!                                                                *      
!----------------------------------------------------------------*      
!                                                                *      
!  required subroutines: RKF43, RKF54, RK546M, RKE54, HIHA5,     *      
!                        RK547S, RK547M, RK547C, RK658M,         *      
!                        RK658S, RK658C, RKV65, RKF65A,          *      
!                        RKF65B, RKC65, RKV65A, RKV65B,          *      
!                        RKV76, RK8713M, RKF87, RKV87, RKV98     *      
!                                                                *      
!*****************************************************************      
!                                                                *      
!  Author   : Volker Krger                                      *      
!  Date     : 18.05.1993                                         *      
!  Source   : FORTRAN 77                                         *      
!                                                                *      
!*****************************************************************      
!                                                                       
! Declarations                                                          
!                                                                       
      DOUBLEPRECISION COEFF (16, M), QG 
!                                                                       
! Set the matrix COEFF equal to zero initially                          
!                                                                       
      DO 10 I = 1, M 
         DO 10 J = 1, M 
            COEFF (I, J) = 0.0D0 
   10 CONTINUE 
!                                                                       
! If IFLAG=1, use the method : RKF4(3)                                  
!                                                                       
      IF (IFLAG.EQ.1) THEN 
         CALL RKF43 (M, COEFF, QG) 
!                                                                       
! If IFLAG=2, use the method : RKF5(4)                                  
!                                                                       
      ELSEIF (IFLAG.EQ.2) THEN 
         CALL RKF54 (M, COEFF, QG) 
!                                                                       
! If IFLAG=3, use the method : RK5(4)6M                                 
!                                                                       
      ELSEIF (IFLAG.EQ.3) THEN 
         CALL RK546M (M, COEFF, QG) 
!                                                                       
! If IFLAG=4, use the method : RKE5(4)                                  
!                                                                       
      ELSEIF (IFLAG.EQ.4) THEN 
         CALL RKE54 (M, COEFF, QG) 
!                                                                       
! If IFLAG=5, use the method : RK5(4)7S                                 
!                                                                       
      ELSEIF (IFLAG.EQ.5) THEN 
         CALL HIHA5 (M, COEFF, QG) 
!                                                                       
! If IFLAG=6, use the method : RK5(4)7S                                 
!                                                                       
      ELSEIF (IFLAG.EQ.6) THEN 
         CALL RK547S (M, COEFF, QG) 
!                                                                       
! If IFLAG=7, use the method : RK5(4)7M                                 
!                                                                       
      ELSEIF (IFLAG.EQ.7) THEN 
         CALL RK547M (M, COEFF, QG) 
!                                                                       
! If IFLAG=8, use the method : RK5(4)7C                                 
!                                                                       
      ELSEIF (IFLAG.EQ.8) THEN 
         CALL RK547C (M, COEFF, QG) 
!                                                                       
! If IFLAG=9, use the method : RK6(5)8M                                 
!                                                                       
      ELSEIF (IFLAG.EQ.9) THEN 
         CALL RK658M (M, COEFF, QG) 
!                                                                       
! If IFLAG=10, use the method : RK6(5)8S                                
!                                                                       
      ELSEIF (IFLAG.EQ.10) THEN 
         CALL RK658S (M, COEFF, QG) 
!                                                                       
! If IFLAG=11, use the method : RK6(5)8C                                
!                                                                       
      ELSEIF (IFLAG.EQ.11) THEN 
         CALL RK658C (M, COEFF, QG) 
!                                                                       
! If IFLAG=12, use the method : RKV6(5)                                 
!                                                                       
      ELSEIF (IFLAG.EQ.12) THEN 
         CALL RKV65 (M, COEFF, QG) 
!                                                                       
! If IFLAG=13, use the method : RKF6(5)A                                
!                                                                       
      ELSEIF (IFLAG.EQ.13) THEN 
         CALL RKF65A (M, COEFF, QG) 
!                                                                       
! If IFLAG=14, use the method : RKF6(5)B                                
!                                                                       
      ELSEIF (IFLAG.EQ.14) THEN 
         CALL RKF65B (M, COEFF, QG) 
!                                                                       
! If IFLAG=15, use the method : RKC6(5)                                 
!                                                                       
      ELSEIF (IFLAG.EQ.15) THEN 
         CALL RKC65 (M, COEFF, QG) 
!                                                                       
! If IFLAG=16, use the method : RKV6(5)9A                               
!                                                                       
      ELSEIF (IFLAG.EQ.16) THEN 
         CALL RKV65A (M, COEFF, QG) 
!                                                                       
! If IFLAG=17, use the method : RKV6(5)9B                               
!                                                                       
      ELSEIF (IFLAG.EQ.17) THEN 
         CALL RKV65B (M, COEFF, QG) 
!                                                                       
! If IFLAG=18, use the method : RKV7(6)                                 
!                                                                       
      ELSEIF (IFLAG.EQ.18) THEN 
         CALL RKV76 (M, COEFF, QG) 
!                                                                       
! If IFLAG=19, use the method : RK8(7)13M                               
!                                                                       
      ELSEIF (IFLAG.EQ.19) THEN 
         CALL RK8713 (M, COEFF, QG) 
!                                                                       
! If IFLAG=20, use the method : RKF8(7)                                 
!                                                                       
      ELSEIF (IFLAG.EQ.20) THEN 
         CALL RKF87 (M, COEFF, QG) 
!                                                                       
! If IFLAG=21, use the method : RKV8(7)                                 
!                                                                       
      ELSEIF (IFLAG.EQ.21) THEN 
         CALL RKV87 (M, COEFF, QG) 
!                                                                       
! If IFLAG=22, use the method : RKV9(8)                                 
!                                                                       
      ELSEIF (IFLAG.EQ.22) THEN 
         CALL RKV98 (M, COEFF, QG) 
      ENDIF 
      RETURN 
      END SUBROUTINE COEFFI                         
