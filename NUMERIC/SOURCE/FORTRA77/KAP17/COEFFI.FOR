C[BA*)
C[LE*)
      SUBROUTINE COEFFI (M,IFLAG,COEFF,QG)
C[IX{COEFFI}*)
C
C*****************************************************************
C                                                                *
C This subroutine determines the coefficients for the designated *
C method.                                                        *
C[BE*)
C                                                                *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C M       : Dimension of the matrix COEFF depending on the chosen*
C           embedding formula                                    *
C IFLAG   : classifies the various embedding formulas,           *
C           1 <= IFLAG <= 22                                     *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C COEFF   : 2-dimensional DOUBLE PRECISION array COEFF(1:16,1:M) *
C           with the coefficients for the embedding formula      *
C QG      : DOUBLE PRECISION value for the global error order    *
C                                                                *
C                                                                *
C LOCAL VARIABLES:                                               *
C ================                                               *
C I       : ] loop                                               *
C J       : ]    counters                                        *
C                                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  required subroutines: RKF43, RKF54, RK546M, RKE54, HIHA5,     *
C                        RK547S, RK547M, RK547C, RK658M,         *
C                        RK658S, RK658C, RKV65, RKF65A,          *
C                        RKF65B, RKC65, RKV65A, RKV65B,          *
C                        RKV76, RK8713M, RKF87, RKV87, RKV98     *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : Volker Kr�ger                                      *
C  Date     : 18.05.1993                                         *
C  Source   : FORTRAN 77                                         *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C Declarations
C
      DOUBLE PRECISION COEFF(16,M),QG
C
C Set the matrix COEFF equal to zero initially
C
      DO 10 I=1,M
       DO 10 J=1,M
        COEFF(I,J)=0.0D0
10    CONTINUE
C
C If IFLAG=1, use the method : RKF4(3)
C
      IF(IFLAG .EQ. 1) THEN
        CALL RKF43(M,COEFF,QG)
C
C If IFLAG=2, use the method : RKF5(4)
C
      ELSEIF(IFLAG .EQ. 2) THEN
        CALL RKF54(M,COEFF,QG)
C
C If IFLAG=3, use the method : RK5(4)6M
C
      ELSEIF(IFLAG .EQ. 3) THEN
        CALL RK546M(M,COEFF,QG)
C
C If IFLAG=4, use the method : RKE5(4)
C
      ELSEIF(IFLAG .EQ. 4) THEN
        CALL RKE54(M,COEFF,QG)
C
C If IFLAG=5, use the method : RK5(4)7S
C
      ELSEIF(IFLAG .EQ. 5) THEN
        CALL HIHA5(M,COEFF,QG)
C
C If IFLAG=6, use the method : RK5(4)7S
C
      ELSEIF(IFLAG .EQ. 6) THEN
        CALL RK547S(M,COEFF,QG)
C
C If IFLAG=7, use the method : RK5(4)7M
C
      ELSEIF(IFLAG .EQ. 7) THEN
        CALL RK547M(M,COEFF,QG)
C
C If IFLAG=8, use the method : RK5(4)7C
C
      ELSEIF(IFLAG .EQ. 8) THEN
        CALL RK547C(M,COEFF,QG)
C
C If IFLAG=9, use the method : RK6(5)8M
C
      ELSEIF(IFLAG .EQ. 9) THEN
        CALL RK658M(M,COEFF,QG)
C
C If IFLAG=10, use the method : RK6(5)8S
C
      ELSEIF(IFLAG .EQ. 10) THEN
        CALL RK658S(M,COEFF,QG)
C
C If IFLAG=11, use the method : RK6(5)8C
C
      ELSEIF(IFLAG .EQ. 11) THEN
        CALL RK658C(M,COEFF,QG)
C
C If IFLAG=12, use the method : RKV6(5)
C
      ELSEIF(IFLAG .EQ. 12) THEN
        CALL RKV65(M,COEFF,QG)
C
C If IFLAG=13, use the method : RKF6(5)A
C
      ELSEIF(IFLAG .EQ. 13) THEN
        CALL RKF65A(M,COEFF,QG)
C
C If IFLAG=14, use the method : RKF6(5)B
C
      ELSEIF(IFLAG .EQ. 14) THEN
        CALL RKF65B(M,COEFF,QG)
C
C If IFLAG=15, use the method : RKC6(5)
C
      ELSEIF(IFLAG .EQ. 15) THEN
        CALL RKC65(M,COEFF,QG)
C
C If IFLAG=16, use the method : RKV6(5)9A
C
      ELSEIF(IFLAG .EQ. 16) THEN
        CALL RKV65A(M,COEFF,QG)
C
C If IFLAG=17, use the method : RKV6(5)9B
C
      ELSEIF(IFLAG .EQ. 17) THEN
        CALL RKV65B(M,COEFF,QG)
C
C If IFLAG=18, use the method : RKV7(6)
C
      ELSEIF(IFLAG .EQ. 18) THEN
        CALL RKV76(M,COEFF,QG)
C
C If IFLAG=19, use the method : RK8(7)13M
C
      ELSEIF(IFLAG .EQ. 19) THEN
        CALL RK8713(M,COEFF,QG)
C
C If IFLAG=20, use the method : RKF8(7)
C
      ELSEIF(IFLAG .EQ. 20) THEN
        CALL RKF87(M,COEFF,QG)
C
C If IFLAG=21, use the method : RKV8(7)
C
      ELSEIF(IFLAG .EQ. 21) THEN
        CALL RKV87(M,COEFF,QG)
C
C If IFLAG=22, use the method : RKV9(8)
C
      ELSEIF(IFLAG .EQ. 22) THEN
        CALL RKV98(M,COEFF,QG)
      ENDIF
      RETURN
      END
