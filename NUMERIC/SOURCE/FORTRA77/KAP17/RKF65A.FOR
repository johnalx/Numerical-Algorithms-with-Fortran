C[BA*)
C[LE*)
      SUBROUTINE RKF65A (M,COEFF,QG)
C[IX{RKF65A}*)
C
C*****************************************************************
C                                                                *
C This subroutine determines the coefficients for the embedding  *
C formula RKF6(5)A (method of Fehlberg).                          *
C[BE*)
C                                                                *
C                                                                *
C INPUT PARAMETERS:                                              *
C =================                                              *
C M       : Dimension of the matrix COEFF depending on the chosen*
C           embedding formula                                    *
C                                                                *
C                                                                *
C OUTPUT PARAMETERS:                                             *
C ==================                                             *
C COEFF   : 2-dimensional DOUBLE PRECISION array COEFF(1:16,1:M) *
C           with the coefficients for the embedding formula      *
C QG      : DOUBLE PRECISION value for the global error order    *
C                                                                *
C                                                                *
C----------------------------------------------------------------*
C                                                                *
C  required subroutines: none                                    *
C                                                                *
C*****************************************************************
C                                                                *
C  Author   : Volker KrÅger                                      *
C  Date     : 26.04.1993                                         *
C  Source   : FORTRAN 77                                         *
C  Sources  : E. Fehlberg :                                      *
C             Klassische Runge-Kutta-Formeln fÅnfter und         *
C             siebenter Ordnung mit Schrittweiten-Kontrolle      *
C             page 93-106                                        *
C             Computing 4 1969                                   *
C                                                                *
C[BA*)
C*****************************************************************
C[BE*)
C
C Declarations
C
      DOUBLE PRECISION COEFF(16,M),QG
C
C Initialize QG
C
      QG=5.0D0
C
C Determine the matrix elements in COEFF
C
C         a - values (see scheme)
C
      COEFF(2,1)=1.0D0/6.0D0
      COEFF(3,1)=4.0D0/15.0D0
      COEFF(4,1)=2.0D0/3.0D0
      COEFF(5,1)=4.0D0/5.0D0
      COEFF(6,1)=1.0D0
      COEFF(8,1)=1.0D0
C
C         b - values (see scheme)
C
      COEFF(2,2)=1.0D0/6.0D0
      COEFF(3,2)=4.0D0/75.0D0
      COEFF(4,2)=5.0D0/6.0D0
      COEFF(5,2)=-8.0D0/5.0D0
      COEFF(6,2)=361.0D0/320.0D0
      COEFF(7,2)=-11.0D0/640.0D0
      COEFF(8,2)=93.0D0/640.0D0
      COEFF(3,3)=16.0D0/75.0D0
      COEFF(4,3)=-8.0D0/3.0D0
      COEFF(5,3)=144.0D0/25.0D0
      COEFF(6,3)=-18.0D0/5.0D0
      COEFF(8,3)=-18.0D0/5.0D0
      COEFF(4,4)=5.0D0/2.0D0
      COEFF(5,4)=-4.0D0
      COEFF(6,4)=407.0D0/128.0D0
      COEFF(7,4)=11.0D0/256.0D0
      COEFF(8,4)=803.0D0/256.0D0
      COEFF(5,5)=16.0D0/25.0D0
      COEFF(6,5)=-11.0D0/80.0D0
      COEFF(7,5)=-11.0D0/160.0D0
      COEFF(8,5)=-11.0D0/160.0D0
      COEFF(6,6)=55.0D0/128.0D0
      COEFF(7,6)=11.0D0/256.0D0
      COEFF(8,6)=99.0D0/256.0D0
      COEFF(8,8)=1.0D0
C
C         A tilde values (see scheme)
C
      COEFF(2,3)=7.0D0/1408.0D0
      COEFF(2,5)=1125.0D0/2816.0D0
      COEFF(2,6)=9.0D0/32.0D0
      COEFF(2,7)=125.0D0/768.0D0
      COEFF(3,4)=5.0D0/66.0D0
      COEFF(3,5)=5.0D0/66.0D0
C
C          A - values (see scheme)
C
      COEFF(1,1)=31.0D0/384.0D0
      COEFF(1,3)=1125.0D0/2816.0D0
      COEFF(1,4)=9.0D0/32.0D0
      COEFF(1,5)=125.0D0/768.0D0
      COEFF(1,6)=5.0D0/66.0D0
      RETURN
      END
