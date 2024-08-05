C[BA*)
C[LE*)
      SUBROUTINE RKV65B (M,COEFF,QG)
C[IX{RKV65B}*)
C
C*****************************************************************
C                                                                *
C This subroutine determines the coefficients for the embedding  *
C formula RKV6(5)9B (method of Verner).                          *
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
C  Sources  : J. H. Verner :                                     *
C             Some Runge-Kutta formular pairs                    *
C             page 496-511                                       *
C             SIAM J. NUMER. ANAL., Vol. 28, No. 2, 1991         *
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
C         a - values (see scheme )
C
      COEFF(2,1)=1.0D0/8.0D0
      COEFF(3,1)=1.0D0/6.0D0
      COEFF(4,1)=1.0D0/4.0D0
      COEFF(5,1)=1.0D0/2.0D0
      COEFF(6,1)=3.0D0/5.0D0
      COEFF(7,1)=4.0D0/5.0D0
      COEFF(8,1)=1.0D0
      COEFF(9,1)=1.0D0
C
C         b - values (see scheme)
C
      COEFF(2,2)=1.0D0/8.0D0
      COEFF(3,2)=1.0D0/18.0D0
      COEFF(4,2)=1.0D0/16.0D0
      COEFF(5,2)=1.0D0/4.0D0
      COEFF(6,2)=134.0D0/625.0D0
      COEFF(7,2)=-98.0D0/1875.0D0
      COEFF(8,2)=9.0D0/50.0D0
      COEFF(9,2)=11.0D0/144.0D0
      COEFF(3,3)=1.0D0/9.0D0
      COEFF(4,4)=3.0D0/16.0D0
      COEFF(5,4)=-3.0D0/4.0D0
      COEFF(6,4)=-333.0D0/625.0D0
      COEFF(7,4)=12.0D0/625.0D0
      COEFF(8,4)=21.0D0/25.0D0
      COEFF(5,5)=1.0D0
      COEFF(6,5)=476.0D0/625.0D0
      COEFF(7,5)=107.36D2/131.25D2
      COEFF(8,5)=-2924.0D0/1925.0D0
      COEFF(9,5)=256.0D0/693.0D0
      COEFF(6,6)=98.0D0/625.0D0
      COEFF(7,6)=-1936.0D0/1875.0D0
      COEFF(8,6)=74.0D0/25.0D0
      COEFF(7,7)=22.0D0/21.0D0
      COEFF(8,7)=-15.0D0/7.0D0
      COEFF(9,7)=125.0D0/504.0D0
      COEFF(8,8)=15.0D0/22.0D0
      COEFF(9,8)=125.0D0/528.0D0
      COEFF(9,9)=5.0D0/72.0D0
C
C         A tilde values (see scheme)
C
      COEFF(2,3)=11.0D0/144.0D0
      COEFF(2,6)=256.0D0/693.0D0
      COEFF(2,8)=125.0D0/504.0D0
      COEFF(2,9)=125.0D0/528.0D0
      COEFF(3,4)=5.0D0/72.0D0
C
C          A - values (see scheme)
C
      COEFF(1,1)=1.0D0/18.0D0
      COEFF(1,4)=32.0D0/63.0D0
      COEFF(1,5)=-2.0D0/3.0D0
      COEFF(1,6)=125.0D0/126.0D0
      COEFF(1,8)=-5.0D0/63.0D0
      COEFF(1,9)=4.0D0/21.0D0
      RETURN
      END
