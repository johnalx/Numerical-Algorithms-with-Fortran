C[BA*)
C[LE*)
      SUBROUTINE RK658C (M,COEFF,QG)
C[IX{RK658C}*)
C
C*****************************************************************
C                                                                *
C This subroutine determines the coefficients for the embedding  *
C formula RK6(5)8C (method of Dormand und Prince).               *
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
C  Sources  : J. R. Dormand P. J. Prince :                       *
C             A reconsideration of some embedded Runge Kutta     *
C             formulae                                           *
C             page 203-211                                       *
C             Journal of Computational and Applied Mathematics,  *
C             15, 1986, North-Holland                            *
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
      COEFF(2,1)=1.0D0/10.0D0
      COEFF(3,1)=1.0D0/6.0D0
      COEFF(4,1)=2.0D0/9.0D0
      COEFF(5,1)=3.0D0/5.0D0
      COEFF(6,1)=4.0D0/5.0D0
      COEFF(7,1)=1.0D0
      COEFF(8,1)=1.0D0
C
C         b - values (see scheme)
C
      COEFF(2,2)=1.0D0/10.0D0
      COEFF(3,2)=1.0D0/36.0D0
      COEFF(4,2)=10.0D0/243.0D0
      COEFF(5,2)=4047.0D0/5500.0D0
      COEFF(6,2)=-5587.0D0/4125.0D0
      COEFF(7,2)=12961.0D0/2376.0D0
      COEFF(8,2)=702799.0D0/199584.0D0
      COEFF(3,3)=5.0D0/36.0D0
      COEFF(4,3)=20.0D0/243.0D0
      COEFF(5,3)=-18.0D0/55.0D0
      COEFF(6,3)=24.0D0/55.0D0
      COEFF(7,3)=-35.0D0/33.0D0
      COEFF(8,3)=-1865.0D0/2772.0D0
      COEFF(4,4)=8.0D0/81.0D0
      COEFF(5,4)=-4212.0D0/1375.0D0
      COEFF(6,4)=9576.0D0/1375.0D0
      COEFF(7,4)=-160845.0D0/5434.0D0
      COEFF(8,4)=-2891375.0D0/152152.0D0
      COEFF(5,5)=17901.0D0/5500.0D0
      COEFF(6,5)=-140049.0D0/23375.0D0
      COEFF(7,5)=1067565.0D0/38896.0D0
      COEFF(8,5)=19332955.0D0/1089088.0D0
      COEFF(6,6)=38.0D0/51.0D0
      COEFF(7,6)=-103375.0D0/47736.0D0
      COEFF(8,6)=-5356375.0D0/4009824.0D0
      COEFF(7,7)=32875.0D0/35568.0D0
      COEFF(8,7)=2207875.0D0/2987712.0D0
C
C         A tilde values (see scheme)
C
      COEFF(2,3)=1.0D0/12.0D0
      COEFF(2,5)=-216.0D0/1235.0D0
      COEFF(2,6)=6561.0D0/12376.0D0
      COEFF(2,7)=1375.0D0/5304.0D0
      COEFF(2,8)=1375.0D0/5928.0D0
      COEFF(3,4)=-5.0D0/168.0D0
      COEFF(3,5)=1.0D0/10.0D0
C
C          A - values (see scheme)
C
      COEFF(1,1)=163.0D0/1440.0D0
      COEFF(1,3)=-2628.0D0/6175.0D0
      COEFF(1,4)=13851.0D0/17680.0D0
      COEFF(1,5)=1525.0D0/7956.0D0
      COEFF(1,6)=6575.0D0/23712.0D0
      COEFF(1,7)=3.0D0/50.0D0
      RETURN
      END
