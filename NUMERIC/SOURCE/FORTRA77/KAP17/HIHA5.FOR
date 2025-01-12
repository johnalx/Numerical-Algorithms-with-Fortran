C[BA*)
C[LE*)
      SUBROUTINE HIHA5 (M,COEFF,QG)
C[IX{HIHA5}*)
C
C*****************************************************************
C                                                                *
C This subroutine determines the coefficients for the embedding  *
C formula HIHA5 (method of Higham und Hall).                     *
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
C  Author   : Volker Kr�ger                                      *
C  Date     : 26.04.1993                                         *
C  Source   : FORTRAN 77                                         *
C  sources  : E. Hairer G. Wanner :                              *
C             Solving Ordinary Differential Equations II         *
C             page 28-31                                         *
C             Springer-Verlag Berlin Heidelberg 1991             *
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
      QG=4.0D0
C
C Determine the matrix elements in COEFF
C
C         a - values (see scheme)
C
      COEFF(2,1)=2.0D0/9.0D0
      COEFF(3,1)=1.0D0/3.0D0
      COEFF(4,1)=1.0D0/2.0D0
      COEFF(5,1)=3.0D0/5.0D0
      COEFF(6,1)=1.0D0
      COEFF(7,1)=1.0D0
C
C         b - values (see scheme)
C
      COEFF(2,2)=2.0D0/9.0D0
      COEFF(3,2)=1.0D0/12.0D0
      COEFF(4,2)=1.0D0/8.0D0
      COEFF(5,2)=91.0D0/500.0D0
      COEFF(6,2)=-11.0D0/20.0D0
      COEFF(7,2)=1.0D0/12.0D0
      COEFF(3,3)=1.0D0/4.0D0
      COEFF(5,3)=-27.0D0/100.0D0
      COEFF(6,3)=27.0D0/20.0D0
      COEFF(4,4)=3.0D0/8.0D0
      COEFF(5,4)=78.0D0/125.0D0
      COEFF(6,4)=12.0D0/5.0D0
      COEFF(7,4)=27.0D0/32.0D0
      COEFF(5,5)=8.0D0/125.0D0
      COEFF(6,5)=-36.0D0/5.0D0
      COEFF(7,5)=-4.0D0/3.0D0
      COEFF(6,6)=5.0D0
      COEFF(7,6)=125.0D0/96.0D0
      COEFF(7,7)=5.0D0/48.0D0
C
C         A tilde values (see scheme)
C
      COEFF(2,3)=1.0D0/12.0D0
      COEFF(2,5)=27.0D0/32.0D0
      COEFF(2,6)=-4.0D0/3.0D0
      COEFF(2,7)=125.0D0/96.0D0
      COEFF(3,4)=5.0D0/48.0D0
C
C          A - values (see scheme)
C

      COEFF(1,1)=2.0D0/15.0D0
      COEFF(1,3)=27.0D0/80.0D0
      COEFF(1,4)=-2.0D0/15.0D0
      COEFF(1,5)=25.0D0/48.0D0
      COEFF(1,6)=1.0D0/24.0D0
      COEFF(1,7)=1.0D0/10.0D0
      RETURN
      END
