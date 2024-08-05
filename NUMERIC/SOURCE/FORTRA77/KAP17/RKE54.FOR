C[BA*)
C[LE*)
      SUBROUTINE RKE54 (M,COEFF,QG)
C[IX{RKE54}*)
C
C*****************************************************************
C                                                                *
C This subroutine determines the coefficients for the embedding  *
C formula RKE5(4) (method of England).                           *
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
C  Sources  : G. Engeln-MÅllges F. Reuter :                      *
C             Formelsammlung zur Numerischen Mathematik mit      *
C             Standard FORTRAN 77 - Programmen                   *
C             page 334                                           *
C             Wissenschaftsverlag  1988                          *
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
C         a - values (see scheme )
C
      COEFF(2,1)=1.0D0/2.0D0
      COEFF(3,1)=1.0D0/2.0D0
      COEFF(4,1)=1.0D0
      COEFF(5,1)=2.0D0/3.0D0
      COEFF(6,1)=1.0D0/5.0D0
C
C         b - values (see scheme)
C
      COEFF(2,2)=1.0D0/2.0D0
      COEFF(3,2)=1.0D0/4.0D0
      COEFF(5,2)=7.0D0/27.0D0
      COEFF(6,2)=28.0D0/625.0D0
      COEFF(3,3)=1.0D0/4.0D0
      COEFF(4,3)=-1.0D0
      COEFF(5,3)=10.0D0/27.0D0
      COEFF(6,3)=-125.0D0/625.0D0
      COEFF(4,4)=2.0D0
      COEFF(6,4)=546.0D0/625.0D0
      COEFF(5,5)=1.0D0/27.0D0
      COEFF(6,5)=54.0D0/625.0D0
      COEFF(6,6)=-378.0D0/625.0D0
C
C         A tilde values (see scheme)
C
      COEFF(2,3)=14.0D0/336.0D0
      COEFF(2,6)=35.0D0/336.0D0
      COEFF(3,4)=162.0D0/336.0D0
      COEFF(3,5)=125.0D0/336.0D0
C
C          A - values (see scheme)
C
      COEFF(1,1)=1.0D0/6.0D0
      COEFF(1,3)=4.0D0/6.0D0
      COEFF(1,4)=1.0D0/6.0D0
      RETURN
      END
