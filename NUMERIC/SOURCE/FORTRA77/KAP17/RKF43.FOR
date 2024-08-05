C[BA*)
C[LE*)
      SUBROUTINE RKF43 (M,COEFF,QG)
C[IX{RKF43}*)
C
C*****************************************************************
C                                                                *
C This subroutine determines the coefficients for the embedding  *
C formula RKF4(3) (method of Fehlberg).                          *
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
C             Klassische Runge-Kutta-Formeln vierter und         *
C             niedrigerer Ordnung mit Schrittweiten-Kontrolle    *
C             und ihrer Anwendung auf WÑrmeleitungsprobleme      *
C             page 61-71                                         *
C             Computing 6 1970                                   *
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
      QG=3.0D0
C
C Determine the matrix elements in COEFF
C
C         a - values (refer to the coefficient schemes in chapter 17.3.4.4)
C                    ( similarly below )
C
      COEFF(2,1)=2.0D0/7.0D0
      COEFF(3,1)=7.0D0/15.0D0
      COEFF(4,1)=35.0D0/38.0D0
      COEFF(5,1)=1.0D0
C
C         b - values (see scheme)
C
      COEFF(2,2)=2.0D0/7.0D0
      COEFF(3,2)=77.0D0/900.0D0
      COEFF(4,2)=805.0D0/1444.0D0
      COEFF(5,2)=79.0D0/490.0D0
      COEFF(3,3)=343.0D0/900.0D0
      COEFF(4,3)=-77175.0D0/54872.0D0
      COEFF(4,4)=97125.0D0/54872.0D0
      COEFF(5,4)=2175.0D0/3626.0D0
      COEFF(5,5)=2166.0D0/9065.0D0
C
C         A tilde values (see scheme)
C
      COEFF(2,3)=79.0D0/490.0D0
      COEFF(2,5)=2175.0D0/3626.0D0
      COEFF(3,4)=2166.0D0/9065.0D0
C
C          A - values (see scheme)
C
      COEFF(1,1)=229.0D0/1470.0D0
      COEFF(1,3)=1125.0D0/1813.0D0
      COEFF(1,4)=13718.0D0/81585.0D0
      COEFF(1,5)=1.0D0/18.0D0
      RETURN
      END
