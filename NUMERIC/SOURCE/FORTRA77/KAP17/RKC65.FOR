C[BA*)
C[LE*)
      SUBROUTINE RKC65 (M,COEFF,QG)
C[IX{RKC65}*)
C
C*****************************************************************
C                                                                *
C This subroutine determines the coefficients for the embedding  *
C formula RKC6(5) (method of Calvo).                             *
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
C  Sources  : M. Calvo, J. I. Montijano, L. Randez :             *
C             A new embedded pair of Runge-Kutta formulas of     *
C             orders 5 and 6                                     *
C             page 15-24                                         *
C             Computers Math. Applic., Vol. 20, No. 1, 1990      *
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
      COEFF(2,1)=2.0D0/15.0D0
      COEFF(3,1)=1.0D0/5.0D0
      COEFF(4,1)=3.0D0/10.0D0
      COEFF(5,1)=14.0D0/25.0D0
      COEFF(6,1)=19.0D0/25.0D0
      COEFF(7,1)=3.5226607D7/3.5688279D7
      COEFF(8,1)=1.0D0
      COEFF(9,1)=1.0D0
C
C         b - values (see scheme)
C
      COEFF(2,2)=2.0D0/15.0D0
      COEFF(3,2)=1.0D0/20.0D0
      COEFF(4,2)=3.0D0/40.0D0
      COEFF(5,2)=8.6727015D7/19.6851553D7
      COEFF(6,2)=-8.6860849D7/4.5628967D7
      COEFF(7,2)=7.7759591D7/1.6096467D7
      COEFF(8,2)=23.7564263D7/3.9280295D7
      COEFF(9,2)=1.7572349D7/28.9262523D7
      COEFF(3,3)=3.0D0/20.0D0
      COEFF(5,3)=-6.0129073D7/5.2624712D7
      COEFF(6,3)=11.1022885D7/2.5716487D7
      COEFF(7,3)=-4.9252809D7/6.452555D6
      COEFF(8,3)=-10.0523239D7/1.0677940D7
      COEFF(4,4)=9.0D0/40.0D0
      COEFF(5,4)=95.7436434D7/137.8352377D7
      COEFF(6,4)=10.8046682D7/10.1167669D7
      COEFF(7,4)=-38.1680111D7/5.1572984D7
      COEFF(8,4)=-26.5574846D7/2.7330247D7
      COEFF(9,4)=5.7513011D7/20.1864250D7
      COEFF(5,5)=8.3886832D7/14.7842441D7
      COEFF(6,5)=-14.1756746D7/3.6005461D7
      COEFF(7,5)=87.9269579D7/6.6788831D7
      COEFF(8,5)=31.7978411D7/1.8988713D7
      COEFF(9,5)=1.5587306D7/35.4501571D7
      COEFF(6,6)=7.3139862D7/6.0170633D7
      COEFF(7,6)=-9.0453121D7/3.3722162D7
      COEFF(8,6)=-12.4494385D7/3.5453627D7
      COEFF(9,6)=7.1783021D7/23.4982865D7
      COEFF(7,7)=11.1179552D7/15.7155827D7
      COEFF(8,7)=8.6822444D7/10.0138635D7
      COEFF(9,7)=2.9672000D7/18.0480167D7
      COEFF(8,8)=-1.2873523D7/72.4232625D7
      COEFF(9,8)=6.5567621D7/12.7060952D7
      COEFF(9,9)=-7.9074570D7/21.0557597D7
C
C         A tilde values (see scheme)
C
      COEFF(2,3)=1.7572349D7/28.9262523D7
      COEFF(2,5)=5.7513011D7/20.1864250D7
      COEFF(2,6)=1.5587306D7/35.4501571D7
      COEFF(2,7)=7.1783021D7/23.4982865D7
      COEFF(2,8)=2.9672000D7/18.0480167D7
      COEFF(2,9)=6.5567621D7/12.7060952D7
      COEFF(3,4)=-7.9074570D7/21.0557597D7
C
C          A - values (see scheme)
C
      COEFF(1,1)=1.5231665D7/51.0830334D7
      COEFF(1,3)=5.9452991D7/11.6050448D7
      COEFF(1,4)=-2.8398517D7/12.2437738D7
      COEFF(1,5)=5.6673824D7/13.7010559D7
      COEFF(1,6)=6.8003849D7/42.6673583D7
      COEFF(1,7)=7.097631D6/3.7564021D7
      COEFF(1,8)=-7.1226429D7/58.3093742D7
      COEFF(1,9)=1.0D0/20.0D0
      RETURN
      END
